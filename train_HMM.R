# Clear the console
cat("\014")

# Vérifier et installer les packages manquants
required_packages <- c("depmixS4", "caret", "pROC", "openxlsx")

# Fonction pour installer les packages manquants
install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# Installer les packages
install_if_missing(required_packages)


# Charger les bibliothèques nécessaires
library(depmixS4)
library(caret)
library(dplyr)

# Charger les données
data <- readxl::read_xlsx("data/output/module_rankings_with_risk.xlsx")

prepare_data <- function(data) {
  # Pivot les colonnes Rank_1, Rank_2, ..., en une seule colonne nommée "module"
  data_long <- data %>%
    pivot_longer(
      cols = starts_with("Rank"),
      names_to = "Rank",
      values_to = "module"
    )
  
  print("Colonnes après pivot_longer :")
  print(colnames(data_long))  # Vérifiez que Sample, module, et Risk sont présents
  
  print("Aperçu des données après pivot_longer :")
  print(head(data_long))
  
  data_long <- data_long %>%
    dplyr::select(Sample, module, Risk)
  
  print("Aperçu final des données prêtes :")
  print(head(data_long))
  
  return(data_long)
}


data_long <- prepare_data(data)

# Diviser les données en ensembles d'entraînement et de test
set.seed(123)  # Pour la reproductibilité
train_idx <- createDataPartition(data_long$Risk, p = 0.8, list = FALSE)
train_data <- data_long[train_idx, ]
test_data <- data_long[-train_idx, ]

# Fonction pour entraîner le modèle HMM
train_hmm <- function(train_data, N) {
  # Début du compteur pour l'apprentissage
  start_time <- Sys.time()
  
  # Création et entraînement du modèle HMM
  hmm_model <- depmix(response = module ~ 1, data = train_data, nstates = N, family = multinomial())
  hmm_model <- fit(hmm_model)
  
  # Fin du compteur pour l'apprentissage
  end_time <- Sys.time()
  train_duration <- end_time - start_time
  cat("Temps d'apprentissage :", train_duration, "\n")
  
  return(hmm_model)
}

# Fonction de validation croisée
cross_validate <- function(data, N, k = 5) {
  folds <- createFolds(data$Risk, k = k, list = TRUE, returnTrain = TRUE)
  results <- lapply(folds, function(train_idx) {
    # Données d'entraînement pour ce fold
    train_fold <- data[train_idx, ]
    train_fold$module <- as.factor(train_fold$module)
    
    # Entraîner le modèle
    hmm_fit <- train_hmm(train_fold, N)
    
    # Données de validation pour ce fold
    val_fold <- data[-train_idx, ]
    val_fold$module <- as.factor(val_fold$module)
    
    # Calculer la log-vraisemblance sur les données de validation
    log_likelihood <- forward_algorithm(hmm_fit, val_fold)
    return(log_likelihood)
  })
  return(results)
}

# Fonction pour calculer la probabilité avec l'algorithme forward
forward_algorithm <- function(model, data) {
  seq_data <- data$module
  seq_data <- as.factor(seq_data)  # S'assurer que c'est un facteur
  forward_probs <- forwardbackward(model)$logLik
  return(forward_probs)
}

# Effectuer la validation croisée sur les données d'entraînement
N <- 2  # Nombre d'états ajustable
cv_results <- cross_validate(train_data, N, k = 2)

# Évaluation sur les données de test
evaluate_hmm <- function(train_data, test_data, N) {
  # Début du compteur pour la prédiction
  start_time <- Sys.time()
  
  # Entraîner les HMM sur chaque catégorie
  hmm_high <- train_hmm(train_data %>% filter(Risk == "High"), N)
  hmm_low <- train_hmm(train_data %>% filter(Risk == "Low"), N)
  
  # Calcul des probabilités pour les données de test
  test_data$high_prob <- sapply(split(test_data, test_data$Sample), 
                                function(sample) forward_algorithm(hmm_high, sample))
  test_data$low_prob <- sapply(split(test_data, test_data$Sample), 
                               function(sample) forward_algorithm(hmm_low, sample))
  
  # Prédiction basée sur les probabilités maximales
  test_data$predicted_risk <- ifelse(test_data$high_prob > test_data$low_prob, "High", "Low")
  
  # Fin du compteur pour la prédiction
  end_time <- Sys.time()
  prediction_duration <- end_time - start_time
  cat("Temps de prédiction :", prediction_duration, "\n")
  
  # Calcul des métriques de performance
  confusion_matrix <- caret::confusionMatrix(as.factor(test_data$predicted_risk), 
                                             as.factor(test_data$Risk))
  return(list(confusion_matrix = confusion_matrix, test_data = test_data))
}

# Évaluer sur l'ensemble de test
evaluation_results <- evaluate_hmm(train_data, test_data, N)

# Résultats
print(evaluation_results$confusion_matrix)
