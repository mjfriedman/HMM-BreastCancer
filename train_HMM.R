# Nettoyer l'environnement et la console
rm(list = ls())
cat("\014")

# Vérifier et installer les packages manquants
required_packages <- c("HMM", "caret", "pROC", "Metrics", 
                       "tidyverse", "readxl", "openxlsx", "timeDate", "dplyr")

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

install_if_missing(required_packages)

# Charger les bibliothèques nécessaires
lapply(required_packages, library, character.only = TRUE)

# Charger les données
data <- read.xlsx("data/output/module_rankings_with_risk.xlsx", sheet = 1)

# Vérification des données
if (is.null(data) || nrow(data) == 0) {
  stop("Les données sont vides ou non valides.")
}

# Séparation des données en deux groupes : train et test
set.seed(123)
train_index <- sample(1:nrow(data), size = 0.8 * nrow(data))
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

# Fonction pour préparer les données pour HMM
prepare_data_for_hmm <- function(data) {
  return(as.character(as.vector(as.matrix(data))))  # Convertir en vecteur de caractères
}

# Préparer les données d'entraînement et de test
train_high <- train_data %>% filter(Risk == "High") %>% dplyr::select(starts_with("Rank_"))
train_high_seq <- prepare_data_for_hmm(train_high)

test_high <- test_data %>% filter(Risk == "High") %>% dplyr::select(starts_with("Rank_"))
test_high_seq <- prepare_data_for_hmm(test_high)

train_low <- train_data %>% filter(Risk == "Low") %>% dplyr::select(starts_with("Rank_"))
train_low_seq <- prepare_data_for_hmm(train_low)

test_low <- test_data %>% filter(Risk == "Low") %>% dplyr::select(starts_with("Rank_"))
test_low_seq <- prepare_data_for_hmm(test_low)

# Identifier les modules uniques
modules <- unique(c(unlist(train_high_seq), unlist(train_low_seq)))
n_modules <- length(modules)

# Initialisation des matrices pour HMM
initialize_hmm <- function(train_data, modules) {
  # Transition gauche-droite
  transition_matrix <- matrix(runif(n_modules * n_modules), nrow = n_modules, ncol = n_modules)
  #transition_matrix <- matrix(0, nrow = n_modules, ncol = n_modules)
  #for (i in 1:(n_modules - 1)) {
  #  transition_matrix[i, i + 1] <- 1
  #}
  #transition_matrix[n_modules, n_modules] <- 1
  # 
  # # Probabilités initiales
  #initial_state_probabilities <- rep(0, n_modules)
  #initial_state_probabilities[1] <- 1
  # 
  # # Matrice d'émission spécifique à l'ensemble d'entraînement
  # emission_matrix <- matrix(0, nrow = n_modules, ncol = n_modules)
  # for (i in 1:n_modules) {
  #   for (j in 1:n_modules) {
  #     module <- modules[j]
  #     occurrence_count <- sum(train_data[, i] == module, na.rm = TRUE)
  #     emission_matrix[i, j] <- occurrence_count / nrow(train_data)
  #   }
  # }
  # # Normalisation
  # #emission_matrix <- t(apply(emission_matrix, 1, function(row) row / sum(row, na.rm = TRUE)))
  # #emission_matrix[is.na(emission_matrix)] <- 1 / n_modules
  ########
  initial_state_probabilities <- runif(n_modules)
  emission_matrix <- matrix(runif(n_modules * length(modules)), 
                           nrow = n_modules, ncol = length(modules))
  
  # Normalize rows to ensure probabilities sum to 1
  initial_state_probabilities <- initial_state_probabilities / sum(initial_state_probabilities)
  transition_matrix <- transition_matrix / rowSums(transition_matrix)
  emission_matrix <- emission_matrix / rowSums(emission_matrix)
  ##########
  # Créer le modèle HMM
  initHMM(
    States = modules,
    Symbols = modules,
    startProbs = initial_state_probabilities,
    transProbs = transition_matrix,
    emissionProbs = emission_matrix
  )
}


# Trainer
trainer <- function(hmm, data, maxIterations = 300) {
  log_likelihoods <- numeric()
  
  start_time <- Sys.time()
  
  # Vérifier les dimensions des matrices et des données
  cat("Transition matrix dimensions: ", dim(hmm$transProbs), "\n")
  cat("Emission matrix dimensions: ", dim(hmm$emissionProbs), "\n")
  cat("Number of training sequences: ", length(data), "\n")
  
  # Baum-Welch sur les données d'entraînement
  bw <- baumWelch(hmm, data, maxIterations = maxIterations)
  
  # Extraire les matrices nécessaires pour l'algorithme Forward
  A <- bw$hmm$transProbs    # Matrice de transition
  B <- bw$hmm$emissionProbs # Matrice d'émission
  pi <- bw$hmm$startProbs   # Probabilités initiales
  
  # Récupérer les log-vraisemblances
  log_likelihoods <- bw$difference
  end_time <- Sys.time()
  
  # Retourner les matrices et les log-vraisemblances
  return(list(A = A,  B = B, pi = pi, log_likelihoods = log_likelihoods
  ))
}


# Initialiser et entraîner le modèle pour "High Risk"
cat("\n\nTraining HMM : High risk")
hmm_high <- initialize_hmm(train_high, modules)
train_start_time <- Sys.time()
bw_high <- trainer(hmm_high, train_high_seq, maxIterations = 300)
# Accéder aux matrices nécessaires pour l'algorithme Forward
A_high <- bw_high$A
B_high <- bw_high$B
pi_high <- bw_high$pi
log_likelihoods_high <- bw_high$log_likelihoods
train_end_time <- Sys.time()
delta_time <- as.numeric(difftime(train_end_time, train_start_time, units = "secs"))
cat("Temps total d'entrainement:", delta_time, "seconds\n")

# Initialiser et entraîner le modèle pour "Low Risk"
cat("\n\nTraining HMM : Low risk")
hmm_low <- initialize_hmm(train_low, modules)
train_start_time <- Sys.time()
bw_low <- trainer(hmm_low, train_low_seq, maxIterations = 700)
# Accéder aux matrices
A_low <- bw_low$A
B_low <- bw_low$B
pi_low <- bw_low$pi
log_likelihoods_low <- bw_low$log_likelihoods
train_end_time <- Sys.time()
delta_time <- as.numeric(difftime(train_end_time, train_start_time, units = "secs"))
cat("Temps total d'entrainement:", delta_time, "seconds\n")

# Plot des log-vraisemblances
plot(log_likelihoods_high, type = "o", col = "blue", xlab = "Itération", ylab = "Log vraisemblance", main = "Log vraisemblance par Itération")
lines(log_likelihoods_low, type = "o", col = "red")
legend("topright", legend = c("Haut Risque", "Faible Risque"), col = c("blue", "red"), lty = 1)

# Fonction pour l'algorithme Forward
forward_algorithm <- function(sequence, A, B, pi) {
  n_states <- nrow(A)  # Nombre d'états
  T <- length(sequence)  # Longueur de la séquence d'observations
  
  # Initialiser la matrice alpha (Forward probabilities)
  alpha <- matrix(0, nrow = n_states, ncol = T)
  
  # Initialisation au temps t=1
  for (i in 1:n_states) {
    alpha[i, 1] <- pi[i] * B[i, sequence[1]]
  }
  
  # Normalisation (optionnelle)
  #alpha[, 1] <- alpha[, 1] / sum(alpha[, 1])
  
  # Récurrence pour t = 2 à T
  for (t in 2:T) {
    for (j in 1:n_states) {
      alpha[j, t] <- sum(alpha[, t - 1] * A[, j]) * B[j, sequence[t]]
    }
    # Normalisation (optionnelle)
    #alpha[, t] <- alpha[, t] / sum(alpha[, t])
  }
  
  # Probabilité totale de la séquence
  prob_sequence <- sum(alpha[, T])
  return(prob_sequence)
}


############################################
library(caret)

# Fonction d'évaluation
evaluate_test_samples <- function(test_data, A_high, B_high, pi_high, A_low, B_low, pi_low) {
  predictions <- c()
  
  true_labels <- test_data$Risk  # Les étiquettes réelles

  # Identifier les colonnes correspondant aux séquences
  rank_columns <- colnames(test_data)[grepl("^Rank_", colnames(test_data))]

  for (i in 1:nrow(test_data)) {
    # Extraire une séquence
    cat("sample =", i, "\n")
    sample <- as.character(unlist(test_data[i, rank_columns]))  # Extraire comme vecteur de caractères
    cat("sample : ", sample, "\n")
    
    # Calculer les log-vraisemblances pour les deux modèles
    log_prob_high <- forward_algorithm(sample, A_high, B_high, pi_high)
    log_prob_low <- forward_algorithm(sample, A_low, B_low, pi_low)
    cat(" log_prob high =", log_prob_high, "and log_prob low=", log_prob_low, "\n")
    
    # Prédiction basée sur les log-vraisemblances
    predicted_label <- ifelse(log_prob_high > log_prob_low, "High", "Low")
    predictions <- c(predictions, predicted_label)
  }
  
  return(predictions)
}

# Appliquer l'évaluation
cat("\nÉvaluation des données de test...\n")
predictions <- evaluate_test_samples(test_data, A_high, B_high, pi_high, A_low, B_low, pi_low)

# Créer une matrice de confusion
confusion_matrix <- table(Predicted = predictions, Actual = test_data$Risk)
print("Matrice de confusion :")
print(confusion_matrix)

# Calcul des métriques : Sensibilité, Spécificité, MCC, F1-Score
true_positive <- confusion_matrix["High", "High"]
false_negative <- confusion_matrix["Low", "High"]
true_negative <- confusion_matrix["Low", "Low"]
false_positive <- confusion_matrix["High", "Low"]

# Sensibilité (TPR)
sensitivity <- ifelse((true_positive + false_negative) > 0, 
                      true_positive / (true_positive + false_negative), 
                      NA)

# Spécificité (TNR)
specificity <- ifelse((true_negative + false_positive) > 0, 
                      true_negative / (true_negative + false_positive), 
                      NA)

# MCC (Matthew's Correlation Coefficient)
numerator <- (true_positive * true_negative) - (false_positive * false_negative)
denominator <- sqrt((true_positive + false_positive) * 
                      (true_positive + false_negative) * 
                      (true_negative + false_positive) * 
                      (true_negative + false_negative))

mcc <- ifelse(denominator > 0, numerator / denominator, NA)

# F1-Score
f1 <- F1_Score(predictions, test_data$Risk, positive = "High")

# Afficher les résultats
cat("\nRésultats de l'évaluation :\n")
cat("Sensibilité (TPR) :", round(sensitivity, 4), "\n")
cat("Spécificité (TNR) :", round(specificity, 4), "\n")
cat("MCC :", round(mcc, 4), "\n")
cat("F1-Score :", round(f1, 4), "\n")
