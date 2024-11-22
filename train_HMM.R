# Nettoyer l'environnement et la console
rm(list = ls())
cat("\014")

# Vérifier et installer les packages manquants
required_packages <- c("depmixS4", "caret", "pROC", "Metrics", 
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

# Charger les bibliothèques
library(depmixS4)
library(caret)
library(pROC)
library(Metrics)
library(tidyverse)
library(readxl)
library(openxlsx)
library(timeDate)
library(readxl)
library(dplyr)

# Charger les bibliothèques nécessaires


# Charger le jeu de données (adapte le chemin en fonction de ton environnement)
data <- read.xlsx("data/output/module_rankings_with_risk.xlsx", sheet = 1)

# Séparation des données en deux groupes : train et test
set.seed(123)  # Pour reproductibilité
train_index <- sample(1:nrow(data), size = 0.8 * nrow(data))  # 80% pour l'entraînement
train_data <- data[train_index, ]
test_data <- data[-train_index, ]

# Extraire les observations "Rank" pour l'ensemble "High"
train_high <- train_data %>% filter(Risk == "High")
train_high <- train_high[, grep("^Rank_", colnames(train_high))]

test_high <- test_data %>% filter(Risk == "High")
test_high <- test_high[, grep("^Rank_", colnames(test_high))]


# Préparer les éléments de base pour le modèle HMM
modules <- unique(unlist(train_high))  # Modules uniques dans les données d'entraînement
n_modules <- length(modules)  # Nombre de modules distincts

# Définir le nombre d'états (qui correspond au nombre de modules)
n_states <- n_modules

# Initialisation de la matrice de transition (modèle gauche-droite)
transition_matrix <- matrix(0, nrow = n_states, ncol = n_states)
for (i in 1:(n_states - 1)) {
  transition_matrix[i, i + 1] <- 1  # Transition de i vers i+1
}

# Initialisation de la distribution des états initiaux
initial_state_probabilities <- rep(0, n_states)
initial_state_probabilities[1] <- 1  # 1 pour le premier état

# Calcul de la matrice d'émission (probabilité d'émettre un module dans un état)
emission_matrix <- matrix(0, nrow = n_states, ncol = n_modules)
for (i in 1:n_states) {
  for (j in 1:n_modules) {
    module <- modules[j]
    # Nombre d'occurrences du module dans l'état i
    occurrence_count <- sum(train_high[, i] == module)
    # Remplir la matrice d'émission
    emission_matrix[i, j] <- occurrence_count / nrow(train_high)
  }
}

# Fonction de propagation avant (Forward algorithm)
forward_algorithm <- function(O, pi, A, B) {
  T <- length(O)  # Longueur de la séquence d'observations
  n_states <- length(pi)
  
  # Initialisation des alpha_t
  alpha <- matrix(0, nrow = T, ncol = n_states)
  
  # Initialisation de l'étape 1
  alpha[1, ] <- pi * B[, O[1]]  # Pour chaque état, calculer alpha_1(i)
  
  # Récurrence
  for (t in 2:T) {
    for (j in 1:n_states) {
      alpha[t, j] <- sum(alpha[t - 1, ] * A[, j]) * B[j, O[t]]  # Transition et émission
    }
  }
  
  # Terminaison
  P_O_given_lambda <- sum(alpha[T, ])  # Probabilité totale
  return(P_O_given_lambda)
}

# Exemple d'une séquence d'observations pour l'évaluation
# Ici, par exemple, la séquence correspond à un échantillon du train
# Nous utilisons les indices des modules (chaque module étant un indice unique dans la matrice)
# Boucle sur les indices de 1 à 29
for (i in 1:29) {
  # Convertir les modules en indices pour chaque échantillon du train
  O <- match(as.character(test_high[i,]), modules)  # Convertir les modules en indices
  
  # Calculer la probabilité que cette séquence soit générée par le modèle
  P_O_given_lambda <- forward_algorithm(O, initial_state_probabilities, transition_matrix, emission_matrix)
  
  # Afficher la probabilité
  cat("La probabilité que la séquence", i, "soit générée par le modèle HMM est:", P_O_given_lambda, "\n")
}
