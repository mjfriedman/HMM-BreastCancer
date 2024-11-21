# Charger les bibliothèques nécessaires
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
library(readxl)

# Charger les données du fichier Excel
file_path <- "data/top_30_gene_sets.xlsx"  # Remplacez par le fichier voulu (top_N_gene_sets.xlsx)
gene_sets <- read_excel(file_path)

# Vérification rapide des colonnes
cat("Colonnes disponibles :\n")
print(colnames(gene_sets))

# Extraire les gènes de chaque module
gene_modules <- lapply(gene_sets$genes, function(genes_string) {
  unlist(strsplit(genes_string, split = ","))  # Séparer les gènes par des virgules
})

# Nommer chaque module selon la colonne "annotation_id" ou "description"
names(gene_modules) <- gene_sets$annotation_id

# Matrice d'expression nettoyée après l'étape précédente
expr_matrix <- exprs(data_clean)  # Matrice d'expression des gènes

# Initialiser une liste pour stocker les représentants des modules
module_representatives <- list()

# Fonction pour calculer le t-statistique pour un module de gènes
calculate_t_stat <- function(gene_indices, expr_matrix) {
  # Sous-matrice des gènes du module
  module_data <- expr_matrix[gene_indices, , drop = FALSE]
  
  # Taille du module (nombre de gènes)
  n_genes <- nrow(module_data)
  
  # Moyenne et écart-type pour chaque échantillon
  mu <- colMeans(module_data)
  sigma <- apply(module_data, 2, sd)
  
  # Calcul du t-statistique pour chaque échantillon
  t_stat <- (sqrt(n_genes) * mu) / sigma
  return(t_stat)
}

# Calcul des représentants pour chaque module
for (module_name in names(gene_modules)) {
  cat("Processing module:", module_name, "\n")
  
  # Identifier les lignes correspondant aux gènes dans ce module
  module_genes <- gene_modules[[module_name]]
  gene_indices <- which(rownames(expr_matrix) %in% module_genes)
  
  # Vérifier si des gènes du module sont présents dans la matrice
  if (length(gene_indices) == 0) {
    cat("Aucun gène trouvé pour le module :", module_name, "\n")
    next
  }
  
  # Calculer les t-statistiques pour ce module
  module_t_stat <- calculate_t_stat(gene_indices, expr_matrix)
  
  # Stocker les t-statistiques dans la liste
  module_representatives[[module_name]] <- module_t_stat
}

# Combiner les représentants des modules en une matrice
module_representatives_matrix <- do.call(cbind, module_representatives)

# Nommer les colonnes par les noms des modules
colnames(module_representatives_matrix) <- names(module_representatives)

# Résultat : Matrice où chaque colonne est un module et chaque ligne un échantillon
cat("\nMatrice des représentants des modules :\n")
print(module_representatives_matrix)

# Cette matrice peut maintenant être utilisée comme entrée pour le modèle HMM
