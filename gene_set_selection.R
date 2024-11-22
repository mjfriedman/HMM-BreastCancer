# Installer les packages nécessaires si ce n'est pas déjà fait
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# Charger les packages
library(readxl)
library(dplyr)
library(openxlsx)

# Lire le fichier Excel
file_path <- 'data/enrich-374GeneSetResults-GO_BP.xlsx'
df <- read_excel(file_path)

# Convertir les colonnes numériques en types numériques
numeric_columns <- c('input_size', 'term_genes', 'universe', 'pval', 'pval_adj', 'relative_enrichment', 'annotsbias')
df[numeric_columns] <- lapply(df[numeric_columns], as.numeric)

# Filtrer les ensembles de gènes par valeur p ajustée
threshold_pval_adj <- 0.05
filtered_df <- df %>% filter(pval_adj < threshold_pval_adj)

# Trier par pval_adj (croissant)
sorted_df <- filtered_df %>% arrange(pval_adj)

# Choisir le nombre de gènes à sélectionner
num_genes <- 30  # Changez cette valeur selon vos besoins

# Sélectionner les meilleurs ensembles de gènes
top_gene_sets <- head(sorted_df, num_genes)

# Créer le nom de fichier de sortie en fonction du nombre de gènes
output_file_name <- paste0('data/top_', num_genes, '_gene_sets.xlsx')

# Exporter les résultats vers un nouveau fichier Excel
write.xlsx(top_gene_sets, output_file_name, rowNames = FALSE)

print(paste("Les meilleurs ensembles de gènes ont été sélectionnés et sauvegardés dans", output_file_name))
