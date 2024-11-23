# Charger les librairies nécessaires
library(readxl)
library(dplyr)
library(writexl)  
library(tidyr)    

# Fonction pour traiter chaque fichier et générer les résultats
process_gene_sets <- function(gene_sets_file, sample_data_files, metastasis_info_file) {
  
  # Charger les modules de gènes
  gene_sets <- read_excel(gene_sets_file)  # Fichier des modules de gènes
  
  # Initialiser des listes pour stocker les résultats
  t_stat_results <- data.frame(Sample = character(), Module = character(), t_stat = numeric(), stringsAsFactors = FALSE)
  
  # Parcourir les fichiers d'échantillons
  for (sample_data_file in sample_data_files) {
    
    # Charger les données d'expression pour cet échantillon
    sample_data <- read_excel(sample_data_file)  # Fichier avec les expressions des gènes
    
    # Extraire les colonnes d'échantillons (commencent par 'GSM')
    sample_columns <- colnames(sample_data)[grepl("^GSM", colnames(sample_data))]
    
    # Parcourir chaque module dans gene_sets
    for (module in unique(gene_sets$annotation_id)) {
      
      # Extraire les gènes pour le module
      genes_in_module <- gene_sets %>%
        filter(annotation_id == module) %>%
        pull(genes) %>%
        strsplit(", ") %>%
        unlist()  # Liste des gènes dans le module
      
      # Vérifier que les gènes sont présents dans le fichier des échantillons
      module_data <- sample_data %>%
        filter(SYMBOL %in% genes_in_module) %>%
        dplyr::select(SYMBOL, all_of(sample_columns))  # Sélectionner les colonnes des échantillons
      
      # Si des gènes sont trouvés pour ce module
      if (nrow(module_data) > 0) {
        # Convertir les expressions des gènes en matrice (lignes = gènes, colonnes = échantillons)
        expression_matrix <- module_data %>%
          dplyr::select(all_of(sample_columns)) %>%
          as.matrix()
        
        # Calculer la moyenne et l'écart type pour chaque échantillon
        mean_expr <- apply(expression_matrix, 2, mean, na.rm = TRUE)
        sd_expr <- apply(expression_matrix, 2, sd, na.rm = TRUE)
        n_genes <- nrow(expression_matrix)
        
        # Calculer la statistique t pour chaque échantillon
        t_stat <- (sqrt(n_genes) * mean_expr) / sd_expr
        
        # Ajouter les résultats pour le module et l'échantillon dans t_stat_results
        for (i in 1:length(t_stat)) {
          t_stat_results <- rbind(t_stat_results, data.frame(
            Sample = sample_columns[i],
            Module = module,
            t_stat = t_stat[i]
          ))
        }
      }
    }
  }
  
  # Trier les résultats par échantillon et t-statistique (ordre décroissant)
  t_stat_results <- t_stat_results %>%
    arrange(Sample, desc(t_stat))
  
  # Enregistrer les résultats des t-statistiques dans un fichier Excel
  write_xlsx(t_stat_results, path = "data/output/t_stat_results.xlsx")
  
  # Créer le premier fichier avec les modules triés pour chaque échantillon
  module_rankings <- t_stat_results %>%
    group_by(Sample) %>%
    mutate(Rank = rank(-t_stat, ties.method = "first")) %>%
    ungroup() %>%
    dplyr::select(Sample, Rank, Module) %>%
    pivot_wider(names_from = Rank, values_from = Module, names_prefix = "Rank_")
  
  # Charger les informations de métastasis
  metastasis_info <- read_excel(metastasis_info_file)
  
  # Ajouter une colonne 'Risk' selon les critères donnés
  metastasis_info$Risk <- ifelse(metastasis_info$e.dmfs == 1 & metastasis_info$t.dmfs <= 60, "High", "Low")
  
  # Joindre les informations 'Risk' à 'module_rankings'
  module_rankings <- module_rankings %>%
    left_join(metastasis_info %>%
                dplyr::select(geo_accn, Risk), by = c("Sample" = "geo_accn"))
  
  # Enregistrer les classements des modules dans un autre fichier Excel
  write_xlsx(module_rankings, path = "data/output/module_rankings_with_risk.xlsx")
  
  message("Les résultats ont été enregistrés dans 'data/output/module_rankings_with_risk.xlsx' et 'data/output/t_stat_results.xlsx'.")
}

# Liste des fichiers d'échantillons à traiter
sample_data_files <- list(
  "data/output/GSE2034_cleaned_annotated.xlsx",
  "data/output/GSE7390_cleaned_annotated.xlsx",
  "data/output/GSE11121_cleaned_annotated.xlsx"
)

# Nom du fichier des modules de gènes
gene_sets_file <- "data/input/top_30_gene_sets.xlsx"

# Nom du fichier d'informations sur les métastasis
metastasis_info_file <- "data/input/metastasis_info.xlsx"

# Appeler la fonction de traitement
process_gene_sets(gene_sets_file, sample_data_files, metastasis_info_file)
