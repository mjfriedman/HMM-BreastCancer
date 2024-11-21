# Clear the environment
rm(list = ls())
# Clear the console
cat("\014")

##### Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install necessary Bioconductor packages
required_packages <- c("affy", "limma", "hgu133a.db")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

library(affy)
library(limma)
library(hgu133a.db)

# Define the list of dataset names and the main data directory
dataset_names <- c("GSE2034", "GSE7390", "GSE11121")
data_dir <- "data/"

# Function for Log Transformation, Normalization, and Arctangent Application
normalize_and_transform <- function(data) {
  # 1: Log Transformation (base-2), adding a small constant (e.g., 1) to avoid log(0)
  expr_matrix <- exprs(data)
  log_transformed_data <- log2(expr_matrix + 1)  # Log2 transformation
  
  # 2: Normalization - Subtract the mean and divide by the standard deviation for each sample
  normalized_data <- t(scale(t(log_transformed_data)))  # Normalize across genes (columns)
  
  # 3: Apply Arctangent function to the normalized data
  transformed_data <- atan(normalized_data)  # Arctangent transformation
  
  # Assign the transformed data back to the expression set
  exprs(data) <- transformed_data
  
  return(data)  # Return the transformed dataset
}

# Function to annotate probes
annotate_probes <- function(data) {
  # Get probeset IDs
  probesets <- featureNames(data)
  
  # Map probesets to gene symbols using hgu133a.db (for HG-U133A)
  gene_symbols <- mapIds(
    hgu133a.db, 
    keys = probesets, 
    column = "SYMBOL", 
    keytype = "PROBEID", 
    multiVals = "first"
  )
  
  # Replace missing annotations with the original probe IDs
  gene_symbols[is.na(gene_symbols)] <- probesets[is.na(gene_symbols)]
  
  # Ensure unique feature names
  featureNames(data) <- make.unique(gene_symbols)
  
  return(data)
}

# Loop to process each dataset
for (dataset_name in dataset_names) {
  cat("\nProcessing dataset:", dataset_name, "\n")
  
  # Check if the dataset directory exists
  dataset_path <- file.path(data_dir, dataset_name)
  if (!dir.exists(dataset_path)) {
    cat("Dataset directory does not exist:", dataset_path, "\n")
    next
  }
  
  # Load the CEL files
  data <- ReadAffy(celfile.path = dataset_path)
  
  # Print initial sample count
  initial_sample_count <- dim(exprs(data))[2]
  cat("Initial number of samples:", initial_sample_count, "\n")
  
  # STEP 1: Remove quality control probes (e.g., those containing "AFFX")
  control_probes <- grep("AFFX", featureNames(data))
  num_control_probes <- length(control_probes)
  cat("Number of control probes identified:", num_control_probes, "\n")
  
  # Remove control probes if found
  if (num_control_probes > 0) {
    data_clean <- data[-control_probes, ]
  } else {
    data_clean <- data
  }
  
  # STEP 2: Remove probes with variance values close to zero
  variance_values <- apply(exprs(data_clean), 1, var)
  low_variance_threshold <- 1e-3  # Threshold for low variance
  low_variance_probes <- which(variance_values < low_variance_threshold)
  
  if (length(low_variance_probes) > 0) {
    data_clean <- data_clean[-low_variance_probes, ]
    cat("Number of probes removed due to low variance:", length(low_variance_probes), "\n")
  } else {
    cat("No probes removed due to low variance.\n")
  }
  
  # STEP 3: Missing Value Handling (Remove probes with > 15% missing data)
  missing_data_proportion <- rowMeans(is.na(exprs(data_clean)))  # NA = missing values
  probes_to_remove <- which(missing_data_proportion > 0.15)
  
  if (length(probes_to_remove) > 0) {
    data_clean <- data_clean[-probes_to_remove, ]
    cat("Number of probes removed due to missing data:", length(probes_to_remove), "\n")
  } else {
    cat("No probes removed due to missing data.\n")
  }
  
  # STEP 4: Apply Log Transformation, Normalization, and Arctangent Transformation
  data_clean <- normalize_and_transform(data_clean)
  
  # STEP 5: Annotate probes
  data_clean <- annotate_probes(data_clean)
  
  # Final sample count after cleaning
  final_sample_count <- dim(exprs(data_clean))[2]
  cat("Final number of samples after cleaning:", final_sample_count, "\n")
  
  # Save the cleaned and annotated data
  save_file <- paste0(dataset_name, "_cleaned_annotated.RData")
  save(data_clean, file = save_file)
  cat("Saved cleaned data to:", save_file, "\n")
}

cat("\nAll datasets have been processed and analyzed.\n")
