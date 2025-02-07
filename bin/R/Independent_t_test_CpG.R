# Load necessary library
library(dplyr)

# Function to compile methylation percentages and perform t-tests
compile_and_test <- function(data) {
  # Initialize a list to store methylation percentages
  compiled_data <- list()
  
  # Compile methylation percentages for each gene and tissue
  unique_genes <- unique(data$GeneName)
  for (gene in unique_genes) {
    gene_data <- data %>% filter(GeneName == gene)
    unique_tissues <- unique(gene_data$Tissue)
    
    compiled_data[[gene]] <- lapply(unique_tissues, function(tissue) {
      tissue_data <- gene_data %>% filter(Tissue == tissue) %>% pull(Methylation_Percentage)
      list(Tissue = tissue, Methylation_Percentage = tissue_data)
    })
    names(compiled_data[[gene]]) <- unique_tissues
  }
  
  # Print compiled data for diagnostics
  print("Compiled methylation percentages for each gene and tissue:")
  print(compiled_data)
  
  # Initialize results data frame for t-test results
  results <- data.frame(Gene = character(), Tissues = character(), Difference = numeric(), T_Statistic = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
  
  # Perform t-tests between tissues for each gene
  for (gene in names(compiled_data)) {
    gene_tissues <- compiled_data[[gene]]
    tissue_combinations <- combn(names(gene_tissues), 2, simplify = TRUE)
    
    for (i in seq_len(ncol(tissue_combinations))) {
      tissue1 <- tissue_combinations[1, i]
      tissue2 <- tissue_combinations[2, i]
      
      tissue1_data <- gene_tissues[[tissue1]]$Methylation_Percentage
      tissue2_data <- gene_tissues[[tissue2]]$Methylation_Percentage
      
      # Perform t-test only if both tissues have more than one observation
      if (length(tissue1_data) > 1 && length(tissue2_data) > 1) {
        t_test <- t.test(tissue1_data, tissue2_data, var.equal = FALSE)
        
        # Calculate mean difference and extract t-statistic
        mean_diff <- mean(tissue1_data) - mean(tissue2_data)
        t_stat <- t_test$statistic
        
        # Add results to output
        results <- rbind(results, data.frame(
          Gene = gene,
          Tissues = paste(tissue1, "vs", tissue2),
          Difference = mean_diff,
          T_Statistic = t_stat,
          P_Value = t_test$p.value,
          stringsAsFactors = FALSE
        ))
      } else {
        message(paste("Skipping t-test for gene", gene, "between tissues", tissue1, "and", tissue2, 
                      "- insufficient data (tissue1:", length(tissue1_data), "tissue2:", length(tissue2_data), ")"))
      }
    }
    
    # Perform t-tests for each tissue vs the rest (combined)
    for (tissue in names(gene_tissues)) {
      tissue_data <- gene_tissues[[tissue]]$Methylation_Percentage
      other_data <- unlist(lapply(gene_tissues[names(gene_tissues) != tissue], function(x) x$Methylation_Percentage))
      
      # Perform t-test only if both tissue and other combined have more than one observation
      if (length(tissue_data) > 1 && length(other_data) > 1) {
        t_test <- t.test(tissue_data, other_data, var.equal = FALSE)
        
        # Calculate mean difference and extract t-statistic
        mean_diff <- mean(tissue_data) - mean(other_data)
        t_stat <- t_test$statistic
        
        # Add results to output
        results <- rbind(results, data.frame(
          Gene = gene,
          Tissues = paste(tissue, "vs others"),
          Difference = mean_diff,
          T_Statistic = t_stat,
          P_Value = t_test$p.value,
          stringsAsFactors = FALSE
        ))
      } else {
        message(paste("Skipping t-test for gene", gene, "for tissue", tissue, "vs others", 
                      "- insufficient data (tissue:", length(tissue_data), "others:", length(other_data), ")"))
      }
    }
  }
  
  return(results)
}

# Read in the CSV file
data <- read.csv("../../data/CpG_methylation_data.csv", stringsAsFactors = FALSE)

# Compile and perform the analysis
output <- compile_and_test(data)

# Write output to a CSV file
write.csv(output, "../../results/CpG_difference.csv", row.names = FALSE)
