# Load required libraries
library(dplyr)
library(ggplot2)
library(stats)  # For KS test

# Load data
data <- read.csv("../../data/CHH_methylation_data.csv")

# Convert Methylation_Percentage to numeric
data$Methylation_Percentage <- as.numeric(data$Methylation_Percentage)

# Check structure of the data
print(str(data))

# Create a dataframe to store KS test results
ks_results <- data.frame(Tissue = character(), KS_Statistic = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Part 1: Tissue-based analysis
# Get the unique tissue types
tissues <- unique(data$Tissue)

for (tissue in tissues) {
  tissue_data <- filter(data, Tissue == tissue)
  
  # Create a GeneType column to differentiate genes and randoms
  tissue_data <- tissue_data %>%
    mutate(GeneType = ifelse(grepl("^Random\\d+", GeneName), "Random", "Gene"))
  
  print(paste("Processing Tissue:", tissue))
  
  gene_methylation <- tissue_data$Methylation_Percentage[tissue_data$GeneType == "Gene"]
  random_methylation <- tissue_data$Methylation_Percentage[tissue_data$GeneType == "Random"]
  
  if (length(gene_methylation) > 0 && length(random_methylation) > 0) {
    ks_test <- ks.test(gene_methylation, random_methylation)
    
    ks_results <- rbind(ks_results, data.frame(Tissue = tissue, KS_Statistic = ks_test$statistic, P_Value = ks_test$p.value))
    
    summary_stats <- tissue_data %>% group_by(GeneType) %>% summarise(
      Median = median(Methylation_Percentage, na.rm = TRUE),
      Mean = mean(Methylation_Percentage, na.rm = TRUE)
    )
    
    plot <- ggplot(tissue_data, aes(x = GeneType, y = Methylation_Percentage, color = GeneType)) +
      geom_jitter(alpha = 0.7, size = 3, width = 0.1) +
      geom_point(data = summary_stats, aes(x = GeneType, y = Median), color = "black", size = 4, shape = 18) +
      geom_point(data = summary_stats, aes(x = GeneType, y = Mean), color = "red", size = 4, shape = 17) +
      geom_text(data = summary_stats, aes(x = GeneType, y = Mean, label = round(log10(Mean), 2)), vjust = -1, color = "red") +
      scale_y_continuous(trans = "log10", labels = scales::log10_format()) +
      labs(
        title = paste("Methylation Percentage in", tissue),
        x = "Gene Type",
        y = "log10(Methylation Percentage)"
      ) +
      theme_minimal()
    
    ggsave(paste0("../../results/DNAme/", tissue, "_CHH_scatter_jitter_plot.pdf"), plot)
  }
}

write.csv(ks_results, "../../results/DNAme/_CHH_ks_test_results.csv", row.names = FALSE)



# Part 2: Individual gene analysis
# Function to filter randoms based on gene's length and CG content
filter_randoms <- function(gene, tissue_data) {
  randoms <- tissue_data %>%
    filter(
      grepl("^Random\\d+", GeneName),
      abs(Length - gene$Length) <= 200,
      abs(Median_CG_Content - gene$Median_CG_Content) <= 20
    )
  print(paste("Filtering Randoms for Gene:", gene$GeneName[1]))
  print(paste("Random Rows Found:", nrow(randoms)))
  return(randoms)
}

# Get unique gene names excluding the "Random#" genes
genes <- unique(data$GeneName[!grepl("^Random\\d+", data$GeneName)])

# Define the desired order for GeneType on the x-axis
desired_order <- c("blood", "random blood", "embryo", "random embryo", 
                   "brain", "random brain", "epithelium", "random epithelium")

# Loop through each gene for individual analysis
for (gene in genes) {
  gene_data <- filter(data, GeneName == gene)
  
  # Create an empty data frame to store the results for each gene and tissue
  all_data <- data.frame(
    Tissue = character(),
    Methylation_Percentage = numeric(),
    GeneType = character()
  )
  
  print(paste("Processing Gene:", gene))
  
  # Loop through each tissue to create random comparisons for the gene
  for (tissue in tissues) {
    tissue_data <- filter(data, Tissue == tissue)
    randoms <- filter_randoms(gene_data[1, ], tissue_data)
    
    # Combine gene data with the corresponding randoms
    all_data <- bind_rows(
      all_data,
      tissue_data %>% filter(GeneName == gene) %>% mutate(GeneType = tissue),
      randoms %>% mutate(GeneType = paste("random", tissue, sep = " "))
    )
  }
  
  # If there is data for the gene, create and save the scatter plot
  if (nrow(all_data) > 0) {
    print(paste("Rows in Final Plot Data for Gene", gene, ":", nrow(all_data)))
    
    # Ensure the order of GeneType on the x-axis
    all_data$GeneType <- factor(all_data$GeneType, levels = desired_order)
    
    # Calculate summary statistics
    summary_stats <- all_data %>% group_by(GeneType) %>% summarise(
      Median = median(Methylation_Percentage, na.rm = TRUE),
      Mean = mean(Methylation_Percentage, na.rm = TRUE)
    )
    
    plot <- ggplot(all_data, aes(x = GeneType, y = Methylation_Percentage, color = GeneType)) +
      geom_jitter(width = 0.2, alpha = 0.7) +
      geom_point(data = summary_stats, aes(x = GeneType, y = Median), color = "black", size = 4, shape = 18) +
      geom_point(data = summary_stats, aes(x = GeneType, y = Mean), color = "red", size = 4, shape = 17) +
      geom_text(data = summary_stats, aes(x = GeneType, y = Mean, label = round(Mean, 2)), vjust = -1, color = "red") +
      scale_y_continuous(trans = "log10") +
      labs(
        title = paste(gene, "CHH vs Random"),
        x = "Gene Type",
        y = "log10(Methylation Percentage)"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
    ggsave(paste0("../../results/DNAme/", gene, "_CHH_individual_scatter_plot.pdf"), plot)
  }
}
