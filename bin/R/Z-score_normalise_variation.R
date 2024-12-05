# Load necessary libraries
library(ggplot2)

# Read the CSV file
data <- read.csv("../../results/EXTENDED_pangenome_variation_summary.csv", stringsAsFactors = FALSE)

# Extract base gene names (remove "_flank" if present)
data$BASE_GENE <- sub("_flank$", "", data$GENE)

# Summarize the total variations for each base gene
total_variations <- aggregate(NUM_VARIATIONS ~ BASE_GENE, data, sum)

# Merge the total variations back with the original data
data <- merge(data, total_variations, by.x = "BASE_GENE", by.y = "BASE_GENE", suffixes = c("", "_TOTAL"))

# Function to calculate Z-score for each gene individually
calculate_z_score <- function(gene_data) {
  # Calculate the mean and standard deviation for each gene
  mu <- mean(gene_data$NUM_VARIATIONS_TOTAL)
  sigma <- sd(gene_data$NUM_VARIATIONS_TOTAL)
  
  # Normalize and calculate Z-score for each row in the gene's data
  gene_data$NORMALIZED_VAR <- gene_data$NUM_VARIATIONS / gene_data$NUM_VARIATIONS_TOTAL
  gene_data$Z_SCORE <- (gene_data$NORMALIZED_VAR - mu) / sigma
  return(gene_data)
}

# Apply the function to each gene separately using dplyr's group_by
library(dplyr)
data <- data %>%
  group_by(BASE_GENE) %>%
  do(calculate_z_score(.))

# Filter out flank rows for better visualization
non_flank_data <- data[!grepl("_flank$", data$GENE), ]

# Plot Z-scores against genes
ggplot(non_flank_data, aes(x = BASE_GENE, y = Z_SCORE)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Z-scores of Normalized Gene Variations", x = "Gene", y = "Z-score") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave("../../results/z_scores_plot.pdf")