# Load necessary libraries
library(dplyr)   # For data manipulation
library(ggplot2) # For plotting
library(scales)   # For scaling/normalization
library(ggtext)   # For colored text in ggplot
library(ggrepel)  # For better text labeling in scatter plots

# Read the CSV file
df <- read.csv('../../results/EXTENDED_SNP_variation_summary.csv')

# Ensure NUM_VARIATIONS and NUM_FLANKS are in the dataset and handle missing values
df$NUM_VARIATIONS[is.na(df$NUM_VARIATIONS)] <- 0
df$NUM_FLANKS[is.na(df$NUM_FLANKS)] <- 0  # Handle missing values in NUM_FLANKS

# Calculate the SNP Enrichment
df$Enrichment <- (df$NUM_VARIATIONS / df$GENE_LENGTH) / 
  (df$NUM_FLANKS / (df$TOTAL_LENGTH - df$GENE_LENGTH))

# Calculate metrics for scatter plot
df$Gene_Var_Metric <- df$NUM_VARIATIONS / df$GENE_LENGTH
df$Flank_Var_Metric <- df$NUM_FLANKS / (df$TOTAL_LENGTH - df$GENE_LENGTH)

# List of special genes
special_genes <- c(
  "FAM30A", "LINC01671", "lnc-SLCO4A1-8", "MIR4538", "SNAR-A1", "SNAR-B2", 
  "SNAR-C3", "SNAR-C4", "SNAR-G1", "SNAR-G2", "TRE-TTC5-1", "TRG-CCC6-1", 
  "TRG-CCC4-1", "TRV-CAC5-1"
)

# Modify the gene names to include HTML-style color tags for special genes
df$GENE_LABEL <- ifelse(df$GENE %in% special_genes,
                        paste0("<span style='color:red;'>", df$GENE, "</span>"),
                        df$GENE)

# View the resulting dataframe with enrichment
head(df)

# Plot a histogram with unsorted data
plot_unsorted <- ggplot(df, aes(x = GENE, y = Enrichment)) +
  geom_bar(stat = 'identity', fill = 'indianred') +
  theme_minimal() +
  labs(title = 'Enrichment of Variations in Gene Length vs Total Length',
       x = 'Gene',
       y = 'Enrichment') +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +  # Apply markdown for special genes
  scale_x_discrete(labels = df$GENE_LABEL)  # Use modified labels for x-axis

# Save the unsorted plot
ggsave('../../results/SNP_Enrichment_histogram.pdf', plot = plot_unsorted, width = 15)

# Sort the dataframe based on Enrichment (high to low)
df_sorted <- df %>% arrange(desc(Enrichment))

# Reorder the GENE factor based on the sorted Enrichment
df_sorted$GENE <- factor(df_sorted$GENE, levels = df_sorted$GENE)

# Plot a sorted histogram based on Z-scores with special gene labels
plot_sorted <- ggplot(df_sorted, aes(x = reorder(GENE, -Enrichment), y = Enrichment)) +
  geom_bar(stat = 'identity', fill = 'indianred') +
  theme_minimal() +
  labs(title = 'Sorted Enrichment of Variations in Gene Length vs Total Length',
       x = 'Gene',
       y = 'Enrichment') +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +  # Apply markdown for special genes
  scale_x_discrete(labels = df_sorted$GENE_LABEL)  # Use modified labels for x-axis

# Save the sorted plot
ggsave('../../results/SNP_Enrichment_sorted_histogram.pdf', plot = plot_sorted, width = 15)

# Create a scatter plot for (Gene Variations / Gene Length) vs. (Flank Variations / (Total Length - Gene Length))
scatter_plot <- ggplot(df, aes(x = Gene_Var_Metric, y = Flank_Var_Metric, label = GENE)) +
  geom_point(aes(color = GENE %in% special_genes), size = 3) +
  scale_color_manual(values = c("indianred", "black"), labels = c("Other Genes", "Special Genes")) +
  geom_text_repel(aes(label = GENE), size = 3) +
  theme_minimal() +
  labs(title = "Scatter Plot of Gene Variations vs Flank Variations",
       x = "Gene Variations / Gene Length",
       y = "Flank Variations / (Total Length - Gene Length)",
       color = "Gene Category")

# Save the scatter plot
ggsave('../../results/Gene_Var_vs_Flank_Var_Scatter_All_Labels.pdf', plot = scatter_plot, width = 10, height = 8)
