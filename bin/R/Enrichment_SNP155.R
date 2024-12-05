# Load necessary libraries
library(dplyr)   # For data manipulation
library(ggplot2) # For plotting
library(scales)   # For scaling/normalization
library(ggtext)   # For colored text in ggplot

# Read the CSV file
df <- read.csv('../../results/EXTENDED_SNP_variation_summary.csv')

# Ensure NUM_VARIATIONS is in the dataset and handle missing values
df$NUM_VARIATIONS[is.na(df$NUM_VARIATIONS)] <- 0
df$NUM_FLANKS[is.na(df$NUM_FLANKS)] <- 0  # Handle missing values in NUM_FLANKS

# Calculate the SNP_Enrichment for each gene
df$Enrichment <- (df$NUM_VARIATIONS / df$GENE_LENGTH) / 
  ((df$NUM_VARIATIONS + df$NUM_FLANKS) / df$TOTAL_LENGTH)

# Normalize the enrichment values (Z-score normalization)
df$Enrichment_ZScore <- scale(df$Enrichment)

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

# View the resulting dataframe with enrichment and Z-scores
head(df)

# Plot a histogram with unsorted data (raw enrichment Z-scores)
plot_unsorted <- ggplot(df, aes(x = GENE, y = Enrichment_ZScore)) +
  geom_bar(stat = 'identity', fill = 'indianred') +
  theme_minimal() +
  labs(title = 'Enrichment of Variations in Gene Length vs Total Length',
       x = 'Gene',
       y = 'Z-score of Enrichment') +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +  # Apply markdown for special genes
  scale_x_discrete(labels = df$GENE_LABEL)  # Use modified labels for x-axis

# Save the unsorted plot
ggsave('../../results/SNP_Enrichment_histogram.pdf', plot = plot_unsorted, width = 15)

# Sort the dataframe based on Enrichment_ZScore (high to low)
df_sorted <- df %>% arrange(Enrichment_ZScore)

# Reorder the GENE factor based on the sorted Enrichment_ZScore
df_sorted$GENE <- factor(df_sorted$GENE, levels = df_sorted$GENE)

# Plot a sorted histogram based on Z-scores with special gene labels
plot_sorted <- ggplot(df_sorted, aes(x = reorder(GENE, Enrichment_ZScore), y = Enrichment_ZScore)) +
  geom_bar(stat = 'identity', fill = 'indianred') +
  theme_minimal() +
  labs(title = 'Sorted Enrichment of Variations in Gene Length vs Total Length',
       x = 'Gene',
       y = 'Z-score of Enrichment') +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +  # Apply markdown for special genes
  scale_x_discrete(labels = df_sorted$GENE_LABEL)  # Use modified labels for x-axis

# Save the sorted plot
ggsave('../../results/SNP_Enrichment_sorted_histogram.pdf', plot = plot_sorted, width = 15)
