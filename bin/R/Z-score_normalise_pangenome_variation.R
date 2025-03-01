# Load necessary libraries
library(dplyr)   # For data manipulation
library(ggplot2) # For plotting
library(scales)  # For scaling/normalization
library(ggtext)  # For colored text in ggplot
library(ggrepel) # For better text labeling in scatter plots

# Read the CSV file
df <- read.csv('../../results/EXTENDED_pangenome_variation_summary.csv')

# Split variation types from the 'VARIATION_TYPES' field to get the number of variations
# Assuming the number of variations is extracted from the string (e.g., "snp(4)")
df$NUM_VARIATIONS <- sapply(df$VARIATION_TYPES, function(x) {
  sum(as.numeric(gsub(".*\\((\\d+)\\).*", "\\1", unlist(strsplit(x, ",\\s*")))))
})

# Fill NA values in NUM_VARIATIONS with 0
df$NUM_VARIATIONS[is.na(df$NUM_VARIATIONS)] <- 0

# Separate genes and flanks
genes <- df[!grepl("_flank", df$GENE), ]       # Rows without "_flank" (genes)
flanks <- df[grepl("_flank", df$GENE), ]       # Rows with "_flank" (flanks)

# Merge the gene and flank data based on the base gene name
# Assuming the gene name is the part before '_flank' (e.g., 'BLACE' for 'BLACE_flank')
genes$GENE_BASE <- sub("_flank$", "", genes$GENE)
flanks$GENE_BASE <- sub("_flank$", "", flanks$GENE)

# Merge gene data with flank data by the base gene name
merged_df <- merge(genes, flanks, by = "GENE_BASE", suffixes = c("_gene", "_flank"))

# Calculate the enrichment ratio for each gene
merged_df$Enrichment <- (merged_df$NUM_VARIATIONS_gene / merged_df$GENE_LENGTH_gene) /
  (merged_df$NUM_VARIATIONS_flank / (merged_df$TOTAL_LENGTH_gene - merged_df$GENE_LENGTH_gene))

# Calculate metrics for scatter plot
merged_df$Gene_Var_Metric <- merged_df$NUM_VARIATIONS_gene / merged_df$GENE_LENGTH_gene
merged_df$Flank_Var_Metric <- merged_df$NUM_VARIATIONS_flank / 
  (merged_df$TOTAL_LENGTH_gene - merged_df$GENE_LENGTH_gene)

# View the resulting dataframe with enrichment
head(merged_df)

# List of special genes
special_genes <- c(
  "FAM30A", "LINC01671", "lnc-SLCO4A1-8", "MIR4538", "SNAR-A1", "SNAR-B2", 
  "SNAR-C3", "SNAR-C4", "SNAR-G1", "SNAR-G2", "TRE-TTC5-1", "TRG-CCC6-1", 
  "TRG-CCC4-1", "TRV-CAC5-1"
)

# Modify the gene names to include HTML-style color tags for special genes
merged_df$GENE_LABEL <- ifelse(merged_df$GENE_BASE %in% special_genes,
                               paste0("<span style='color:red;'>", merged_df$GENE_BASE, "</span>"),
                               merged_df$GENE_BASE)

# Plot with colored x-axis labels using ggtext
plot <- ggplot(merged_df, aes(x = GENE_BASE, y = Enrichment)) +
  geom_bar(stat = 'identity', fill = 'steelblue') +
  theme_minimal() +
  labs(title = 'Enrichment of Variations in Gene Length vs Total Length',
       x = 'Gene',
       y = 'Enrichment') +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +  # This applies markdown to the x-axis text
  scale_x_discrete(labels = merged_df$GENE_LABEL)   # Use modified labels for x-axis

# Save the plot
ggsave("../../results/Enrichment_pangenome_var.pdf", plot = plot, width = 15)

# Sort the dataframe based on Enrichment (high to low)
merged_df_sorted <- merged_df %>% arrange(desc(Enrichment))

# Reorder the GENE_BASE factor based on the sorted Enrichment
merged_df_sorted$GENE_BASE <- factor(merged_df_sorted$GENE_BASE, levels = merged_df_sorted$GENE_BASE)

# Plot a sorted bar plot based on enrichment with special gene labels
plot_sorted <- ggplot(merged_df_sorted, aes(x = reorder(GENE_BASE, -Enrichment), y = Enrichment)) +
  geom_bar(stat = 'identity', fill = 'steelblue') +
  theme_minimal() +
  labs(title = 'Enrichment of Variations in Gene Length vs Total Length',
       x = 'Gene',
       y = 'Enrichment') +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +  # Apply markdown for special genes
  scale_x_discrete(labels = merged_df_sorted$GENE_LABEL)  # Use modified labels for x-axis

# Save the sorted plot
ggsave("../../results/Enrichment_pangenome_var_sorted.pdf", plot = plot_sorted, width = 15)

# Create a scatter plot for (Gene Variations / Gene Length) vs. (Flank Variations / Flank Length)
scatter_plot <- ggplot(merged_df, aes(x = Gene_Var_Metric, y = Flank_Var_Metric, label = GENE_BASE)) +
  geom_point(aes(color = GENE_BASE %in% special_genes), size = 3) +
  scale_color_manual(values = c("steelblue", "red"), labels = c("Other Genes", "Special Genes")) +
  geom_text_repel(aes(label = GENE_BASE), size = 3) +
  theme_minimal() +
  labs(title = "Scatter Plot of Gene Variations vs Flank Variations",
       x = "Gene Variations / Gene Length",
       y = "Flank Variations / Flank Length",
       color = "Gene Category")

# Save the scatter plot
ggsave("../../results/Pangenome_Var_vs_Flank_Var_Scatter_All_Labels.pdf", plot = scatter_plot, width = 10, height = 8)
