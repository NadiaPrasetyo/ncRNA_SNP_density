# Load necessary libraries
library(dplyr)   # For data manipulation
library(ggplot2) # For plotting
library(scales)   # For scaling/normalization
library(ggtext)  # For colored text in ggplot

# Read the CSV file
df <- read.csv('../../results/EXTENDED_pangenome_variation_summary.csv')

# Split variation types from the 'VARIATION_TYPES' field to get the number of variations
# Assuming the number of variations is extracted from the string (e.g., "snp(4)")
df$NUM_VARIATIONS <- sapply(df$VARIATION_TYPES, function(x) {
  sum(as.numeric(gsub(".*\\((\\d+)\\).*", "\\1", unlist(strsplit(x, ",\\s*")))))
})

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
  ((merged_df$NUM_VARIATIONS_gene + merged_df$NUM_VARIATIONS_flank) / merged_df$TOTAL_LENGTH_gene)

# Normalize the enrichment values (Z-score normalization)
merged_df$Enrichment_ZScore <- scale(merged_df$Enrichment)

# View the resulting dataframe with enrichment and Z-scores
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
plot <- ggplot(merged_df, aes(x = GENE_BASE, y = Enrichment_ZScore)) +
  geom_bar(stat = 'identity', fill = 'steelblue') +
  theme_minimal() +
  labs(title = 'Enrichment of Variations in Gene Length vs Total Length',
       x = 'Gene',
       y = 'Z-score of Enrichment') +
  theme(axis.text.x = ggtext::element_markdown(angle = 90, hjust = 1)) +  # This applies markdown to the x-axis text
  scale_x_discrete(labels = merged_df$GENE_LABEL)   # Use modified labels for x-axis

ggsave("../../results/Enrichment_pangenome_var.pdf", plot = plot)