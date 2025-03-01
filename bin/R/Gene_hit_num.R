library(ggplot2)
library(reshape2)
library(readr)

# Load the data
csv_file <- "../../results/gene_hits_checkboard.csv"
gene_data <- read_csv(csv_file)

# Convert to long format for ggplot
melted_data <- melt(gene_data, id.vars = "Gene name", variable.name = "Genome", value.name = "Hits")

# Plot heatmap
heatmap_plot <- ggplot(melted_data, aes(x = `Gene name`, y = Genome, fill = Hits)) +
  geom_tile() +
  scale_fill_gradient(low = "cornflowerblue", high = "white") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Gene Hits Heatmap", x = "Genes", y = "Genomes", fill = "Hits")

# Show plot
ggsave("../../results/Gene_hits_heatmap.pdf", plot = heatmap_plot, width = 25, height = 15, limitsize = FALSE)