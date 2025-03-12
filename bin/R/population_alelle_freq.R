library(ggplot2)
library(reshape2)
library(pheatmap)

# Load data
data <- read.csv("../../data/gnomad_region_data.csv")

# Ensure results directory exists
if (!dir.exists("../../results/pop_freq")) {
  dir.create("../../results/pop_freq")
}

# Mapping of population abbreviations to full names
population_map <- c(
  afr = "African/African American",
  ami = "Amish",
  amr = "Latino/Admixed American",
  asj = "Ashkenazi Jewish",
  eas = "East Asian",
  fin = "Finnish",
  nfe = "Non-Finnish European",
  mid = "Middle Eastern",
  sas = "South Asian",
  remaining = "Other (population not assigned)"
)

# List of population prefixes matching column names
populations <- names(population_map)
existing_populations <- populations[sapply(paste0(populations, "_af"), function(col) col %in% colnames(data))]

# Compute log10(frequency) safely
for (pop in existing_populations) {
  af_col <- paste0(pop, "_af")
  data[[paste0("log10_", pop, "_af")]] <- log10(ifelse(is.na(data[[af_col]]), 1e-10, data[[af_col]]) + 1e-10)
}

unique_genes <- unique(data$gene_name)
all_genes_avg_freq <- data.frame()

for (gene in unique_genes) {
  print(paste("Processing gene:", gene))
  gene_data <- subset(data, gene_name == gene)
  gene_data <- subset(gene_data, grepl("^n\\.\\d+[ACGT]>[ACGT]$", gene_data$variation_consequence))
  if (nrow(gene_data) == 0) next
  
  gene_data_melted <- melt(gene_data, id.vars = c("variation_id", "variation_consequence"),
                           measure.vars = paste0("log10_", existing_populations, "_af"),
                           variable.name = "Population", value.name = "Log10_AF")
  
  gene_data_melted$Population <- population_map[gsub("log10_|_af", "", gene_data_melted$Population)]
  
  # Cumulative Distribution Plot
  cdf_plot <- ggplot(gene_data_melted, aes(x = Log10_AF, color = Population)) +
    stat_ecdf(size = 1.2) +
    labs(title = paste("Cumulative Distribution of log10 Allele Frequencies for", gene),
         x = "log10(Allele Frequency)", y = "Cumulative Proportion") +
    theme_minimal(base_size = 22) +
    theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22),
          plot.title = element_text(size = 24, face = "bold"))
  
  ggsave(filename = paste0("../../results/pop_freq/", gene, "_population_cdf.pdf"), plot = cdf_plot, width = 25, height = 10)
  
  # Heatmap Preparation
  heatmap_data <- reshape2::dcast(gene_data_melted, Population ~ variation_consequence, value.var = "Log10_AF")
  rownames(heatmap_data) <- heatmap_data$Population
  heatmap_data$Population <- NULL
  
  heatmap_matrix <- as.matrix(heatmap_data)
  heatmap_matrix[is.na(heatmap_matrix)] <- -10  # Assign low log value for missing variants
  
  color_breaks <- c(-10, seq(min(heatmap_matrix[heatmap_matrix > -10], na.rm = TRUE), max(heatmap_matrix, na.rm = TRUE), length.out = 100))
  color_palette <- c("grey", colorRampPalette(c("blue", "red"))(99))
  
  pdf(paste0("../../results/pop_freq/", gene, "_population_heatmap.pdf"), width = 25, height = 10)
  pheatmap(heatmap_matrix, cluster_rows = FALSE, cluster_cols = FALSE,
           main = paste("Heatmap of log10 Allele Frequencies for", gene),
           fontsize = 25,
           fontsize_row = 20, fontsize_col = 18,
           breaks = color_breaks, color = color_palette,
           labels_col = ifelse(apply(heatmap_matrix, 2, max, na.rm = TRUE) > -1, colnames(heatmap_matrix), ""))
  dev.off()
  
  gene_avg_freq <- aggregate(Log10_AF ~ Population, data = gene_data_melted, FUN = mean, na.rm = TRUE)
  gene_avg_freq$Gene <- gene
  all_genes_avg_freq <- rbind(all_genes_avg_freq, gene_avg_freq)
}

# Cumulative Distribution Plot for Combined Data
combined_cdf_plot <- ggplot(all_genes_avg_freq, aes(x = Log10_AF, color = Population)) +
  stat_ecdf(size = 1.5) +
  labs(title = "Cumulative Distribution of log10 Allele Frequencies Across Populations",
       x = "log10(Allele Frequency)", y = "Cumulative Proportion") +
  theme_minimal(base_size = 22) +
  theme(legend.text = element_text(size = 20), legend.title = element_text(size = 22),
        plot.title = element_text(size = 24, face = "bold"))

ggsave(filename = "../../results/pop_freq/Combined_Population_CDF.pdf", plot = combined_cdf_plot, width = 30, height = 12)
