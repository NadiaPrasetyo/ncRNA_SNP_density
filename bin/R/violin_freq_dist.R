library(ggplot2)

# Load data
file_path <- "../../data/gnomad_gene_rand_freq.csv"  # Replace with your actual file path
data <- read.csv(file_path, stringsAsFactors = FALSE)

# Calculate log10 of allele frequency
data$log_freq <- log10(data$allele_freq)

# Categorize genes vs randoms
data$category <- ifelse(grepl("^Random\\d+", data$gene_name), "Random", "Gene")

# Plot distribution using jitter scatter and violin plot
plot <-ggplot(data, aes(x = category, y = log_freq, fill = category, color = category)) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  geom_violin(alpha = 0.5, position = position_dodge(width = 0.9)) +
  theme_minimal() +
  scale_fill_manual(values = c("Random" = "darkblue", "Gene" = "darkred")) +
  scale_color_manual(values = c("Random" = "cornflowerblue", "Gene" = "firebrick1")) +
  labs(title = "Distribution of log10(Frequency) of Genes vs Randoms",
       x = "Category",
       y = "log10(Frequency)")

ggsave("../../results/frequency_20gene_100rand_dist.pdf", plot = plot, height = 10, width = 8)