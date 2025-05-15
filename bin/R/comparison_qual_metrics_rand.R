# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dplyr)
library(gridExtra)
library(grid)

# Set paths to data files
gene_data_path <- "../../data/variant_qual_metrics.csv"
random_data_path <- "../../data/variant_qual_metrics_random.csv"

# Read data from CSV files
gene_data <- read.csv(gene_data_path, stringsAsFactors = FALSE)
random_data <- read.csv(random_data_path, stringsAsFactors = FALSE)

# Select and clean numeric columns
numeric_columns <- c("Mean_Coverage", "Allele.frequency", "popmax_af", "SiteQuality", "AS_FS", "AS_MQRankSum", "AS_pab_max", "AS_ReadPosRankSum", "AS_SOR", "AS_VarDP")
gene_data[numeric_columns] <- lapply(gene_data[numeric_columns], function(x) suppressWarnings(as.numeric(as.character(x))))
random_data[numeric_columns] <- lapply(random_data[numeric_columns], function(x) suppressWarnings(as.numeric(as.character(x))))

# Remove rows with NA values
gene_data <- gene_data %>% filter(complete.cases(.))
random_data <- random_data %>% filter(complete.cases(.))

# Apply log transformation
log_transform_metrics <- c("Allele.frequency", "SiteQuality", "AS_VarDP", "popmax_af")
gene_data[log_transform_metrics] <- lapply(gene_data[log_transform_metrics], function(x) log10(x + 1))
random_data[log_transform_metrics] <- lapply(random_data[log_transform_metrics], function(x) log10(x + 1))

# Perform KS Test for each metric
ks_results <- data.frame(Metric = character(), Statistic = numeric(), P_Value = numeric())
for (metric in numeric_columns) {
  ks_test <- ks.test(gene_data[[metric]], random_data[[metric]])
  ks_results <- rbind(ks_results, data.frame(Metric = metric, Statistic = ks_test$statistic, P_Value = ks_test$p.value))
}

# Save KS test results to CSV
write.csv(ks_results, "../../results/qual_metrics/ks_test_results.csv", row.names = FALSE)

# Distribution plot comparison
metrics <- numeric_columns
for (metric in metrics) {
  combined_data <- rbind(
    data.frame(Value = gene_data[[metric]], Type = "Gene"),
    data.frame(Value = random_data[[metric]], Type = "Random")
  )
  
  # Adjust plot title for log-transformed metrics
  x_label <- ifelse(metric %in% log_transform_metrics, paste("Log10-transformed", metric), metric)
  
  # Adjust for density plot when most values are zeros
  p <- ggplot(combined_data, aes(x = Value, fill = Type)) +
    geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.6, bins = 30) +
    geom_density(alpha = 0.7, adjust = 1.5) +  # Smoothing factor added
    labs(title = paste("Distribution of", metric, "(Gene vs Random)"), x = x_label, y = "Density") +
    theme_minimal()
  
  ggsave(paste0("../../results/qual_metrics/comparison_", metric, ".pdf"), plot = p)
}
