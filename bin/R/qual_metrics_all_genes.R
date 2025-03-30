# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dplyr)
library(gridExtra)
library(grid)

# Set the path to your data
data_path <- "../../data/variant_qual_metrics.csv"

# Read the data from the CSV
data <- read.csv(data_path, stringsAsFactors = FALSE)

# Manually select columns that are numeric
numeric_columns <- c("Mean_Coverage", "Allele.frequency", "popmax_af", "SiteQuality", 
                     "AS_FS", "AS_MQRankSum", "AS_pab_max", 
                     "AS_ReadPosRankSum", "AS_SOR", "AS_VarDP")

# Convert all numeric columns to numeric, handling non-numeric values
data[numeric_columns] <- lapply(data[numeric_columns], function(x) {
  suppressWarnings(as.numeric(as.character(x)))
})

# Remove rows with NA values in numeric columns
data_clean <- data %>%
  filter(complete.cases(data[numeric_columns]))

# Convert data to long format
long_data <- data_clean %>%
  select(GeneName, Variation_ID, SiteQuality, AS_FS, AS_MQRankSum, AS_pab_max, 
         AS_ReadPosRankSum, AS_SOR, AS_VarDP) %>%
  melt(id.vars = c("GeneName", "Variation_ID"), variable.name = "Metric", value.name = "Value")

# Log-transform SiteQuality and AS_VarDP
long_data$Value_log <- long_data$Value
long_data$Value_log[long_data$Metric == "SiteQuality"] <- log(long_data$Value[long_data$Metric == "SiteQuality"] + 1)
long_data$Value_log[long_data$Metric == "AS_VarDP"] <- log(long_data$Value[long_data$Metric == "AS_VarDP"] + 1)

# Define consistent metric order
metric_order <- c("SiteQuality", "AS_FS", "AS_MQRankSum", "AS_pab_max", 
                  "AS_ReadPosRankSum", "AS_SOR", "AS_VarDP")

# Heatmap color palette
heatmap_colors <- colorRampPalette(c("blue", "white", "red"))(100)

# Loop through each gene and create combined heatmaps
genes <- unique(long_data$GeneName)

for (gene in genes) {
  # Filter data for the specific gene
  gene_data <- long_data %>% filter(GeneName == gene)
  
  # Set metric order
  gene_data$Metric <- factor(gene_data$Metric, levels = metric_order)
  
  # Create list to store individual heatmaps
  heatmap_list <- list()
  
  # Generate a separate heatmap for each metric
  for (metric in metric_order) {
    # Filter data for the current metric
    metric_data <- gene_data %>% filter(Metric == metric)
    
    # Pivot data to wide format
    metric_data_wide <- metric_data %>%
      select(Variation_ID, Metric, Value_log) %>%
      spread(key = "Variation_ID", value = "Value_log")
    
    # Create matrix
    metric_matrix <- as.matrix(metric_data_wide[, -1])  # Exclude Metric column
    rownames(metric_matrix) <- metric_data_wide$Metric  # Set row names
    
    # Skip if matrix has NA values
    if (any(is.na(metric_matrix))) {
      warning(paste("Skipped gene:", gene, "metric:", metric, "- Contains NA values"))
      next
    }
    
    # Generate heatmap for the metric (removed 'main' argument)
    p <- pheatmap(metric_matrix,
                  cluster_rows = FALSE,
                  cluster_cols = FALSE,
                  display_numbers = FALSE,
                  color = heatmap_colors,
                  na_col = "grey",
                  legend = TRUE,
                  breaks = seq(min(metric_matrix, na.rm = TRUE), 
                               max(metric_matrix, na.rm = TRUE), 
                               length.out = 101),
                  show_rownames = TRUE,
                  show_colnames = FALSE, # Removed X-axis labels
                  silent = TRUE,         # Silent to prevent auto-plotting
                  height = 4,            # Fixed height for each heatmap
                  width = 6)             # Fixed width for each heatmap
    
    # Store heatmap as a grob
    heatmap_list[[metric]] <- p$gtable
  }
  
  # Combine heatmaps vertically with reduced height and width
  combined_plot <- arrangeGrob(grobs = heatmap_list, ncol = 1)
  
  # Add title at the top of the PDF
  combined_plot_with_title <- arrangeGrob(
    textGrob(paste(gene, "Quality Metrics GnoMAD"), gp = gpar(fontsize = 20, fontface = "bold")),
    combined_plot,
    ncol = 1,
    heights = c(0.1, 1)  # Adjust heights for title and plots
  )
  
  # Save the combined heatmap with title
  heatmap_file <- paste0("heatmap_quality_", gene, ".pdf")
  ggsave(heatmap_file, combined_plot_with_title, width = 15, height = 1 * length(heatmap_list) + 1, units = "in") # Adjusted dimensions
}

# Create summary distribution graphs (no gene categorization)
metrics <- c("SiteQuality", "AS_FS", "AS_MQRankSum", "AS_pab_max", "AS_ReadPosRankSum", "AS_SOR", "AS_VarDP")

for (metric in metrics) {
  # Filter data for the specific metric
  metric_data <- data_clean %>%
    select(Variation_ID, all_of(metric)) %>%
    filter(!is.na(get(metric)))  # Remove rows with NA values in the specific metric
  
  # If the metric is either SiteQuality or AS_VarDP, use the log-transformed values
  if (metric == "SiteQuality" || metric == "AS_VarDP") {
    metric_data$Value_log <- long_data$Value_log[long_data$Metric == metric]
    p <- ggplot(metric_data, aes(x = Value_log)) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +  # Histogram
      geom_density(color = "blue", size = 1) +  # Density curve
      geom_rug(sides = "b", color = "gray", alpha = 0.5) +  # Rug plot at the bottom
      labs(title = paste("Log-transformed Distribution of", metric), x = paste("Log-transformed", metric), y = "Density") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 24, face = "bold"),  # Title size
        axis.title.x = element_text(size = 20),               # X-axis label size
        axis.title.y = element_text(size = 20),               # Y-axis label size
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16),  # X-axis text size
        axis.text.y = element_text(size = 16)                 # Y-axis text size
      )
  } else {
    # For other metrics, plot the raw values
    p <- ggplot(metric_data, aes(x = get(metric))) +
      geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +  # Histogram
      geom_density(color = "blue", size = 1) +  # Density curve
      geom_rug(sides = "b", color = "gray", alpha = 0.5) +  # Rug plot at the bottom
      labs(title = paste("Distribution of", metric), x = metric, y = "Density") +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 24, face = "bold"),  # Title size
        axis.title.x = element_text(size = 20),               # X-axis label size
        axis.title.y = element_text(size = 20),               # Y-axis label size
        axis.text.x = element_text(angle = 90, hjust = 1, size = 16),  # X-axis text size
        axis.text.y = element_text(size = 16)                 # Y-axis text size
      )
  }
  
  # Save plot
  ggsave(paste0("distribution_", metric, ".pdf"), plot = p)
}

