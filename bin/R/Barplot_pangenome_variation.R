# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggtext)  # For colored text in ggplot

# Step 1: Load the CSV data for both datasets
data1 <- read.csv("../../results/gene_variation_summary.csv")
data2 <- read.csv("../../results/SNP_variation_summary.csv")

# Add a new column to indicate the dataset source
data1$Dataset <- "Dataset1"
data2$Dataset <- "Dataset2"

# Combine the two datasets
combined_data <- bind_rows(data1, data2)

# Step 2: Process the variation types and extract counts
# Split 'VARIATION_TYPES' into separate rows and extract counts and variation types
data_processed <- combined_data %>%
  separate_rows(VARIATION_TYPES, sep = ", ") %>%  # Split by commas
  mutate(
    Variation_Type = str_extract(VARIATION_TYPES, "^[a-zA-Z]+"),  # Extract variation type (e.g., snp, mnp)
    Count = as.numeric(str_extract(VARIATION_TYPES, "(?<=\\()\\d+(?=\\))"))  # Extract the count inside parentheses
  ) %>%
  select(GENE, Variation_Type, Count, Dataset)

# Step 3: Summarize the total count for each variation type per gene and dataset
variation_counts <- data_processed %>%
  group_by(GENE, Variation_Type, Dataset) %>%
  summarise(Total_Variation_Count = sum(Count), .groups = 'drop')

# List of genes with special font color
special_genes <- c(
  "FAM30A", "LINC01671", "lnc-SLCO4A1-8", "MIR4538", "SNAR-A1", "SNAR-B2", 
  "SNAR-C3", "SNAR-C4", "SNAR-G1", "SNAR-G2", "TRE-TTC5-1", "TRG-CCC6-1", 
  "TRG-CCC4-1", "TRV-CAC5-1"
)

# Modify the genes in the special list to use a markdown format for color
variation_counts$GENE <- ifelse(variation_counts$GENE %in% special_genes, 
                                paste0("<span style='color:red;'>", variation_counts$GENE, "</span>"),
                                variation_counts$GENE)

# Step 4: Create the combined stacked bar plot (both datasets stacked)
combined_plot <- ggplot(variation_counts, aes(x = reorder(GENE, -Total_Variation_Count), 
                                              y = Total_Variation_Count, fill = Variation_Type)) +
  geom_bar(stat = "identity", position = "stack") +  # Stacked bars for datasets
  labs(x = "Gene", y = "Total Variation Count (Log scale)", title = "Variation Counts by Type for Each Gene across Datasets") +
  scale_y_continuous(trans = "log10") +  # Set y-axis to log2 scale
  theme_minimal() +
  theme(axis.text.x = ggtext::element_markdown(size = 7, angle = 45, hjust = 1))  # Use markdown to allow coloring

# Save the combined plot
ggsave("../../results/Variation-combined-counts-stacked-bar.pdf", plot = combined_plot, width = 15)

# Step 5: Create separate stacked bar plots for each dataset
# Dataset 1 Plot
data1_plot <- variation_counts %>%
  filter(Dataset == "Dataset1") %>%
  ggplot(aes(x = reorder(GENE, -Total_Variation_Count), y = Total_Variation_Count, fill = Variation_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Gene", y = "Total Variation Count", title = "Variation Counts by Type for Each Gene (Pangenome)") +
  scale_fill_brewer(palette = "YlGnBu") +
  theme_minimal() +
  theme(axis.text.x = ggtext::element_markdown(size = 7, angle = 45, hjust = 1))

# Save Dataset 1 Plot
ggsave("../../results/Variation-counts-pangenome-bar.pdf", plot = data1_plot, width = 15)

# Dataset 2 Plot
data2_plot <- variation_counts %>%
  filter(Dataset == "Dataset2") %>%
  ggplot(aes(x = reorder(GENE, -Total_Variation_Count), y = Total_Variation_Count, fill = Variation_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Gene", y = "Total Variation Count", title = "Variation Counts by Type for Each Gene (SNP155)") +
  scale_fill_brewer(palette = "YlOrRd") +
  theme_minimal() +
  theme(axis.text.x = ggtext::element_markdown(size = 7, angle = 45, hjust = 1))

# Save Dataset 2 Plot
ggsave("../../results/Variation-counts-SNP155-bar.pdf", plot = data2_plot, width = 15)
