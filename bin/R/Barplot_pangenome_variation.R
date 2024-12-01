# Load necessary libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(ggtext)  # For colored text in ggplot

# Step 1: Load the CSV data
data <- read.csv("../../results/gene_variation_summary.csv")

# Step 2: Process the variation types and extract counts
# Split 'VARIATION_TYPES' into separate rows and extract counts and variation types
data_processed <- data %>%
  separate_rows(VARIATION_TYPES, sep = ", ") %>%  # Split by commas
  mutate(
    Variation_Type = str_extract(VARIATION_TYPES, "^[a-zA-Z]+"),  # Extract variation type (e.g., snp, mnp)
    Count = as.numeric(str_extract(VARIATION_TYPES, "(?<=\\()\\d+(?=\\))"))  # Extract the count inside parentheses
  ) %>%
  select(GENE, Variation_Type, Count)

# Step 3: Summarize the total count for each variation type per gene
variation_counts <- data_processed %>%
  group_by(GENE, Variation_Type) %>%
  summarise(Total_Variation_Count = sum(Count), .groups = 'drop')


# List of genes with special font color
special_genes <- c(
  "FAM30A", "LINC01671", "lnc-SLCO4A1-8", "MIR4538", "SNAR-A1", "SNAR-B2", 
  "SNAR-C3", "SNAR-C4", "SNAR-G1", "SNAR-G2", "TRE-TTC5-1", "TRG-CCC6-1", 
  "TRG-CCC4-1", "TRV-CAC5-1"
)


# Step 4: Modify the genes in the special list to use a markdown format for color
variation_counts$GENE <- ifelse(variation_counts$GENE %in% special_genes, 
                                paste0("<span style='color:red;'>", variation_counts$GENE, "</span>"),
                                variation_counts$GENE)

# Step 5: Create the bar plot with variation counts labeled on bars
plot <- ggplot(variation_counts, aes(x = reorder(GENE, -Total_Variation_Count), y = Total_Variation_Count, fill = Variation_Type)) +
  geom_bar(stat = "identity", position = "stack") +
  #geom_text(aes(label = Total_Variation_Count), position = position_stack(vjust = 0.5), color = "black") +  # Add labels to bars
  labs(x = "Gene", y = "Total Variation Count", title = "Variation Counts by Type for Each Gene") +
  scale_fill_brewer(palette = "YlGnBu") +  # Optional: Choose color palette for fill
  theme_minimal() +
  theme(axis.text.x = ggtext::element_markdown(size = 7, angle = 45, hjust = 1))  # Use markdown to allow coloring


#Step 6: save the plot
ggsave("../../results/Variation-counts-pangenome-bar.pdf", plot = plot, width = 15)