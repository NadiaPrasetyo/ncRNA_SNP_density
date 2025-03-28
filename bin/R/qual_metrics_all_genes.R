library(tidyverse)
library(pheatmap)

# Load data
data <- read.csv("data/variant_qual_metrics.csv")

# Function to check if a variant is a SNP
is_snp <- function(variation_id) {
  parts <- unlist(strsplit(variation_id, "-"))
  if (length(parts) < 4) return(FALSE)  # Ensure the format is correct
  ref <- parts[length(parts) - 1]  # Second to last part
  alt <- parts[length(parts)]  # Last part
  return(nchar(ref) == 1 & nchar(alt) == 1)  # SNPs have single character ref and alt
}

# Filter for SNPs only
data_filtered <- data %>% 
  filter(sapply(Variation_ID, is_snp)) %>% 
  select(GeneName, Variation_ID, SiteQuality, AS_MQ, AS_FS, AS_MQRankSum, AS_pab_max, 
         AS_ReadPosRankSum, AS_SOR, AS_VarDP)

# Normalize data for heatmap visualization
scaled_data <- data_filtered %>% 
  select(-GeneName, -Variation_ID) %>% 
  scale() %>% 
  as.data.frame()

# Add row names for identification
rownames(scaled_data) <- data_filtered$Variation_ID

# Generate heatmap
pheatmap(scaled_data, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         main = "Variant Quality Metrics Heatmap")
