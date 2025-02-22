# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)  # Load viridis for color palettes

# Read the CSV file
data <- read.csv("../../data/snp_frequencies.csv")  # Replace with your CSV file path

# Convert the data to long format for plotting
long_data <- data %>%
  pivot_longer(cols = c("SNP155_Count", "Pangenome_Count"), 
               names_to = "CountType", 
               values_to = "Count") %>%
  mutate(Mutation = as.factor(Mutation))  # Ensure Mutation is a factor

# Summarize the data to get the total counts for each gene and mutation type
total_data <- long_data %>%
  group_by(Gene, Mutation) %>%
  summarize(Total_Count = sum(Count), .groups = 'drop')

# Histogram for SNP155 Count
snp155_plot <- ggplot(data, aes(x = Gene, fill = Mutation)) +
  geom_bar(aes(y = SNP155_Count), stat = "identity") +
  theme_minimal() +
  labs(title = "SNP155 Count by Gene and Mutation",
       x = "Gene",
       y = "SNP155 Count") +
  scale_fill_viridis(discrete = TRUE) +  # Use viridis colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the SNP155 Count plot
ggsave("../../results/snp155_count_frequency.pdf", plot = snp155_plot, width = 20, height = 7)

# Histogram for Pangenome Count
pangenome_plot <- ggplot(data, aes(x = Gene, fill = Mutation)) +
  geom_bar(aes(y = Pangenome_Count), stat = "identity") +
  theme_minimal() +
  labs(title = "Pangenome Count by Gene and Mutation",
       x = "Gene",
       y = "Pangenome Count") +
  scale_fill_viridis(discrete = TRUE) +  # Use viridis colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the Pangenome Count plot
ggsave("../../results/pangenome_count_frequency.pdf", plot = pangenome_plot, width = 20, height = 7)

# Combined histogram for total counts of SNP155 and Pangenome
combined_plot <- ggplot(total_data, aes(x = Gene, y = Total_Count, fill = Mutation)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Combined SNP frequency for SNP155 and Pangenome Count",
       x = "Gene",
       y = "Total Count") +
  scale_fill_viridis(discrete = TRUE) +  # Use viridis colors
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the combined histogram plot
ggsave("../../results/combined_frequency.pdf", plot = combined_plot, width = 20, height = 7)
