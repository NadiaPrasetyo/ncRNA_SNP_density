# Load required libraries
library(dplyr)
library(ggplot2)
library(stats)  # For KS test

# Load CHH and CpG data
data_chh <- read.csv("../../data/CHH_methylation_data_temp.csv")
data_cpg <- read.csv("../../data/CpG_methylation_data.csv")

data_chh$Context <- "CHH"
data_cpg$Context <- "CpG"

# Combine both datasets
all_data <- bind_rows(data_chh, data_cpg)

# Convert Methylation_Percentage to numeric
all_data$Methylation_Percentage <- as.numeric(all_data$Methylation_Percentage)

# Check structure of the data
print(str(all_data))

# Create a dataframe to store KS test results
ks_results_chh <- data.frame(Tissue = character(), KS_Statistic = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)
ks_results_cpg <- data.frame(Tissue = character(), KS_Statistic = numeric(), P_Value = numeric(), stringsAsFactors = FALSE)

# Part 1: Tissue-based KS tests
# Get the unique tissue types
tissues <- unique(all_data$Tissue)

for (tissue in tissues) {
  tissue_data <- filter(all_data, Tissue == tissue)
  tissue_data <- tissue_data %>%
    mutate(GeneType = ifelse(grepl("^Random\\d+", GeneName), "Random", "Gene"))
  
  for (context in c("CHH", "CpG")) {
    context_data <- filter(tissue_data, Context == context)
    
    gene_methylation <- context_data$Methylation_Percentage[context_data$GeneType == "Gene"]
    random_methylation <- context_data$Methylation_Percentage[context_data$GeneType == "Random"]
    
    if (length(gene_methylation) > 0 && length(random_methylation) > 0) {
      ks_test <- ks.test(gene_methylation, random_methylation)
      
      if (context == "CHH") {
        ks_results_chh <- rbind(ks_results_chh, data.frame(Tissue = tissue, KS_Statistic = ks_test$statistic, P_Value = ks_test$p.value))
      } else {
        ks_results_cpg <- rbind(ks_results_cpg, data.frame(Tissue = tissue, KS_Statistic = ks_test$statistic, P_Value = ks_test$p.value))
      }
    }
  }
}

write.csv(ks_results_chh, "../../results/DNAme/_CHH_ks_test_results.csv", row.names = FALSE)
write.csv(ks_results_cpg, "../../results/DNAme/_CpG_ks_test_results.csv", row.names = FALSE)

# Part 2: Tissue-based plotting
for (tissue in tissues) {
  tissue_data <- filter(all_data, Tissue == tissue)
  tissue_data <- tissue_data %>%
    mutate(GeneType = ifelse(grepl("^Random\\d+", GeneName), "Random", "Gene")) %>%
    mutate(GeneType = paste(Context, GeneType))
  
  summary_stats <- tissue_data %>% group_by(GeneType) %>% summarise(
    Median = median(Methylation_Percentage, na.rm = TRUE),
    Mean = mean(Methylation_Percentage, na.rm = TRUE)
  )
  
  plot <- ggplot(tissue_data, aes(x = GeneType, y = Methylation_Percentage, color = Context)) +
    geom_jitter(alpha = 0.7, size = 3, width = 0.1) +
    geom_text(data = summary_stats, aes(x = GeneType, y = Mean, label = round(Mean, 2)), vjust = -1, color = "red") +
    labs(
      title = paste("Methylation Percentage in", tissue, "(CHH & CpG)"),
      x = "Gene Type",
      y = "Methylation Percentage"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0("../../results/DNAme/", tissue, "_CpG_CHH_scatter_jitter_plot.pdf"), plot = plot, width = 5, height = 8)
}

# Part 3: Individual gene analysis
filter_randoms <- function(gene, tissue_data) {
  randoms <- tissue_data %>%
    filter(
      grepl("^Random\\d+", GeneName),
      abs(Length - gene$Length) <= 200,
      abs(Median_CG_Content - gene$Median_CG_Content) <= 20
    )
  return(randoms)
}

# Get unique gene names excluding "Random#" genes
genes <- unique(all_data$GeneName[!grepl("^Random\\d+", all_data$GeneName)])

for (gene in genes) {
  gene_data <- filter(all_data, GeneName == gene)
  
  all_gene_data <- data.frame(Tissue = character(), Methylation_Percentage = numeric(), GeneType = character(), Context = character())
  
  for (tissue in tissues) {
    tissue_data <- filter(all_data, Tissue == tissue)
    
    for (context in c("CHH", "CpG")) {
      context_data <- filter(tissue_data, Context == context)
      randoms <- filter_randoms(gene_data[1, ], context_data)
      
      all_gene_data <- bind_rows(
        all_gene_data,
        context_data %>% filter(GeneName == gene) %>% mutate(GeneType = paste(tissue, context)),
        randoms %>% mutate(GeneType = paste(tissue, context, "random"))
      )
    }
  }
  
  if (nrow(all_gene_data) > 0) {
    summary_stats <- all_gene_data %>% group_by(GeneType) %>% summarise(
      Median = median(Methylation_Percentage, na.rm = TRUE),
      Mean = mean(Methylation_Percentage, na.rm = TRUE)
    )
    
    plot <- ggplot(all_gene_data, aes(x = GeneType, y = Methylation_Percentage, color = GeneType)) +
      geom_jitter(width = 0.2, alpha = 0.7) +
      geom_text(data = summary_stats, aes(x = GeneType, y = Mean, label = round(Mean, 2)), vjust = -1, color = "red") +
      labs(
        title = paste(gene, "CHH & CpG vs Random"),
        x = "Gene Type",
        y = "Methylation Percentage"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(paste0("../../results/DNAme/", gene, "_CHH_CpG_individual_scatter_plot.pdf"), plot = plot, width = 8, height = 8)
  }
}
