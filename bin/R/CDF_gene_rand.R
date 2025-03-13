library(ggplot2)

# Load data
data <- read.csv("../../data/gnomad_gene_rand_freq.csv")

# Handle zero allele frequencies to avoid non-finite values
data$allele_freq[data$allele_freq == 0] <- 10^(-10)

# Compute log allele frequency
data$log_allele_freq <- log10(data$allele_freq)

# Identify gene vs. random
data$category <- ifelse(grepl("^Random", data$gene_name), "Random", "Gene")

# Population mapping
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

# Reshape data for population-specific allele frequencies
populations <- names(population_map)
pop_data <- lapply(populations, function(pop) {
  subset <- data[, c("gene_name", paste0(pop, "_af"))]
  colnames(subset)[2] <- "allele_freq"
  subset$population <- population_map[[pop]]
  subset$allele_freq[subset$allele_freq == 0] <- 10^(-10)
  subset$log_allele_freq <- log10(subset$allele_freq)
  subset$category <- ifelse(grepl("^Random", subset$gene_name), "Random", "Gene")
  return(subset)
})
pop_data <- do.call(rbind, pop_data)

# Increase font size
custom_theme <- theme_minimal(base_size = 20)

# Plot overall allele frequency CDF
plot_overall_ecdf <- ggplot(data, aes(x = log_allele_freq, color = category)) +
  stat_ecdf(geom = "step", size = 1) +
  custom_theme +
  labs(title = "CDF of Log(Allele Frequency)",
       x = "Log10(Allele Frequency)",
       y = "Cumulative Probability",
       color = "Category")

# Plot population-specific allele frequency CDF with gene vs. random distinction
plot_pop_ecdf <- ggplot(pop_data, aes(x = log_allele_freq, color = population, linetype = category)) +
  stat_ecdf(geom = "step", size = 1) +
  custom_theme +
  labs(title = "CDF of Log(Allele Frequency) by Population",
       x = "Log10(Allele Frequency)",
       y = "Cumulative Probability",
       color = "Population",
       linetype = "Category")
# Save plots
ggsave("../../results/gene_rand_CDF.pdf", plot = plot_overall_ecdf, width = 15, height = 10)
ggsave("../../results/gene_rand_population_CDF.pdf", plot = plot_pop_ecdf, width = 15, height = 10)

