  library(ggplot2)
  library(reshape2)
  
  # Read data
  data <- read.csv("../../results/var_freq.csv")
  
  # Define base transition matrix
  bases <- c("A", "C", "G", "T")
  mutation_matrix <- matrix(NA, nrow = 4, ncol = 4, dimnames = list(bases, bases))
  
  # Fill matrix with mutation counts
  total_mutations <- sum(data$Mutation_Count)
  for (i in 1:nrow(data)) {
    from_base <- substr(data$Mutation[i], 1, 1)
    to_base <- substr(data$Mutation[i], 4, 4)
    mutation_matrix[from_base, to_base] <- data$Mutation_Count[i] / total_mutations
  }
  
  # Convert to log odds
  expected_freq <- 1/12  # 0.083333... as expected baseline
  df <- as.data.frame(as.table(mutation_matrix))
  df$Freq[is.na(df$Freq)] <- 0  # Replace NA with 0
  
  df$LogOdds <- log2(df$Freq / expected_freq)
  
  df$Var1 <- factor(df$Var1, levels = bases)
  df$Var2 <- factor(df$Var2, levels = bases)
  
  # Plot heatmap
  plot<-ggplot(df, aes(Var2, Var1, fill = LogOdds)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = "Mutation Log Odds Heatmap", x = "To Base", y = "From Base", fill = "Log2 Odds") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"))
  
  ggsave("../../results/var_freq.pdf", plot = plot, height=5, width = 6)