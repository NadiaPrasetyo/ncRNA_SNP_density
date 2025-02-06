# Load necessary libraries
library(dplyr)
library(tidyr)
library(readr)

# Read the data
data <- read_csv("../../data/methylation_signals.csv")

# Function to clean methylation signal data by removing NaNs, "NA" strings, and spaces
clean_methylation_data <- function(methylation_signal_str) {
  # Remove the brackets and split the string by commas
  methylation_values <- gsub("[\\[\\]]", "", methylation_signal_str) %>%
    strsplit(",") %>%
    unlist() %>%
    trimws()  # Remove any leading or trailing whitespace
  
  # Remove entries that are "NA" (as a string) or empty
  methylation_values <- methylation_values[!(tolower(methylation_values) == "na" | methylation_values == "")]
  
  # Strip spaces within the values (e.g., " 9" becomes "9")
  methylation_values <- gsub("\\s", "", methylation_values)
  
  # Convert the cleaned data to numeric, while suppressing warnings and ensuring only valid numbers are retained
  methylation_values <- suppressWarnings(as.numeric(methylation_values))
  
  # Change all Nan to 0
  methylation_values[is.nan(methylation_values)] <- 0
  
  # Remove the NA values from the start and end of the array
  methylation_values <- methylation_values[!is.na(methylation_values)]
  
  return(methylation_values)
}

# Apply the cleaning function
data$MethylationSignal <- lapply(data$MethylationSignal, clean_methylation_data)

# Check if the data is cleaned properly
head(data$MethylationSignal, 10)  # View first few cleaned signals

# Reshape the data for comparison
data <- data %>%
  mutate(TissueType = gsub(".*_(.*)\\.WGBS.*", "\\1", BigWigFile)) %>%
  select(GeneName, TissueType, MethylationSignal)


