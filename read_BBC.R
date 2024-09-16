library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)


remove_suffix <- function(string, suffix) {
  sub(paste0(suffix, "$"), "", string)
}

transform_lower_triangular <- function(df) {
  
  elements <- c()

  # Loop through each row and column to extract lower triangular elements
  for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
      if (!is.na(df[i, j]))
      {elements[length(elements) + 1] <- df[i, j]}
    }
  }
  
  # Combine the elements into a single-column data frame
  result <- data.frame(value = elements)
  
  return(result)
}

compute_row_differences <- function(df) {
  # Compute differences between consecutive rows for each column
  differences <- lapply(df, function(col) c(NA, diff(col)))
  
  # Combine differences into a dataframe
  differences_df <- as.data.frame(differences)
  
  # Rename columns
  colnames(differences_df) <- paste0("diff_", colnames(df))
  
  return(differences_df[-1,])
  # return(differences_df[-1,-ncol(differences_df)])
}

extract_numeric <- function(column_names, string) {
     # Extract numeric part using regular expression
     numeric_part <- gsub(string, "", column_names)
     # Convert to numeric
     numeric_values <- as.numeric(numeric_part)
     return(numeric_values)
}

# Set the working directory

setwd("C:/Users/Federico/Desktop/BayesClustering/final")

# mean en

file_list <- list.files(pattern = "mean_entropy.*\\.csv$")


data_list <- list()

# Loop through each file and read the second column
for (file in file_list) {
  # Read the CSV file
  data <- read.csv2(file)
  # Keep only the second column
  second_column <- data[, 2, drop = FALSE]
  # Add the second column to the list
  data_list[[file]] <- second_column
}

# Combine all second columns into one data frame
combined_data <- do.call(cbind, data_list)

# Apply the function to each element in the list
cleaned_list <- sapply(file_list, remove_suffix, suffix = ".csv")

colnames(combined_data) <- cleaned_list
row.names(combined_data) <- 2:(nrow(combined_data)+1)

final_mean <- combined_data

# numeric part of column names
numeric_values <- extract_numeric(colnames(final_mean), "mean_entropy_var")

colnames(final_mean) <- numeric_values

# Reorder columns based on numeric values
final_mean <- final_mean[, order(numeric_values)]

write.csv2(combined_data, "final_mean.csv", row.names = FALSE)

#### rel en

file_list <- list.files(pattern = "rel_entropy_pbbs_var.*\\.csv$")

data_list <- list()

# Loop through each file and flatten
for (file in file_list) {
  # Read and transform
  data <- transform_lower_triangular(read.csv2(file))
  # Add to the list
  data_list[[file]] <- data
}

combined_data <- do.call(cbind, data_list)

cleaned_list <- sapply(file_list, remove_suffix, suffix = ".csv")

colnames(combined_data) <- cleaned_list

# Extract numeric part of column names
numeric_values <- extract_numeric(colnames(combined_data), "rel_entropy_pbbs_var")

colnames(combined_data) <- numeric_values

# Reorder dataframe columns based on numeric values
combined_data <- combined_data[, order(numeric_values)]

final_rel <- combined_data

# Create an empty vector
result_vector <- c()

# Generate the vector for n = 2 to 6
for (n in 2:6) {
  # Calculate the number of times to repeat 'n'
  repetitions <- n * (n - 1) / 2
  
  # Append 'n' to the result vector 'repetitions' times
  result_vector <- c(result_vector, rep(n, times = repetitions))
}

final_rel$index <- result_vector

write.csv2(final_rel, "final_rel.csv", row.names = FALSE)

maxes_final_rel <- final_rel %>%
  group_by(index) %>%
  summarise(across(everything(), max))

##### visual

df <- final_mean %>%
  mutate(index = row_number()+1)

df_long <- df %>%
  pivot_longer(cols = -index, names_to = "s", values_to = "value")

# Create the line plot
p <- ggplot(df_long, aes(x = index, y = value, color = s)) +
  geom_line() +
  labs(title = "",
       x = "K",
       y = expression(bar(S))) +
  theme_bw() +
  theme(text = element_text(size = 20))

p

ggsave("plot_m1_k.png", plot = p, width = 10, height = 6)

##### visual

# plot m1 against various s

df3 <- as.data.frame(t(final_mean))

df4 <- df3 %>%
  mutate(index = as.numeric(row.names(df3)))

df_long <- df4 %>%
  pivot_longer(cols = -index, names_to = "variable", values_to = "value")

# Create the line plot
p <- ggplot(df_long, aes(x = index, y = value, color = variable)) +
  geom_line() +
  labs(title = "Line Plot of m1 for each s, different k",
       x = "s of the variance",
       y = "Mean entropy") +
  theme_bw()

ggsave("plot_m1_s.png", plot = p, width = 10, height = 6)

##### visual

df_long_maxes <- maxes_final_rel %>%
  pivot_longer(cols = -index, names_to = "s", values_to = "value")

# Create the line plot
p <- ggplot(df_long_maxes, aes(x = index, y = value, color = s)) +
  geom_line() +
  labs(title = "",
       x = "K",
       y = expression(bar(S)[l][,][m])) +
  theme_bw() +  
theme(text = element_text(size = 20))

p


ggsave("plot_m2_max.png", plot = p, width = 10, height = 6)

setwd("C:/Users/Federico/Desktop/BayesClustering")