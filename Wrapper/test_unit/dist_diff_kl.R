calculate_categorical_difference <- function(population1, population2) {
  # Identify categorical columns
  categorical_cols <- sapply(population1, function(x) is.factor(x) || is.character(x))
  
  # Initialize the total statistic value
  total_statistic <- 0
  total_p_sig <- 0
  
  # Loop through each categorical column
  for (col in names(population1)[categorical_cols]) {
    table1 <- table(population1[[col]])
    table2 <- table(population2[[col]])
    
    # Perform Chi-squared test and sum up the statistic values
    chi_squared_result <- chisq.test(table1, table2)
    total_statistic <- total_statistic + chi_squared_result$statistic
    print(chi_squared_result$p.value)
    if (chi_squared_result$p.value < 0.05) {
      total_p_sig <- total_p_sig + 1
    }
  }
  
  return(c(total_statistic, total_p_sig))
}

# Example usage
population1 <- data.frame(
  numeric_var1 = c(1.5, 2.7, 3.2, 4.8),
  categorical_var1 = c("A", "B", "A", "C"),
  categorical_var2 = c("X", "Y", "X", "Z")
)

population2 <- data.frame(
  numeric_var1 = c(1.7, 2.5, 3.0, 4.5),
  categorical_var1 = c("B", "A", "B", "C"),
  categorical_var2 = c("X", "Y", "Z", "X")
)

difference_score <- calculate_categorical_difference(population1, population2)
print(difference_score)
