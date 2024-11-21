encode_categorical <- function(input_df) {
  # Create a copy of the input dataframe
  encoded_df <- input_df
  
  # Iterate through columns
  for (col in colnames(encoded_df)) {
    col_data <- encoded_df[[col]]
    
    # Check if the column is categorical (contains text data)
    if (is.character(col_data) && length(unique(col_data)) < length(col_data)) {
      # Create a mapping of unique categories to numeric values
      unique_categories <- unique(col_data)
      category_mapping <- as.numeric(factor(unique_categories))
      
      # Replace text categories with numeric values
      encoded_df[[col]] <- category_mapping[match(col_data, unique_categories)]
    }
  }
  
  return(encoded_df)
}

data <- data.frame(
  numeric_var = c(10, 20, 30, 40),
  binary_var = c("A", "B", "A", "A"),
  categorical_var = c("A", "B", "A", "C")
)

# Call the encoding function
encoded_data <- encode_categorical(data)

# Display the encoded 
print(encoded_data)


