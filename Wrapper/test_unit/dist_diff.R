#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

calculate_categorical_difference <- function(population1, population2) {
  # Identify categorical columns
  categorical_cols <- sapply(population1, function(x) is.factor(x) || is.character(x))
  
  # Initialize the total statistic value
  total_statistic <- 0
  total_p_sig <- 0
  
  # Loop through each categorical column
  for (col in names(population1)[categorical_cols]) {
 
    obs_freq_table <- table(population1[[col]], population2[[col]])
    
    # Perform Chi-squared test and sum up the statistic values
    #chi_squared_result <- chisq.test(table(population1[[col]], population2[[col]]))
    chi_squared_result <- chisq.test(table(population1[[col]]), table(population2[[col]]))

    cat(col, "Categorical Stat:", chi_squared_result$statistic, "\n")
    cat(col, "P val", chi_squared_result$p.value,"\n")
    
    total_statistic <- total_statistic + chi_squared_result$statistic
    if (chi_squared_result$p.value < 0.05) {
      cat("Sig in categorical: ", col, "\n")
      total_p_sig <- total_p_sig + 1
    }
  }
  
  return(c(total_statistic, total_p_sig))
}

calculate_numeric_difference <- function(population1, population2) {
  # Identify numeric columns
  numeric_cols <- sapply(population1, is.numeric)
  
  # Initialize the total statistic value
  total_statistic <- 0
  total_p_sig <- 0
  
  # Loop through each numeric column
  for (col in names(population1)[numeric_cols]) {
    ks_result <- ks.test(population1[[col]], population2[[col]])
    cat(col, "Numeric Stat: ", ks_result$statistic, "\n")
    total_statistic <- total_statistic + ks_result$statistic
    if (ks_result$p.value < 0.05) {
      cat("Sig in numeric:", col, "\n")
      total_p_sig <- total_p_sig + 1
    }
  }
  
  return(c(total_statistic, total_p_sig))
}


calculate_distributional_difference <- function(population1, population2) {
  
  # Calculate the distributional difference for categorical variables 
  chi_squared_res <- calculate_categorical_difference(population1, population2)
  chi_squared_categorical <- chi_squared_res[1]
  chi_squared_pcount <- chi_squared_res[2]
  
  # Calculate the distributional difference for numeric variables 
  ks_res <- calculate_numeric_difference(population1, population2)
  ks_numeric <- ks_res[1]
  ks_numeric_pcount <- ks_res[2]
  
  # Combine the differences into a single interpretable value
  total_difference <- ks_numeric + chi_squared_categorical 
  total_p_count <- ks_numeric_pcount + chi_squared_pcount
  
  return(c(total_difference, total_p_count))
}

allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed_secondary.csv")
allhat_data$X <- NULL

#create the treated populations
ctrl <- allhat_data %>% filter(RANDASSIGN==2) #Ch
tr1 <- allhat_data %>% filter(RANDASSIGN==3) #Am
tr2 <- allhat_data %>% filter(RANDASSIGN==4) #Li

## Split the Trial population into TA, CC and EC
## Change RANDASSIGN to 0 for controls and 1 for treated
g1 <- rbind(tr1, ctrl)
g1 <- g1 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))
g2 <- rbind(tr2, ctrl)
g2 <- g2 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))


TA_Size <- 4000
CC_Size <- 2000
SC_Size <- TA_Size - CC_Size
IPF_maxiter <- 100
all_data <- partition_data_generic(g1, partition_seed=1)
target_data <- data.frame(all_data[2])
RCT <- data.frame(all_data[[1]][[1]])
BIASED_EC <- data.frame(all_data[[1]][[2]])

new_pop <- tr1 %>% sample_n(size = nrow(BIASED_EC), replace = T)
res1 <- calculate_distributional_difference(new_pop %>% select(Age_Group:BMI), 
                                            BIASED_EC %>% select(Age_Group:BMI))



print(res1)
