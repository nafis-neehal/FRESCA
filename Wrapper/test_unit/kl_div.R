#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

library(philentropy)

calculate_total_kl_divergence <- function(population1, population2) {
  # Calculate KL divergence for each variable and sum them up
  total_kl_divergence <- 0
  
  for (col_name in colnames(population1)) {
    variable1 <- population1[[col_name]]
    variable2 <- population2[[col_name]]
    
    if (is.numeric(variable1) && is.numeric(variable2)) {
      # Convert numerical variables to probability distributions
      ndist1 <- variable1 / sum(variable1)
      ndist2 <- variable2 / sum(variable2)
      
      kl_divergence <- KL(rbind(ndist1, ndist2), unit = "log2")
      cat(col_name, kl_divergence, "\n")
      
      total_kl_divergence <- total_kl_divergence + kl_divergence
    } 
    
    else if (is.factor(variable1) || is.character(variable1) && is.factor(variable2) || is.character(variable2)) {
      # Convert categorical variables to probability distributions
      
      cdist1 <- table(variable1)
      cdist2 <- table(variable2)
      
      print(cdist1)
      print(cdist2)
      
      kl_divergence <- KL(rbind(cdist1, cdist2), est.prob = "empirical")
      cat(col_name, kl_divergence, "\n")
      total_kl_divergence <- total_kl_divergence + kl_divergence
    }
  }
  
  return(total_kl_divergence)
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
result <- calculate_total_kl_divergence(new_pop %>% select(Age_Group:BMI), 
                                        BIASED_EC %>% select(Age_Group:BMI))

print(result)

# # Example usage:
# # Create example populations with numeric and categorical variables
# population1 <- data.frame(
#   numeric_var = c(1.2, 2.3, 3.4, 4.5),
#   categorical_var = factor(c("A", "B", "A", "C"))
# )
# 
# population2 <- data.frame(
#   numeric_var = c(1.5, 2.0, 3.8, 4.0),
#   categorical_var = factor(c("B", "B", "A", "C"))
# )
# 
# # Calculate and print the total KL divergence
# result <- calculate_total_kl_divergence(population1, population2)
# print(result)
