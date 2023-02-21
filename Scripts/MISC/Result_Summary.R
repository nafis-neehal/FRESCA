setwd("/Users/nafisneehal/ESCA")
source("./Scripts/common.R")

######## Data Load
directory <- "./Data/Results/PATT_m6/CC_EC_Vary_Size/matchCCEC/"
summary_dir <- "summary_seed_results/"
filename <- "1500TA_vary_CC_1_5EC_matchCCEC_seed10.csv"

dat<- read.csv(paste(directory, filename, sep = ""))
gt <- read.csv("./Data/Results/PATT_ground.csv")
group2_data <- gt %>% select("PATT_ground")

randomization_ratio <- 1.5
TA_size <- 1500

######## Helper Functions
get_conf_int <- function(data_column){
  mean_value <- mean(data_column)
  n <- length(data_column)
  std <- sd(data_column)
  
  std_error <- std / sqrt(n)
  alpha <- 0.05
  dof <- n-1
  t_score <- qt(p=alpha/2, df=dof, lower.tail = F)
  margin_error <- t_score * std_error
  
  lower_bound <- mean_value - margin_error
  upper_bound <- mean_value + margin_error
  
  return(c(lower_bound, upper_bound))
}

######## Processing
cc_sizes = seq(TA_size, 0, -50)

df <- data.frame(matrix(ncol = 13, nrow = 0))
x <- c("PATT", "VAR_PATT", "P_VAL", "LOWER_CF", "UPPER_CF", "LDM_CC", "VAR_LDM_CC", "LDM_HC", "VAR_LDM_HC", 
       "LDM_IPF_TA", "VAR_LDM_IPF_TA", "LDM_IPF_HC", "VAR_LDM_IPF_HC")
colnames(df) <- x

for (s in cc_sizes){
  #### get slice with particular CC sample size
  temp <- dat %>% filter(CC_Sample_Size==s) %>% select(Mean_Patt, LDM_CC, LDM_HC, LDM_IPF_TA, LDM_IPF_HC)
  
  #### P-VALUE
  group1_data <- temp %>% select(Mean_Patt)
  
  #Fisher F-test to compare variances between two groups
  var_res <- var.test(unlist(group1_data, use.names = FALSE), unlist(group2_data, use.names = FALSE))
  var_p_val <- var_res$p.value
  
  if(var_p_val>=0.05){ #equal variances
    t_test_res <- t.test(group1_data, group2_data, var.equal = TRUE)
  }
  else{ #unequal variances
    t_test_res <- t.test(group1_data, group2_data, var.equal = FALSE)
  }
  
  #### CONFIDENCE INTERVAL
  bounds <- get_conf_int(unlist(group1_data, use.names = FALSE))
  lower_cf <- bounds[1]
  upper_cf <- bounds[2]
  
  #### merge data
  row <- c(mean(temp$Mean_Patt), var(temp$Mean_Patt), t_test_res$p.value, lower_cf, upper_cf,
           mean(temp$LDM_CC), var(temp$LDM_CC), mean(temp$LDM_HC), var(temp$LDM_HC), 
           mean(temp$LDM_IPF_TA), var(temp$LDM_IPF_TA), mean(temp$LDM_IPF_HC), var(temp$LDM_IPF_HC))
  df <- rbind(df, row)
}

colnames(df) <- x
df$CC_Sample_Size <- cc_sizes
df$EC_Sample_Size <- randomization_ratio * TA_size - cc_sizes
df$Diff_From_GT <- abs(df$PATT + 15.174)

df <- df %>% select("CC_Sample_Size", "EC_Sample_Size","PATT", "VAR_PATT", "LOWER_CF", "UPPER_CF", 
                    "Diff_From_GT", "P_VAL", "LDM_CC","VAR_LDM_CC","LDM_HC","VAR_LDM_HC","LDM_IPF_TA", 
                    "VAR_LDM_IPF_TA", "LDM_IPF_HC", "VAR_LDM_IPF_HC")


write.csv(df, paste(directory, paste(summary_dir, filename, sep = ""), sep = ""), row.names = FALSE)
