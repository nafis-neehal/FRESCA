# suppressPackageStartupMessages({
#   library(ipfr)
#   library(dplyr)
# })
# 
# data <- read.csv("./App/th_controls.csv")
# data <- na.omit(data)
# 
# #change hierarchy of Races (Asian, Black, White, Others):(0,1,2,3)
# data$RACE_[is.element(data$RACE_, c(2,3,4,5,6))] <- 3
# data$RACE_[data$RACE_==7] <- 2
# 
# targets <- list()
# targets$GENDER <- tibble(
#   '1' = .50,
#   '2' = .50
# )
# targets$RACE_ <- tibble(
#   '1' = .25,
#   '2' = .25,
#   '0' = .25,
#   '3' = .75
# )

return_ipf <- function(data, targets, var_list, num_samples){
  set.seed(42)
  data$weight <- 1
  dat2 <- data %>% select(unlist(var_list), weight)
  dat2 <- as_tibble(dat2)
  result <- ipu(dat2, targets, max_iterations = 1000)
  sampled_data <- data[sample(seq_len(nrow(data)), replace = TRUE, num_samples, prob = result$weight_tbl$weight),]
  drop_column <- c("weight")
  sampled_data <- sampled_data[, !(names(sampled_data) %in% drop_column)]
  return (sampled_data)
}

# sampled_data <- return_ipf(data, targets, num_samples = 1965)
# rowSums(table(sampled_data%>%select(GENDER, RACE_)))/1965
# colSums(table(sampled_data%>%select(GENDER, RACE_)))/1965
