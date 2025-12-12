setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")

#imports
suppressPackageStartupMessages({
  library(Rlab)
  library(survival)
  library(discreteRV)
  library(dplyr)
})

########## Helper Functions ###########
get_rv <- function(dat,n){
  prob_tab <- table(dat)/n
  names_tab <- names(prob_tab)
  name_list <- c()
  prob_list <- c()
  for (name in names_tab){
    name_list <- c(name_list, name)
    prob_list <- c(prob_list, prob_tab[[name]])
  }
  return(RV(name_list, odds = prob_list))
}

########## Read Data #################
trial_data <- read.csv("./Data/Synthetic/trial_data_ate.csv")
external_data <- read.csv("./Data/Synthetic/external_data_ate.csv")
target_data <- read.csv("./Data/Synthetic/target_data_ate.csv")
hybrid_data <- rbind(trial_data, external_data)

############# SIR Experiment from Trial to Target #############

N <- nrow(trial_data)

available_count_df <- trial_data %>% dplyr::count(x1,x2,x3) %>% arrange(x1,x2,x3)
total_count_df <- data.frame(x1=rep(0:1, each=10),
                             x2=rep(0:1, times=10),
                             x3=rep(0:4, times=4)) %>% 
  arrange(x1,x2,x3)
total_count_df$n <- 0
total_count_df <- total_count_df %>% 
  left_join(available_count_df, by=c("x1","x2","x3")) %>% 
  mutate(n=n.y) %>% 
  select(x1,x2,x3,n) %>% 
  replace(is.na(.),0)
total_count_df$proba <- total_count_df$n/nrow(trial_data)

trial_data$weights <- 1

for (i in 1:nrow(total_count_df)){
  prob <- total_count_df[i,'proba']
  cx1 <- total_count_df[i,'x1']
  cx2 <- total_count_df[i,'x2']
  cx3 <- total_count_df[i,'x3']
  trial_data <- trial_data %>% mutate(weights = ifelse(x1 == cx1,
                                                       ifelse(x2 == cx2, 
                                                              ifelse(x3 == cx3, 1/prob, weights), weights), weights))
}

trial_data$weights <- trial_data$weights/sum(trial_data$weights)

#new probabilities in g
trial_data$new_proba <- 0

t1 <- get_rv(target_data$x1, nrow(target_data))
t2 <- get_rv(target_data$x2, nrow(target_data))
t3 <- get_rv(target_data$x3, nrow(target_data))
tjoint <- joint(t1, joint(t2, t3))

for (i in 1:length(outcomes(tjoint))){
  outcome <- outcomes(tjoint)[i]
  prob <- probs(tjoint)[i]
  out_list <- str_split_1(outcome, ',')
  trial_data <- trial_data %>% mutate(new_proba = ifelse(x1 == out_list[1],
                                                         ifelse(x2 == out_list[2], 
                                                                ifelse(x3 == out_list[3], prob, new_proba), new_proba), new_proba))
}

trial_data$new_weights <- trial_data$weights * trial_data$new_proba
trial_data$new_weights <- trial_data$new_weights/sum(trial_data$new_weights)

trial_model <- lm(y~Tr+x1+x2+x3, data = trial_data, weights = trial_data$new_weights)
summary(trial_model)

# tr_weights <- trial_data$new_weights[trial_data$Tr==1]
# cn_weights <- trial_data$new_weights[trial_data$Tr==0]
# tr_w_mean <- weighted.mean(trial_data$y[trial_data$Tr==1], tr_weights)
# cn_w_mean <- weighted.mean(trial_data$y[trial_data$Tr==0], cn_weights)

# adj_trialATE <- abs(tr_w_mean - cn_w_mean)
# paste("Adjusted Trial ATE is: ", adj_trialATE)


############# SIR Experiment from Hybrid to Target #############
available_count_df <- hybrid_data %>% dplyr::count(x1,x2,x3) %>% arrange(x1,x2,x3)
total_count_df <- data.frame(x1=rep(0:1, each=10),
                             x2=rep(0:1, times=10),
                             x3=rep(0:4, times=4)) %>% 
  arrange(x1,x2,x3)
total_count_df$n <- 0
total_count_df <- total_count_df %>% 
  left_join(available_count_df, by=c("x1","x2","x3")) %>% 
  mutate(n=n.y) %>% 
  select(x1,x2,x3,n) %>% 
  replace(is.na(.),0)
total_count_df$proba <- total_count_df$n/nrow(hybrid_data)

hybrid_data$weights <- 1 


for (i in 1:nrow(total_count_df)){
  prob <- total_count_df[i,'proba']
  cx1 <- total_count_df[i,'x1']
  cx2 <- total_count_df[i,'x2']
  cx3 <- total_count_df[i,'x3']
  hybrid_data <- hybrid_data %>% mutate(weights = ifelse(x1 == cx1,
                                                         ifelse(x2 == cx2, 
                                                                ifelse(x3 == cx3, 1/prob, weights), weights), weights))
}

hybrid_data$weights <- hybrid_data$weights/sum(hybrid_data$weights)

#new probabilities in g
hybrid_data$new_proba <- 0

t1 <- get_rv(target_data$x1, nrow(target_data))
t2 <- get_rv(target_data$x2, nrow(target_data))
t3 <- get_rv(target_data$x3, nrow(target_data))
tjoint <- joint(t1, joint(t2, t3))

for (i in 1:length(outcomes(tjoint))){
  outcome <- outcomes(tjoint)[i]
  prob <- probs(tjoint)[i]
  out_list <- str_split_1(outcome, ',')
  hybrid_data <- hybrid_data %>% mutate(new_proba = ifelse(x1 == out_list[1],
                                                           ifelse(x2 == out_list[2], 
                                                                  ifelse(x3 == out_list[3], prob, new_proba), new_proba), new_proba))
}

hybrid_data$new_weights <- hybrid_data$weights * hybrid_data$new_proba
hybrid_data$new_weights <- hybrid_data$new_weights/sum(hybrid_data$new_weights)

hybrid_model <- lm(y~Tr+x1+x2+x3+x4, data = hybrid_data, weights = hybrid_data$new_weights)
summary(hybrid_model)

# tr_weights <- hybrid_data$new_weights[hybrid_data$Tr==1]
# cn_weights <- hybrid_data$new_weights[hybrid_data$Tr==0]
# tr_w_mean <- weighted.mean(hybrid_data$y[hybrid_data$Tr==1], tr_weights)
# cn_w_mean <- weighted.mean(hybrid_data$y[hybrid_data$Tr==0], cn_weights)
# adj_hybridATE <- abs(tr_w_mean - cn_w_mean)
# paste("Adjusted Hybrid ATE is: ", adj_hybridATE)

new_hybrid_data <- sample_n(hybrid_data, nrow(hybrid_data), replace = T,
                            hybrid_data$new_weights)
new_hybrid_model <- lm(y~Tr+x1+x2+x3+x4, data = new_hybrid_data, weights = new_hybrid_data$new_weights)
summary(new_hybrid_model)
# tr_weights <- new_hybrid_data$new_weights[new_hybrid_data$Tr==1]
# cn_weights <- new_hybrid_data$new_weights[new_hybrid_data$Tr==0]
# tr_w_mean <- weighted.mean(new_hybrid_data$y[new_hybrid_data$Tr==1], tr_weights)
# cn_w_mean <- weighted.mean(new_hybrid_data$y[new_hybrid_data$Tr==0], cn_weights)
# adj_newhybridATE <- abs(tr_w_mean - cn_w_mean)
# paste("Adjusted Hybrid ATE is: ", adj_newhybridATE)



