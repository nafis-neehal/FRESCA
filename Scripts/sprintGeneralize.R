setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")

confidence_interval <- function(vector, interval){
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

###### Create New Target Population by extrapolation using joint probabilities ######

summarized_nhanes <- nhanes %>% 
  group_by(Gender, Age_Group, Race_or_Ethnicity) %>% 
  summarise(n = sum(background_n))

summarized_nhanes$ratio <- summarized_nhanes$n / sum(summarized_nhanes$n)

sum_nhanes <- summarized_nhanes %>% select(Gender, Age_Group, Race_or_Ethnicity, 
                                           ratio) %>% ungroup()

target_size <- 20000

new_nhanes <- sample_n(sum_nhanes, size = target_size, replace = T, 
                       weight = ratio)



new_nhanes %>% group_by(Gender, Age_Group, Race_or_Ethnicity) %>% 
  summarise(n=n()/target_size)

###### Load SPRINT data and split into Trial and External ######

#collect the weights for generalizing
new_nhanes$EVENTDAYS <- NA
new_nhanes$EVENT <- NA

new_nhanes$inTrial <- 0
data$inTrial <- 1

te_list <- c()

for (i in 1:50){
  
  print(i)
  
  set.seed(i)
  
  #new_data <- sample_frac(data, 0.95)
  new_data <- sample_n(data, size= nrow(data), replace = T)
  
  # treat <- new_data %>% filter(RANDASSIGN==1)
  # cont <- new_data %>% filter(RANDASSIGN==0)
  # 
  # new_treat <- sample_n(treat, 2000)
  # new_cont <- sample_n(cont, 1500)
  # new_new_data <- rbind(new_treat, new_cont)

  merged <- rbind(new_data %>% select(Gender, Age_Group, Race_or_Ethnicity, EVENTDAYS, EVENT, inTrial),
                new_nhanes %>% select(Gender, Age_Group, Race_or_Ethnicity, EVENTDAYS, EVENT, inTrial))

  merged <- merged %>% mutate(Gender = recode(Gender, "Male"=0,"Female"=1),
                            Age_Group = recode(Age_Group, "40-59"=0, "59+"=1),
                            Race_or_Ethnicity = recode(Race_or_Ethnicity,
                                                       'NH White' = 0,
                                                       'NH Black' = 1,
                                                       'NH Asian' = 2,
                                                       'Hispanic' = 3,
                                                       'Other' = 4))


  covariates <- c("Age_Group", "Gender", "Race_or_Ethnicity")

  assess_object <- assess(trial="inTrial", selection_covariates = covariates,
                        data = merged, selection_method = "lr",
                        trimpop = F)

  #summary(assess_object)

  all_weights <- assess_object$weights
  weights <- all_weights[1:nrow(new_data)]

  iptw_model <- coxph(Surv(EVENTDAYS, EVENT)~RANDASSIGN, data = new_data, weights = weights)
  result <- summary(iptw_model)

  te <- exp(result$coefficients[1])
  
  te_list <- c(te_list, te)
  
}

print("------------------------------------------------")
cat(paste0("Estimated Mean is ",mean(unlist(te_list)), "\n",
           "Variance is ", var(unlist(te_list)), "\n",
           "Standard Deviation is ", sd(unlist(te_list)), "\n",
           "Confidence Interval is", confidence_interval(unlist(te_list), 0.95)))

# TA_size <- 1500
# TA_World <- data %>% filter(RANDASSIGN==1)
# TA <- sample_n(TA_World, TA_size)
# CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)
# EC_world <- external_population %>% filter(RANDASSIGN==0)


















