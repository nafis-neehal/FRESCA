setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
#seed <- (123, 10, 42) works perfectly 
#seed <- (1,10,42) works perfectly

#imports
suppressPackageStartupMessages({
  library(Rlab)
  library(survival)
  library(discreteRV)
  library(coxed)
})

########## Helper Functions ###########
get_data_joint_marginal <- function(sample_size, seed_value){
  #create joint probability sample based on target marginals
  #assuming x1, x2 and x3 are independent
  x1_marg <- c(0.446, 0.554)
  x2_marg <- c(0.312, 0.688)
  x3_marg <- c(0.693, 0.120, 0.039, 0.100, 0.048)
  
  df <- data.frame(x1=rep(0:1, each=10),
                   x2=rep(0:1, times=10),
                   x3=rep(0:4, times=4)) %>% arrange(x1,x2,x3)
  
  df$probs <- 0
  
  for (i in 1:nrow(df)){
    df[i,4] = x1_marg[df[i,1]+1] * x2_marg[df[i,2]+1] * x3_marg[df[i,3]+1]
  }
  
  set.seed(seed_value)
  data <- sample_n(df, size=sample_size, replace=T, weight=df$probs)
  
  return(data)
}

############################ Trial Data ############################
N <- 10000 #total samples
D <- 5 #features
seed <- 5
set.seed(seed)

#create x1 ~ Bernouli(0.5) ~ GENDER ~ 0/Male 1/Female
x1 <- rbern(N, prob = 0.1)

#create x2 ~ Bernouli(0.6) ~ Age ~ 0/40-59 1/59+
x2 <- rbern(N, prob = 0.6)

#create x3 ~ Categorical(0,1,2,3,4 ~ each 20%) ~ Race
#~ 0/NH-White 1/NH-Black 2/NH-Asian 3/Hispanic 4/Other
x3 <- sample(seq(0,4), N, replace=TRUE, prob=c(0.40,0.23,0.13,0.15,0.09))

#create x4 ~ Normal(60,5)-60 ~ HDL Cholesterol
x4 <- rnorm(n=N, mean=60, sd=5)-60

#create x5 ~ Normal(21,2)-21 ~ BMI
x5 <- rnorm(n=N, mean=21, sd=2)-21

#create Tr ~ Bernouli(t) ~ t={0.67(2:1), 0.75(3:1)} ~ Treatment indicator
Tr <- rbern(N, 0.5) #2:1 randomization ratio

#create Failure time variable T_failure ~ Exp[log(a).Tr + log(b).X]
#a = {0.5, 0.75, 0.875, 1},
#b={(1.25, 0.67, 0.55, 0.98, 1.06),(2.25, 0.4, 0.55, 0.93, 1.21)}

X <- matrix(do.call("cbind",list(x1,x2,x3,x4,x5)), nrow=N, ncol=D)
conf_deg <- c(1.25, 0.67, 0.55, 0.98, 1.06)
#T_failure <- c(exp(log(0.5)*(Tr) + log(conf_deg) %*% t(X))) #avoid confounding in RCT
T_failure <- c(exp((log(3.5)*x1+log(0.5)*(1-x1))*(Tr) + log(conf_deg) %*% t(X)))

# set.seed(seed)
# T_failure <- exp(log(0.5)*(1-Tr)+rnorm(N,4,1)-4)

#create censor time variable T_censor ~ Exp(0.1)
T_censor <- rexp(N, 0.1)

#outcome variable - Survival Time Y = min(T_failure, T_censor)
y <- pmin(T_failure, T_censor)

#Status variable S = 1 if (T_failure_i < T_censor_i), else S = 0
failed <- ifelse(T_failure<T_censor, 1, 0)

#create the final dataframe
trial_data <- data.frame(x1, x2, x3, x4, x5, Tr, y, failed)

#introduce bias in a lowest level subgroup
#trial_data$x3[trial_data$x1=='0' & trial_data$x2=='1' & trial_data$x3=='4']=0

#Hazard Ratio
res <- coxph(Surv(y, failed)~Tr, data = trial_data)
summary(res)


############################ External Control Data ############################
N <- 10000 #total samples
D <- 5 #features
seed <- 20
set.seed(seed)

#create x1 ~ Bernouli(0.5) ~ GENDER ~ 0/Male 1/Female
x1 <- rbern(N, prob = 0.55)

#create x2 ~ Bernouli(0.6) ~ Age ~ 0/40-59 1/59+
x2 <- rbern(N, prob = 0.4)

#create x3 ~ Categorical(0,1,2,3,4 ~ each 20%) ~ Race
#~ 0/NH-White 1/NH-Black 2/NH-Asian 3/Hispanic 4/Other
x3 <- sample(seq(0,4), N, replace=TRUE, prob=c(0.45, 0.19, 0.08, 0.21, 0.07))

#create x4 ~ Normal(60,5)-60 ~ HDL Cholesterol
x4 <- rnorm(n=N, mean=60, sd=10)-60

#create x5 ~ Normal(21,2)-21 ~ BMI
x5 <- rnorm(n=N, mean=23, sd=2)-21

#create Tr ~ Bernouli(t) ~ t={0.67(2:1), 0.75(3:1)} ~ Treatment indicator
Tr <- rbern(N, 0) #All Tr value is 0

#create Failure time variable T_failure ~ Exp[log(a).Tr + log(b).X]
#a = {0.5, 0.75, 0.875, 1}, 
#b={(1.25, 0.67, 0.55, 0.98, 1.06),(2.25, 0.4, 0.55, 0.93, 1.21)}
conf_deg <- c(1.25, 0.67, 0.55, 0.98, 1.06)
X <- matrix(do.call("cbind",list(x1,x2,x3,x4,x5)), nrow=N, ncol=D)
T_failure <- c(exp(log(0.5)*(Tr) + log(conf_deg) %*% t(X))) #avoid confounding in RCT

#T_failure <- exp(log(0.5)*(1-Tr)+rnorm(N,4,1)-4)

#create censor time variable T_censor ~ Exp(0.1)
T_censor <- rexp(N, 0.25)

#outcome variable - Survival Time Y = min(T_failure, T_censor)
y <- pmin(T_failure, T_censor)

#Status variable S = 1 if (T_failure_i < T_censor_i), else S = 0
failed <- ifelse(T_failure<T_censor, 1, 0)

#create the final dataframe
external_data <- data.frame(x1, x2, x3, x4, x5, Tr, y, failed)

#No treatment effect in external data source, because there is no Tr=1

############################ Target Data ############################
N <- 20000 #total samples
D <- 5 #features
seed <- 55
set.seed(seed)

#create x1, x2 and x3 based on pre-specified marginals from ~ NHANES
# temp_data <- get_data_joint_marginal(N, seed)
# x1 <- temp_data$x1
# x2 <- temp_data$x2
# x3 <- temp_data$x3

#create x1 ~ Bernouli(0.5) ~ GENDER ~ 0/Male 1/Female
x1 <- rbern(N, prob = 0.9)

#create x2 ~ Bernouli(0.6) ~ Age ~ 0/40-59 1/59+
x2 <- rbern(N, prob = 0.6)

#create x3 ~ Categorical(0,1,2,3,4 ~ each 20%) ~ Race
#~ 0/NH-White 1/NH-Black 2/NH-Asian 3/Hispanic 4/Other
x3 <- sample(seq(0,4), N, replace=TRUE, prob=c(0.40,0.23,0.13,0.15,0.09))

#create x4 ~ Normal(60,5)-60 ~ HDL Cholesterol
x4 <- rnorm(n=N, mean=45, sd=8)-45

#create x5 ~ Normal(21,2)-21 ~ BMI
x5 <- rnorm(n=N, mean=16, sd=4)-16

#create Tr ~ Bernouli(t) ~ t={0.50} ~ Treatment indicator
Tr <- rbern(N, 0.5) #1:1 randomization ratio

#create Failure time variable T_failure ~ Exp[log(a).Tr + log(b).X]
#a = {0.5, 0.75, 0.875, 1},
#b = {(1.25, 0.67, 0.55, 0.98, 1.06),(2.25, 0.4, 0.55, 0.93, 1.21)}
conf_deg <- c(2.25, 0.4, 0.32, 0.93, 1.21)
X <- matrix(do.call("cbind",list(x1,x2,x3,x4,x5)), nrow=N, ncol=D)
#T_failure <- c(exp(log(3.5)*x1*(Tr) + log(conf_deg) %*% t(X))) 
T_failure <- c(exp((log(3.5)*x1+log(0.5)*(1-x1))*(Tr) + log(conf_deg) %*% t(X))) 

#T_failure <- exp(log(0.875)*(1-Tr)+rnorm(N,4,1)-4)

#create censor time variable T_censor ~ Exp(0.1)
T_censor <- rexp(N, 0.15)

#outcome variable - Survival Time Y = min(T_failure, T_censor)
y <- pmin(T_failure, T_censor)

#Status variable S = 1 if (T_failure_i < T_censor_i), else S = 0
failed <- ifelse(T_failure<T_censor, 1, 0)

#create the final dataframe
target_data <- data.frame(x1, x2, x3, x4, x5, Tr, y, failed)

#Hazard Ratio
target_res <- coxph(Surv(y, failed)~Tr, data = target_data)
summary(target_res)

##################### Save the data files #################
write.csv(trial_data, file="./Data/Synthetic/trial_data_paper.csv", row.names = F)
write.csv(external_data, file="./Data/Synthetic/external_data_paper.csv", row.names = F)
write.csv(target_data, file="./Data/Synthetic/target_data_paper.csv", row.names = F)
###########################

