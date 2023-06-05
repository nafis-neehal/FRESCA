seed <- 123
set.seed(seed) # Set seed for reproducibility

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

############## Trial ###############
N <- 10000 #total samples
D <- 5 #features
seed <- 5
set.seed(seed)

#create x1 ~ Bernouli(0.5) ~ GENDER ~ 0/Male 1/Female
x1 <- rbern(N, prob = 0.5)

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
Tr <- rbern(N, 0.67) #2:1 randomization ratio

# Generate continuous outcome variable using covariates and treatment
# unconfounded because Tr is random, and confounding is balanced in both arm, 
# confounding nullifies each other
# treatment effect ~ 4
#y <- rnorm(N, mean = ifelse(Tr == 0, 2*x1 - 3*x2 + 1*x3, 2*x1 - 3*x2 + 1*x3 + 4*Tr), sd = 2)
y <- 2*x1 - 3*x2 + 1*x3 + 4*Tr

# Create a data frame with the generated data
trial_data <- data.frame(x1,x2,x3,x4,x5, Tr, y)

#introduce bias in a first level subgroup
#trial_data$x1[trial_data$x1=='1']=0

#introduce bias in a second level subgroup
trial_data$x2[trial_data$x1=='0' & trial_data$x2=='1']=0

#introduce bias in lowest-level subgroups
# trial_data$x3[trial_data$x1=='0' & trial_data$x2=='1' & trial_data$x3=='4']=0
# trial_data$x3[trial_data$x1=='0' & trial_data$x2=='0' & trial_data$x3=='4']=0
# trial_data$x3[trial_data$x1=='1' & trial_data$x2=='0' & trial_data$x3=='2']=1

# trialATE <- abs(mean(trial_data$y[trial_data$Tr==1])-mean(trial_data$y[trial_data$Tr==0]))
# paste("Trial ATE is: ", trialATE)
trial_model <- lm(y~Tr+x1+x2+x3, data = trial_data)
summary(trial_model)

############## External ###############
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

# Generate continuous outcome variable using covariates and treatment
# unconfounded because Tr is random, and confounding is balanced in both arm, 
# confounding nullifies each other
# treatment effect ~ NULL
# y <- rnorm(N, mean = ifelse(Tr == 0, 2*x1 - 3*x2 + 1*x3, 
#                             2*x1 - 3*x2 + 1*x3 + 4*Tr), sd = 2)
y <- 2*x1 - 3*x2 + 1*x3 + 4*Tr

# Create a data frame with the generated data
external_data <- data.frame(x1,x2,x3,x4,x5, Tr, y)

############## Target ###############
N <- 100000 #total samples
D <- 5 #features
seed <- 55
set.seed(seed)

#create x1, x2 and x3 based on pre-specified marginals from ~ NHANES
temp_data <- get_data_joint_marginal(N, seed)
x1 <- temp_data$x1
x2 <- temp_data$x2
x3 <- temp_data$x3

#create x4 ~ Normal(60,5)-60 ~ HDL Cholesterol
x4 <- rnorm(n=N, mean=45, sd=8)-45

#create x5 ~ Normal(21,2)-21 ~ BMI
x5 <- rnorm(n=N, mean=16, sd=4)-16

#create Tr ~ Bernouli(t) ~ t={0.50} ~ Treatment indicator
Tr <- rbern(N, 0.5) #1:1 randomization ratio

# Generate continuous outcome variable using covariates and treatment
# unconfounded because Tr is random, and confounding is balanced in both arm, 
# confounding nullifies each other
# treatment effect ~ 12
# y <- rnorm(N, mean = ifelse(Tr == 0, 2*x1 - 3*x2 + 1*x3, 
#                             2*x1 - 3*x2 + 1*x3 + 10*Tr), sd = 2)
y <- 2*x1 - 3*x2 + 1*x3 + 10*Tr

# Create a data frame with the generated data
target_data <- data.frame(x1,x2,x3,x4,x5, Tr, y)

# targetATE <- abs(mean(target_data$y[target_data$Tr==1])-mean(target_data$y[target_data$Tr==0]))
#paste("Target ATE is: ", targetATE)
target_model <- lm(y~Tr+x1+x2+x3, data = target_data)
summary(target_model)

##################### Save the data files #################
write.csv(trial_data, file="./Data/Synthetic/trial_data_ate.csv", row.names = F)
write.csv(external_data, file="./Data/Synthetic/external_data_ate.csv", row.names = F)
write.csv(target_data, file="./Data/Synthetic/target_data_ate.csv", row.names = F)



