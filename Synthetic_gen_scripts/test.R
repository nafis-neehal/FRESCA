setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")

#imports
suppressPackageStartupMessages({
  library(Rlab)
  library(generalize)
})

set.seed(123)

#size
trial_size <- 10000
target_size <- 20000

#covariate generation
trial_X1 <- rbern(trial_size, 0.4)
target_X1 <- rbern(target_size, 0.6)

trial_X2 <- rbern(trial_size, 0.7)
target_X2 <- rbern(target_size, 0.2)

#treatment generation
trial_T <- rbern(trial_size, 0.5)
target_T <- rbern(target_size, 0.5)

#outcome generation
trial_Y <- 5 * trial_X1 * trial_X2 + trial_T + rnorm(trial_size, mean = 0, sd=0.25)
target_Y <- 5 * target_X1 * target_X2 + target_T + rnorm(target_size, mean = 0, sd=0.25)

#dataframe
trial_data <- data.frame(X1 = trial_X1, X2 = trial_X2, Tr = trial_T, Y = trial_Y)
target_data <- data.frame(X1 = target_X1, X2 = target_X2, Tr = target_T, Y = target_Y)

#treatment effect ATE
trial_effect <- mean(trial_data$Y[trial_data$Tr==1]) - mean(trial_data$Y[trial_data$Tr==0])
target_effect <- mean(target_data$Y[target_data$Tr==1]) - mean(target_data$Y[target_data$Tr==0])

#treatment effect Linear Model
trial_fit <- lm(Y~Tr, data=trial_data)
summary(trial_fit)
confint(trial_fit, "Tr", 0.95)

target_fit <- lm(Y~Tr, data=target_data)
summary(target_fit)
confint(target_fit, "Tr", 0.95)

########## Use Generalize Package #########
#trial membership indicator variable
trial_data$trial <- 1
target_data$trial <- 0

target_data$Y <- NA
target_data$Tr <- NA

covariates <- c("X1", "X2")

#generalize effect using trial to target
generalize_object <- generalize(outcome="Y", treatment="Tr", trial="trial",
                                selection_covariates = covariates,
                                data = rbind(trial_data, target_data),
                                method = "weighting",
                                selection_method = "lr",
                                trimpop = TRUE)
summary(generalize_object)



