library(survival)

set.seed(123) # set random seed for reproducibility

# generate 1000 observations
n <- 1000

# generate 3 covariates (X1, X2, X3) as normally distributed random variables
X1 <- rnorm(n)
X2 <- rnorm(n)
X3 <- rnorm(n)

# set the true coefficients for the covariates
beta1 <- 0.5
beta2 <- -0.3
beta3 <- 0.2

# generate the linear predictor using the true coefficients and the covariate values
lp <- beta1*X1 + beta2*X2 + beta3*X3

# set the baseline survival function using the Weibull distribution
alpha <- 1
lambda <- 0.01
S0 <- exp(-lambda * (0:n)^(1/alpha))

# generate the survival times using the 'rweibull' function with right censoring
T <- rweibull(n, shape=alpha, scale=lambda)
C <- ifelse(runif(n) < 0.2, 1, 0)
T[C==1] <- runif(sum(C==1), max(T), max(T)+50)

# generate a binary treatment indicator variable
Z <- ifelse(runif(n) < 0.5, 1, 0)

# set the treatment effect on the linear predictor
delta <- -0.5

# update the linear predictor based on the treatment effect
lp_treat <- lp + delta*Z

# generate the survival times for the treated group
T_treat <- rweibull(n, shape=alpha, scale=lambda)
C_treat <- ifelse(runif(n) < 0.2, 1, 0)
T_treat[C_treat==1] <- runif(sum(C_treat==1), max(T_treat), max(T_treat)+50)

# create a data frame with the survival times, covariates, treatment indicator, and censoring indicator
data <- data.frame(T=ifelse(Z==1, T_treat, T), X1=X1, X2=X2, X3=X3, Z=Z, C=ifelse(Z==1, C_treat, C))

# create a survival object with the survival time, censoring, and treatment information
surv_obj <- with(data, Surv(T, C, type="right") ~ Z)

# fit a Cox proportional hazards model with treatment as a covariate
coxph_obj <- coxph(Surv(T,C, type = 'right') ~ X1 + X2 + X3 + Z, data=data)

# print the summary of the Cox model
summary(coxph_obj)
