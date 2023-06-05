# Load the lpSolve package
library(lpSolve)

set.seed(123)

# Dummy data for trial population and target population
n <- 10000 # number of trial patients
t <- 200 # number of target patients
n_cat <- 8 #number of low-level categories

# Propensity scores for trial population
propensity_scores <- runif(n)

# Distribution of lowest-level subgroups in target population
p_target <- c(0.1, 0.1, 0.15, 0.15, 0.1, 0.05, 0.1, 0.25)
#p_target <- c(0.94, 0.01, 0.01, 0.01, 0.01, 0.01)

# Allowed margin of error for each subgroup
d <- 0.05

# Formulate the linear programming problem
# Objective function: maximize the total propensity score
f.obj <- propensity_scores #1xN

# Inequality constraints: distribution matching with margin of error
f.const.lhs_ratio <- matrix(0, nrow = n_cat, ncol = n) #JxN
for (i in 1:n) {
  # Extract covariate information for patient i
  # Here, we assume patients are randomly assigned to subgroups
  covariates <- sample(1:n_cat, 1, replace = TRUE)
  f.const.lhs_ratio[covariates, i] <- 1/n
}
f.const.rhs_ratio_high <- p_target + d
f.const.rhs_ratio_low  <- p_target - d
f.const.dir_ratio_high <- rep("<=", n_cat)
f.const.dir_ratio_low  <- rep(">=", n_cat)

# # Equality constraint: summation of weights equals 1
# f.const.lhs_eq <- matrix(0, nrow = 1, ncol = n) #1xN
# f.const.lhs_eq[1, ] <- 1
# f.const.rhs_eq <- 1
# f.const.dir_eq <- rep("=",1)

# Bounds on weights: non-negativity and upper-bound
f.const.lhs_nonneg <- diag(n) #NxN
f.const.rhs_nonneg <- rep(0, n)
f.const.dir_nonneg <- rep(">=",n)
# f.const.rhs_upper <- rep(1, n)
# f.const.dir_upper <- rep("<=",n)

f.const.mat <- do.call("rbind", list(f.const.lhs_ratio, f.const.lhs_ratio, 
                                     #f.const.lhs_eq,
                                     f.const.lhs_nonneg))#, f.const.lhs_nonneg))
f.const.dir <- do.call("c", list(f.const.dir_ratio_high, f.const.dir_ratio_low,
                                     #f.const.dir_eq,
                                     f.const.dir_nonneg))#, f.const.dir_upper))
f.const.rhs <- as.numeric(c(f.const.rhs_ratio_high, f.const.rhs_ratio_low,
                                 #f.const.rhs_eq,
                                 f.const.rhs_nonneg))#, f.const.rhs_upper))

#bounds <- list(lower = rep(0, n), upper = rep(1, n))

# Solve the linear programming problem
lp_result <- lp("max", f.obj, f.const.mat, f.const.dir, f.const.rhs)

# Extract the optimal weights
weights <- lp_result$solution

# Print the optimal weights
cat("Optimal weights for trial patients:", weights)

#test if we are getting only strata weights for each strata
indices <- which(weights!=0)
weight_vec <- weights[weights!=0]

mat <- t(f.const.lhs_ratio)
weight_mat <- mat[indices, ]

catg <- data.frame(x1=rep(0:1, each=n_cat/2),
                   x2=rep(0:1, each=n_cat/4),
                   x3=rep(0:1, times=n_cat/2)) %>% arrange(x1,x2,x3)
catg$weight <- 0

for(i in 1:nrow(weight_mat)){
  weight_val <- weight_vec[i]
  cat_comb_index <- which(weight_mat[i,]!=0)
  catg[cat_comb_index,]$weight <- weight_val
}

trial_data <- data.frame(x1=rbern(n, 0.5),
                        x2=rbern(n,0.5),
                        x3=rbern(n,0.5),
                        ps=propensity_scores)

trial_data$w <- 0

for (i in 1:nrow(catg)){
  trial_data$w[(trial_data$x1==catg[i,]$x1) & 
               (trial_data$x2==catg[i,]$x2) &
               (trial_data$x3==catg[i,]$x3)] <- catg[i,]$weight
}

new_trial_data <- sample_n(trial_data, n, replace = T, 
                           weight = trial_data$w)

res <- table(new_trial_data$x1, new_trial_data$x2, new_trial_data$x3)/n

new_trial_data$adj_w <- new_trial_data$w/sum(new_trial_data$w)
new_trial_data$w_ps <- new_trial_data$adj_w * new_trial_data$ps



