# Load the lpSolve library
library(lpSolve)

# Define the inputs to the lpSolve function

# Number of individuals in the trial population
N <- 100

# Number of lowest-level subgroups (combination of age, and gender)
J <- 4

# Propensity score matrix (N x J)
propensity_scores <- matrix(runif(N*J), nrow = N, ncol = J)

# Representation matrix (N x J) where a row represents the sensitive attribute vector of an individual
first_matrix <- diag(J)
representation_matrix <- do.call(rbind, replicate((N/J), first_matrix, simplify=FALSE))
random_rows <- sample(nrow(representation_matrix))
representation_matrix <- representation_matrix[random_rows,]

# Target population's sensitive attribute vector for each subgroup (J x 2)
target_representativeness <- matrix(runif(J*2), nrow = J, ncol = 2)
target_representativeness <- target_representativeness/rowSums(target_representativeness)

# Relaxation parameter for representativeness constraint
epsilon <- 0.2

# Define the decision variable matrix (NxJ) with binary variables
decision_var <- matrix(0, nrow = N, ncol = J)

# Define the objective function coefficients (1 for each propensity score)
objective_coef <- as.vector(t(propensity_scores))

# Define the constraint matrix (including selection, propensity score, and representativeness constraints)
constraint_matrix <- rbind(
  # Selection constraint: Each individual can be selected only once
  #diag(N),
  # Propensity score constraint: The propensity score of each individual should be below a threshold
  propensity_scores,
  # Representativeness constraint: The distribution of sensitive attributes should match the target population
  representation_matrix
)

# Define the constraint right-hand side (including the right-hand side of selection, propensity score, and representativeness constraints)
rhs <- c(rep(0.8, N), (1 + epsilon) * target_representativeness) #rep(1, N), 

# Define the constraint direction (less than or equal to)
direction <- c(rep("<=", N), rep("<=", J)) #rep("<=", N), 

# Define the variable type (binary)
var_types <- rep("binary", N * J)

# Use lpSolve to solve the linear programming problem
lp_result <- lp(direction = "max", objective.in = objective_coef, 
                const.mat = constraint_matrix, const.dir = direction, 
                const.rhs = rhs, all.bin = TRUE) #, var.types = var_types

#Extract the solution from lp_result
selected_decision_var <- matrix(lp_result$solution, nrow=N, ncol = J)

# Extract the objective value from lp_result
objective_value <- lp_result$objval

# Print the selected decision variables and the objective value
cat("Selected Decision Variables:\n")
print(selected_decision_var)
cat("\nObjective Value:", objective_value)
