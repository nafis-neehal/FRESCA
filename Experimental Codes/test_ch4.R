# Load the lpSolve package
library(lpSolve)

set.seed(123)

N <- 20
PS <- seq(1,20,1)
S <- 5

obj_in <- PS 
const_mat <- matrix(data = rep(1,N), nrow = 1, byrow = T)
const_dir <- c("=")
const_rhs <- c(S)

res <- lp("max", objective.in = obj_in, const.mat = const_mat, 
   const.dir = const_dir, all.bin = T, const.rhs = const_rhs)

res$solution
