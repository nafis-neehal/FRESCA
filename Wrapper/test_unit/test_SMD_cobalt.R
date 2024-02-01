library(cobalt)

calculate_smd_cobalt <- function(population1, population2) { #first control then treated
  data_combined <- rbind(population1, population2)
  
  bal_table <- bal.tab(data_combined, treat = data_combined$treat, s.d.denom = "treated")
  
  smd_data <- data.frame(Variable = rownames(bal_table$smd), SMD = bal_table$smd$std.diff)
  return(smd_data)
}

# Example data
population1 <- data.frame(
  Gender = c("Male", "Female", "Male", "Male", "Female"),
  Age = c(25, 30, 28, 35, 40),
  Height = c(165, 170, 175, 160, 180),
  treat = c(0,0,0,0,0)
)

population2 <- data.frame(
  Gender = c("Female", "Male", "Female", "Male", "Female"),
  Age = c(22, 29, 26, 34, 38),
  Height = c(160, 168, 172, 158, 175),
  treat = c(1,1,1,1,1)
)

cmd <- rbind(population1, population2)

res <- bal.tab(cmd, treat = cmd$treat, s.d.denom = "treat", binary="std")

# Calculate standardized mean differences using cobalt
#smd_results <- calculate_smd_cobalt(population1, population2)



# Print the results
#print(smd_results)
