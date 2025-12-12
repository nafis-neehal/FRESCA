# Install and load required packages
# install.packages("forestplot")
# install.packages("survival")
library(forestplot)
library(survival)

# Example data (replace this with your actual data)
# Suppose you have a data frame named "my_data" with variables: "study", "hazard_ratio", "lower_ci", "upper_ci"
my_data <- data.frame(
  study = c("Study 1", "Study 2", "Study 3", "Study 4"),
  hazard_ratio = c(0.75, 0.60, 1.20, 1.50),
  lower_ci = c(0.60, 0.45, 0.95, 1.30),
  upper_ci = c(0.90, 0.80, 1.50, 1.70)
)

# Prepare the data in the correct format for forestplot
forest_data <- data.frame(
  Variables = my_data$study,
  HR = my_data$hazard_ratio,
  lower = my_data$lower_ci,
  upper = my_data$upper_ci
)

# Create the forest plot with labels
forestplot(
  mean = forest_data$HR,
  lower = forest_data$lower,
  upper = forest_data$upper,
  labeltext = forest_data$Variables,
  txt_gp = fpTxtGp(label = gpar(cex = 0.8)),
  xlab = "Hazard Ratio (95% CI)",
  is.summary = c(TRUE, rep(FALSE, nrow(forest_data) - 1)),
  xticks = c(0.2, 0.5, 1, 1.5, 2),
  col = fpColors(box = "blue", lines = "black", summary = "red"),
  boxsize = 0.3,
  cex = 0.8,
  zero = 1,
  hrzl_lines = list("2" = gpar(lwd = 1, lty = 2))
)
