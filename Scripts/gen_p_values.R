setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

confidence_interval <- function(vector, interval) {
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

directory <- "./Data/Results/M6/"
patt_filename <- "PATT_Summary_extrabias.csv"
ldm_filename <- "LDM_Summary_extrabias.csv"
Final_PATT_Summary <- read.csv(paste(directory, patt_filename, sep = ""))
Final_LDM_Summary <- read.csv(paste(directory, ldm_filename, sep = ""))

Final_PATT_Summary <- Final_PATT_Summary %>% filter(CC_Size==500)
Final_LDM_Summary <- Final_LDM_Summary %>% filter(CC_Size==500)

Final_LDM_Summary_agg <- Final_LDM_Summary %>%
  group_by(Seed) %>%
  summarise_at(vars(LDM_TA, LDM_CC, LDM_HC_Prop, LDM_HC_IPF, LDM_TA_Prop_IPF, 
                    LDM_HC_Prop_IPF), funs(median))

patt_g <- read.csv("./Data/Results/PATT_ground50.csv")
satt_g <- read.csv("./Data/Results/SATT_ground50.csv")

pldm_g <- read.csv("./Data/Results/PLDM_ground50.csv")
sldm_g <- read.csv("./Data/Results/SLDM_ground50.csv")

pldm_g_agg <- pldm_g %>% group_by(Seed) %>%
  summarise_at(vars(LDM_CC_IPF), funs(median))

sldm_g_agg <- sldm_g %>% group_by(Seed) %>%
  summarise_at(vars(LDM_CC), funs(median))

#compare between PATT Reference and PATT values
res <- t.test(Final_LDM_Summary_agg$LDM_HC_IPF, pldm_g_agg$LDM_CC_IPF)
res$p.value

# #compare between PATT Reference and PATT values
# res <- t.test(Final_PATT_Summary$TA_HC_Both, patt_g$PATT_ground)
# res$p.value





