#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed.csv")
allhat_data$X <- NULL


############### Influence Test ################
# Perform Cox proportional hazards regression
cox_model <- coxph(Surv(EVENTDAYS, EVENT) ~ ., data = allhat_data)

# Get the summary of the model to see coefficients, p-values, etc.
summary_table <- summary(cox_model)

# Extract the p-values to identify the significance of each input feature
p_values <- summary_table$coefficients[, "Pr(>|z|)"]

# Rank the input features based on their significance (lower p-value indicates higher significance)
significant_features <- names(p_values)[order(p_values)]

# Print the ranked list of significant features with corresponding p-values
significant_features_with_p_values <- data.frame(Feature = significant_features, P_Value = p_values[order(p_values)])

# Round up the entries in the 'column_name' column to 3 decimal points
significant_features_with_p_values <- significant_features_with_p_values %>%
  mutate(P_Value = round(P_Value, digits = 4))
print(significant_features_with_p_values)
#################################################

############## See Distribution of Smoke, Had_MajorST and Had_HDLC ##############

#create the treated populations
ctrl <- allhat_data %>% filter(RANDASSIGN==2) %>% mutate(grp = "ctrl") #Ch
tr1 <- allhat_data %>% filter(RANDASSIGN==3) %>% mutate(grp = "tr1") #Am
tr2 <- allhat_data %>% filter(RANDASSIGN==4) %>% mutate(grp = "tr2") #Li

fullpop <- do.call("rbind", list(ctrl, tr1, tr2))

table1(~ Gender + Age_Group + Race_or_Ethnicity + Smoker + Had_MajorST + Had_HDLC  | grp, data = fullpop, overall=F)








