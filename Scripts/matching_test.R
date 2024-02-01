#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed.csv")
allhat_data$X <- NULL

#create the treated populations
ctrl <- allhat_data %>% filter(RANDASSIGN==2) #Ch
tr1 <- allhat_data %>% filter(RANDASSIGN==3) #Am
tr2 <- allhat_data %>% filter(RANDASSIGN==4) #Li

g1 <- rbind(ctrl, tr1)
g1 <- g1 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))

all_vars <- colnames(allhat_data)
len <- length(all_vars)-2
matching_vars_list <- all_vars[2:len]

All_strings <- lapply(matching_vars_list, as.character)
var_string <- paste(All_strings, collapse = " + ")
outcome <- "RANDASSIGN"
formula <- reformulate(termlabels = var_string, response = outcome)

variables_to_convert <- c("Age_Group", "Gender", "Race_or_Ethnicity", "Education", 
                          "Prior_HypTreat", "Smoker", "Had_MIS", "Had_CRV", "Had_ASCVD", 
                          "Had_MajorST", "Had_T2D", "Had_HDLC", "Had_LVH_ELCT", "Had_LVH_ECHO", 
                          "Has_CHD")
for (var in variables_to_convert) {
  g1[[var]] <- as.factor(g1[[var]])
}

m1 <- matchit(formula, data = g1, method = "nearest", distance = "glm", replace = T)
matched_data <- match.data(m1)




