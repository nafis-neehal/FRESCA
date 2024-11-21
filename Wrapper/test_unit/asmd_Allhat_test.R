#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

allhat_data <- read.csv("./Data/ALLHAT_data/ALLHAT_processed_secondary.csv")
allhat_data$X <- NULL

#create the treated populations
ctrl <- allhat_data %>% filter(RANDASSIGN==2) #Ch
tr1 <- allhat_data %>% filter(RANDASSIGN==3) #Am
tr2 <- allhat_data %>% filter(RANDASSIGN==4) #Li


## Split the Trial population into TA, CC and EC
## Change RANDASSIGN to 0 for controls and 1 for treated
g1 <- rbind(tr1, ctrl)
g1 <- g1 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))
g2 <- rbind(tr2, ctrl)
g2 <- g2 %>% mutate(RANDASSIGN = ifelse(RANDASSIGN==2, 0, 1))

get_asmd_dataframe <- function(TA, CC, EC_world, Biased_EC_world){
  
  TrialPop <- rbind(TA, CC)
  TrialPop$inTrial <- "Yes" 
  
  ExtPop <- EC_world
  ExtPop$inTrial <- "No"
  
  BExtPop <- Biased_EC_world
  BExtPop$inTrial <- "No"
  
  Pop <- rbind(TrialPop, ExtPop)
  BPop <- rbind(TrialPop, BExtPop)
  
  vars <- c("Age_Group", "Gender", "Race_or_Ethnicity", "Education",
            "Smoker", "SBP", "FPG", "TC", "CVDHISTORY", "CVDPOINTS", "SERUMCREAT",
            "GFRESTIMATE")
  
  unbiasedTab <- CreateTableOne(vars = vars, strata = "inTrial", data = Pop, test=FALSE)
  biasedTab <- CreateTableOne(vars = vars, strata = "inTrial", data = BPop, test=FALSE)
  
  dataPlot <- data.frame(variable   = rownames(ExtractSmd(unbiasedTab)),
                         Unbiased   = as.numeric(ExtractSmd(unbiasedTab)),
                         Biased     = as.numeric(ExtractSmd(biasedTab)))
  
  return (dataPlot)
  
}



TA_Size <- 4000
CC_Size <- 2000
SC_Size <- TA_Size - CC_Size
IPF_maxiter <- 100
all_data <- partition_data_generic(g1, partition_seed=1)
target_data <- data.frame(all_data[2])
RCT <- data.frame(all_data[[1]][[1]])
BIASED_EC <- data.frame(all_data[[1]][[2]])

vars <- colnames(tr1 %>% select(Age_Group:BMI))
tr1$group <- "c1"
BIASED_EC$group <- "c2"


save_ASMD <- function(pop1, pop2, seed, controls, var_list){
  pop1$group <- "c1"
  pop2$group <- "c2"
  v <- CreateTableOne(var_list, strata = "group", data = rbind(pop1, pop2), test=F)
  df <- data.frame(Seed        = seed,
                   Controls    = controls,
                   Variable    = rownames(ExtractSmd(v)),
                   SMD         = as.numeric(ExtractSmd(v)))
}




