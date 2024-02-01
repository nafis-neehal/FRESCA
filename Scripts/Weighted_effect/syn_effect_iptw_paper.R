#Reproducing Paper: Generalizing Treatment Effect Estimates from Sample to Population: 
#Stuart et. al. 2016

setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")

# Load necessary packages
suppressPackageStartupMessages({
  library(ipw)
  library(tableone)
  library(ggplot2)
  library(ggfortify)
  library(survminer)
})

########## Read Data #################
trial_data <- read.csv("./Data/Synthetic/trial_data_paper.csv")
external_data <- read.csv("./Data/Synthetic/external_data_paper.csv")
target_data <- read.csv("./Data/Synthetic/target_data_paper.csv")
#hybrid_data <- rbind(trial_data, external_data)


########### Kaplan Meier Curves for Trial/Target Difference ##########
new_trial <- trial_data
new_trial$trial <- "Trial"

new_target <- target_data
new_target$trial <- "Target"

fit <- survfit(Surv(y, failed) ~ trial, data=rbind(new_trial, new_target))
ggsurvplot(
  fit,
  data = rbind(new_trial, new_target),
  size = 1,                 # change line size
  palette =
    c("#E7B800", "#2E9FDF"),# custom color palettes
  conf.int = TRUE,          # Add confidence interval
  pval = TRUE,              # Add p-value
  risk.table = TRUE,        # Add risk table
  risk.table.col = "strata",# Risk table color by groups
  legend.labs =
    c("Target", "Trial"),    # Change legend labels
  risk.table.height = 0.25, # Useful to change when you have multiple groups
  ggtheme = theme_bw()      # Change ggplot2 theme
)

#difference
survdiff(Surv(y, failed) ~ trial, data=rbind(new_trial, new_target))

#######################

# Create a table of covariate balance before weighting to check balance in Tr/Cn
table1 <- CreateTableOne(vars = c("x1", "x2", "x3", "x4", "x5"), 
                         strata = "Tr", data = trial_data)

# Check the covariate balance before weighting
print(table1, test = FALSE)

########## Apply on Trial Data ###########

#### Step 1: Create a combined (stacked) data set with the randomized trial sample and the
####         target population dataset, with a set of covariates X observed in both groups.
#### Step 2: Create an indicator variable for being in the target population (P).

merged_data <- rbind(trial_data %>% mutate(P=0), target_data %>% mutate(P=1))

#### Step 3: Estimate a model of membership in the population (P) as a function of the
####         covariates X, e.g., using logistic regression.

ps_model <- glm(P ~ x1+x2+x3+x4+x5, data = merged_data, family = binomial)
trial_data$ps <- predict(ps_model, newdata = trial_data, type = "response")

#### Step 4: Create weights for individuals in the trial sample: w_i = ps_i/(1-ps_i)
trial_data$weight <- trial_data$ps / (1 - trial_data$ps)

#### Step 5: Estimate treatment effects using the individuals in the randomized trial by
####         running a weighted regression of outcome as a function of treatment status and
####         the covariates, with the weights calculated in Step 4.
iptw_model <- coxph(Surv(y, failed)~Tr, data = trial_data, weights = weight)
summary(iptw_model)

