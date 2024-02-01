setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")

#imports
suppressPackageStartupMessages({
  library(generalize)
})

########## Read Data #################
trial_data <- read.csv("./Data/Synthetic/trial_data_paper.csv")
#external_data <- read.csv("./Data/Synthetic/external_data_ate.csv")
target_data <- read.csv("./Data/Synthetic/target_data_paper.csv")
#hybrid_data <- rbind(trial_data, external_data)

########### Regular Treatment Effect
# trial_fit <- lm(y~assignment+x1+x2+x3+x4+x5, data = trial_data)
# summary(trial_fit)
# confint(trial_fit, "assignment", level = 0.95)
# 
# target_fit <- lm(y~assignment+x1+x2+x3+x4+x5, data = target_data)
# summary(target_fit)
# confint(target_fit, "assignment", level = 0.95)

######### Generalize Package ##########
new_trial <- trial_data
new_target <- target_data

new_trial$trial <- 1
new_target$trial <- 0

new_trial$y <- as.double(new_trial$y)
new_target$y <- NA
new_target$failed <- NA
new_target$Tr <- NA

# new_target$y <- NA
# new_target$assignment <- NA

covariates <- c("x1", "x2", "x3", "x4", "x5")

assess_object <- assess(trial="trial", selection_covariates = covariates,
                        data = rbind(new_trial, new_target), selection_method = "lr",
                        trimpop = TRUE)

summary(assess_object)

all_weights <- assess_object$weights
weights <- all_weights[1:nrow(trial_data)]

iptw_model <- coxph(Surv(y, failed)~Tr, data = trial_data, weights = weights)
summary(iptw_model)
# 
# d <- density(weights)
# plot(d$x, d$y, type = 'l', xlab = "", ylab = "")
# title(main = "Distribution of weights derived from Generalize Package (LR) for Simsurv Data", 
#       xlab = "Weight", ylab = "Density")

########### Kaplan Meier Curves for Trial/Target Difference ##########
# new_trial <- trial_data
# new_trial$trial <- "Trial"
# 
# new_target <- target_data
# new_target$trial <- "Target"
# 
# #fit <- survfit(Surv(y, failed) ~ trial, data=rbind(new_trial, new_target))
# fit <- survfit(Surv(y, failed) ~ Tr, data=target_data)
# ggsurvplot(
#   fit,
#   data = trial_data, #rbind(new_trial, new_target),
#   size = 1,                     # change line size
#   # palette =
#   #   c("#E7B800", "#2E9FDF"),  # custom color palettes
#   conf.int = TRUE,              # Add confidence interval
#   pval = F,                     # Add p-value
#   xlim = c(0, 5),               # Add trim in x-axis
#   break.time.by = 1,            # break X axis in time intervals by 500.
#   risk.table = TRUE,            # Add risk table
#   risk.table.col = "strata",    # Risk table color by groups
#   # legend.labs = c("Target", 
#   #                 "Trial"),     # Change legend labels  
#   legend.labs = c("Controls", 
#                   "Treated"),     # Change legend labels
#   risk.table.height = 0.25,     # Useful to change when you have multiple groups
#   ggtheme = theme_bw(),         # Change ggplot2 theme,
#   title = "Survival of Treated VS Controls in Target Data : Simsurv Package"
# )
# 
# ######################################################################
# 
# # weight_object <- weighting(outcome="y", treatment="Tr", trial="trial",
# #                            selection_covariates = covariates, 
# #                            data = rbind(new_trial, new_target), 
# #                            selection_method = "lr")
# # 
# # summary(weight_object)
# 
generalize_object <- generalize(outcome="y", treatment="assignment", trial="trial",
                                selection_covariates = covariates,
                                data = rbind(new_trial, new_target),
                                method = "weighting",
                                selection_method = "lr",
                                trimpop = TRUE)
summary(generalize_object)

# generalize_bart_object <- generalize_bart(outcome="y", treatment="Tr", trial="trial", 
#                                 selection_covariates = covariates, 
#                                 data = rbind(new_trial, new_target))
# summary(generalize_bart_object)

# generalize_tmle_object <- generalize_tmle(outcome="y", treatment="Tr", trial="trial", 
#                                           selection_covariates = covariates, 
#                                           data = rbind(new_trial, new_target))
# summary(generalize_tmle_object)



