setwd("/home/neehan/data/Nafis/ESCA_Git")
source("./Scripts/common.R")
source("./Scripts/run_patt_scenario.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1000
total_seeds <- 20
IPF_maxiter <- 100
step_size <- 100
randomization_ratio <- 1
num_col <- 7
CC_size_list <- seq(TA_size-step_size, 100, -step_size)
EC_size_list <- randomization_ratio*TA_size - CC_size_list


confidence_interval <- function(vector, interval){
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


get_ATT_val <- function(population, weights){
  
  #run linear model to get the treatment effect
  if (!is.null(weights)){
    fit <- lm(SEATSYS ~ RANDASSIGN, data = population, weights = weights)
  }
  else{
    fit <- lm(SEATSYS ~ RANDASSIGN, data = population)
  }
  
  estimate <- coef(fit)['RANDASSIGN']
  
  return(estimate)
}

get_propensity_model <- function(population){
  m_ps <- glm(inTrial ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                GFRESTIMATE, data = population, family = binomial())
  return (m_ps)
}

get_propensity_score <- function(Model, population){
  prs_df <- data.frame(MASKID = population %>% select(MASKID), 
                       RANDASSIGN = population %>% select(RANDASSIGN), 
                       pr_score = predict(Model, population, type = "response"))
  return (prs_df)
}

get_matched_data <- function(population, pscore){
  p <- population %>% mutate(Age_Group = factor(Age_Group),
               Gender = factor(Gender),
               Race_or_Ethnicity = factor(Race_or_Ethnicity),
               Education = factor(Education),
               Smoker = factor(Smoker),
               FPG = factor(FPG),
               TC = factor(TC),
               CVDPOINTS = factor(CVDPOINTS),
               GFRESTIMATE = factor(GFRESTIMATE))
  
  m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                     Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                     GFRESTIMATE, data=p, distance = pscore, replace = T)
  
  matched_data <- match.data(m.out) 
  
  return(matched_data %>% filter(RANDASSIGN==0)) #contains weights
  
}



######### Generate DataFrame ##########

Final_PATT_Summary <- data.frame(matrix(ncol = num_col, nrow = 0))
colsPATT <- c("Seed", "CC_Size", "TA_CC", "TA_HC", "TA_HC_Prop", "TA_HC_IPF", "TA_HC_Both")
colnames(Final_PATT_Summary) <- colsPATT

cores=detectCores(logical = TRUE)
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoSNOW(cl)

packages <- c("dplyr", "ipfr", "MatchIt", "marginaleffects")

pb <- txtProgressBar(max = length(total_seeds * length(CC_size_list)), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

Final_PATT_Summary <- foreach(seed_value=seq(1,total_seeds,1), .combine = 'rbind', 
                              .packages = packages, 
                              .options.snow = opts) %:%
      foreach (size=CC_size_list, .combine = "rbind", .packages = packages) %dopar%{
                                   run_one_patt_scenario(seed_value, size)
  }


close(pb)
stopCluster(cl)
stopImplicitCluster()

directory <- "./Data/Results/PATT_m6/"
filename <- "PATT_Summary_Seed20.csv"
write.csv(Final_PATT_Summary, paste(directory, filename, sep = ""), row.names = FALSE)
Final_PATT_Summary <- read.csv(paste(directory, filename, sep = ""))

mean_summary <- Final_PATT_Summary %>% 
  group_by(CC_Size) %>% 
  summarise(mean(TA_CC), 
            TA_CC_Low = confidence_interval(TA_CC, 0.95)[1], 
            TA_CC_High = confidence_interval(TA_CC, 0.95)[2],
            mean(TA_HC), 
            TA_HC_Low = confidence_interval(TA_HC, 0.95)[1], 
            TA_HC_High = confidence_interval(TA_HC, 0.95)[2],
            mean(TA_HC_Prop), 
            TA_HC_Prop_Low = confidence_interval(TA_HC_Prop, 0.95)[1], 
            TA_HC_Prop_High = confidence_interval(TA_HC_Prop, 0.95)[2],
            mean(TA_HC_IPF), 
            TA_HC_IPF_Low = confidence_interval(TA_HC_IPF, 0.95)[1], 
            TA_HC_IPF_High = confidence_interval(TA_HC_IPF, 0.95)[2], 
            mean(TA_HC_Both), 
            TA_HC_Both_Low = confidence_interval(TA_HC_Both, 0.95)[1], 
            TA_HC_Both_High = confidence_interval(TA_HC_Both, 0.95)[2]) %>%
  rename(TA_CC = "mean(TA_CC)",
         TA_HC = "mean(TA_HC)",
         TA_HC_Prop = "mean(TA_HC_Prop)",
         TA_HC_IPF = "mean(TA_HC_IPF)",
         TA_HC_Both = "mean(TA_HC_Both)") %>%
  filter(TA_CC != 0)



#mean_summary_long <- gather(mean_summary, Population, Estimate, TA_CC:TA_HC_Both_High)

min_point <- min(mean_summary %>% select(TA_CC:TA_HC_Both_High)) - 0.25
max_point <- max(mean_summary %>% select(TA_CC:TA_HC_Both_High)) + 0.25

p1 <- ggplot(mean_summary, aes(x=CC_Size)) + 
  geom_line(aes(y=TA_CC, colour="Only CC")) + geom_point(aes(y=TA_CC, colour="Only CC")) + 
  geom_errorbar(aes(ymin=TA_CC_Low, ymax=TA_CC_High, colour="Only CC")) + 
  ylim(min_point, max_point) + 
  geom_hline(yintercept=-15.2, linetype="dashed", color = "red") + 
  theme_classic()+
  labs(x="CC_Size", y="PATT Estimate", caption = "TA + CC") + 
  scale_color_manual(values=c("Only CC"="#D55E00")) 
# p5 <- ggplot(mean_summary, aes(x=CC_Size)) +
#   geom_line(aes(y=TA_HC)) + geom_point(aes(y=TA_HC)) +
#   geom_errorbar(aes(ymin=TA_HC_Low, ymax=TA_HC_High)) +
#   ylim(min_point, max_point) +
#   geom_hline(yintercept=-15.2, linetype="dashed", color = "red") +
#   theme_classic()+
#   labs(x="CC_Size", y="PATT Estimate", caption = "TA + HC ")
p2 <- ggplot(mean_summary, aes(x=CC_Size)) +
  geom_line(aes(y=TA_HC_Prop, color="HC+Propensity")) + geom_point(aes(y=TA_HC_Prop, color="HC+Propensity")) + 
  geom_errorbar(aes(ymin=TA_HC_Prop_Low, ymax=TA_HC_Prop_High, color="HC+Propensity")) +
  ylim(min_point, max_point) + 
  geom_hline(yintercept=-15.2, linetype="dashed", color = "red") + 
  theme_classic()+
  labs(x="CC_Size", y="PATT Estimate", caption = "TA + HC + Propensity") + 
  scale_color_manual(values=c("HC+Propensity"="#336699")) 
p3 <- ggplot(mean_summary, aes(x=CC_Size)) +
  geom_line(aes(y=TA_HC_IPF, color="HC+IPF")) + geom_point(aes(y=TA_HC_IPF, color="HC+IPF")) + 
  geom_errorbar(aes(ymin=TA_HC_IPF_Low, ymax=TA_HC_IPF_High, color="HC+IPF")) +
  ylim(min_point, max_point) + 
  geom_hline(yintercept=-15.2, linetype="dashed", color = "red") + 
  theme_classic()+
  labs(x="CC_Size", y="PATT Estimate", caption = "TA + HC + IPF") + 
  scale_color_manual(values=c("HC+IPF"="#CC66FF")) 
p4 <- ggplot(mean_summary, aes(x=CC_Size)) +
  geom_line(aes(y=TA_HC_Both, color="HC+Propensity+IPF")) + geom_point(aes(y=TA_HC_Both, color="HC+Propensity+IPF")) + 
  geom_errorbar(aes(ymin=TA_HC_Both_Low, ymax=TA_HC_Both_High, color="HC+Propensity+IPF")) +
  ylim(min_point, max_point) + 
  geom_hline(yintercept=-15.2, linetype="dashed", color = "red") + 
  theme_classic()+
  labs(x="CC_Size", y="PATT Estimate", caption = "TA + HC + Propensity + IPF") + 
  scale_color_manual(values=c("HC+Propensity+IPF"="#FFCC00")) 

p1 + p2 + p3 + p4 + 
  plot_layout(ncol=2, nrow=2) 
