setwd("/home/neehan/data/Nafis/ESCAESCA_Git")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1000
total_seeds <- 20
IPF_maxiter <- 100
step_size <- 100
randomization_ratio <- 1
num_col <- 10
CC_size_list <- seq(TA_size-step_size, 100, -step_size)
EC_size_list <- randomization_ratio*TA_size - CC_size_list

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
                             GFRESTIMATE = factor(GFRESTIMATE))
  
  m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                     Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                     GFRESTIMATE, data=p, distance = pscore, replace = T)
  
  matched_data <- match.data(m.out) 
  
  return(matched_data %>% filter(RANDASSIGN==0)) #contains weights
  
}

######### Generate DataFrame ##########
#pb <- txtProgressBar(min=0, max=total_seeds, initial = 0)

run_one_scenario <- function(seed_value, CC_Size){
  
  LDM_Summary <- data.frame(matrix(ncol = 10, nrow = 0))
  colsLDM <- c("Seed", "CC_Size", "Var", "Level", "LDM_TA", "LDM_CC", "LDM_HC_Prop", "LDM_HC_IPF", 
               "LDM_TA_Prop_IPF", "LDM_HC_Prop_IPF" )
  colnames(LDM_Summary) <- colsLDM
  
  set.seed(seed_value)
  TA <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_size)
  W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = IPF_maxiter)
  LDM_TA <- get_LDM_from_count(TA) %>% select(LDM) %>% rename(LDM_TA = LDM)
  
  set.seed(seed_value)
  CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)
  
  external_population <- setdiff(data, rbind(TA_World, CC_World))
  EC_world <- external_population %>% filter(RANDASSIGN==0)
  
  W_EC_bias_ipf <- get_ECWorld_Bias_weights(dat = EC_world, maxIter = IPF_maxiter)
  Biased_EC_world <- EC_world[sample(seq_len(nrow(EC_world)), 
                                     size = nrow(EC_world), 
                                     replace = TRUE, prob = W_EC_bias_ipf),]
  
  #create the propensity score model - propensity to get in the trial
  M <- get_propensity_model(do.call("rbind", list(TA %>% mutate(inTrial = 1),
                                                  CC_World %>% mutate(inTrial = 1),
                                                  Biased_EC_world %>% mutate(inTrial = 0))))
  
  LDM_CC_list <- list()
  LDM_HC_Prop_list <- list()
  LDM_HC_IPF_list <- list()
  LDM_TA_Prop_IPF_list <- list()
  LDM_HC_Prop_IPF_list <- list()
  
  flag <- 1
  n_plus_sum <- 1
  success_counter <-1
  total_success <- 10
  inner_seed <- 1
  
  while (success_counter <= total_success){ #loops 10 times
    
    print(paste("n=", success_counter))
    
    set.seed(inner_seed)
    CC <- sample_n(CC_World %>% filter(RANDASSIGN==0), CC_Size)
    
    #<- Collect LDM_CC 
    LDM_CC_three <- get_LDM_from_count(CC) %>% select(VAR, Level, LDM)
    VAR <- LDM_CC_three$VAR
    Level <- LDM_CC_three$Level 
    
    LDM_CC <- LDM_CC_three %>% select(LDM) %>% rename(LDM_CC = LDM)
    LDM_CC_list <- c(LDM_CC_list, LDM_CC)
    
    for (n_plus_sum in 0:4){ #for each inner_seed value I'll try 5 times to get a pop without missing subcategory
      
      set.seed(inner_seed + n_plus_sum)
      unmatched_EC <- sample_n(Biased_EC_world, randomization_ratio * TA_size - CC_Size)
      
      #get propensity score of this unmatched EC sample
      ps <- get_propensity_score(M, rbind(TA, unmatched_EC)) #contains MASKID, RANDASSIGN and PRS
      
      #get matched EC with weights
      matched_EC_df <- get_matched_data(rbind(TA, unmatched_EC), unlist(ps$pr_score))
      matched_EC <- matched_EC_df %>% select(MASKID:SEATSYS)
      matched_EC_weights <- matched_EC_df$weights #contains only EC weights
      
      
      #create matched and unmatched HC
      unmatched_HC <- rbind(CC, unmatched_EC)
      matched_HC <- rbind(CC, matched_EC) #just contains maskid:seatsys
      matched_HC_weights <- c(rep(1, nrow(CC)), matched_EC_weights)
      matched_HC_df <- matched_HC
      matched_HC_df$weights <- matched_HC_weights
      
      flag <- check_subgroup_counts(age_sbgrps = 2, gender_sbgrps = 2, race_sbgrps = 5, 
                                    gfr_sbgrps = 2, TA = TA, HC = matched_HC)
      if (flag==1){
        break
      }
      else{
        flag<- 0
        print("Changing seed inside 1..n loop due to missing subcategory!")
      }
    }
    
    if(flag==0){
      inner_seed <- inner_seed + 5
      next
    }
    
    success_counter <- success_counter + 1
    
    #<- Collect LDM_HC_Prop
    LDM_HC_Prop <- get_LDM_from_count(matched_HC) %>% select(LDM) %>% rename(LDM_HC_Prop = LDM)
    LDM_HC_Prop_list <- c(LDM_HC_Prop_list, LDM_HC_Prop)
    
    #Apply IPF on TA and matched HC
    W_unmatched_HC_ipf <- get_IPF_weights(dat = unmatched_HC, maxIter = IPF_maxiter) #no propensity
    W_matched_HC_ipf <- get_IPF_weights(dat = matched_HC, maxIter = IPF_maxiter) #propensity
    
    #<- Collect TA_HC_IPF c[0]
    #<- Collect TA_HC_both c[1]
    set.seed(inner_seed)
    TA_IPF_sample <- TA[sample(seq_len(nrow(TA)), size = nrow(TA), replace = TRUE, prob = W_TA_ipf),]
    
    set.seed(inner_seed)
    unmatched_HC_IPF_sample <- unmatched_HC[sample(seq_len(nrow(unmatched_HC)), 
                                                   size = nrow(TA), replace = TRUE, prob = W_unmatched_HC_ipf),]
    
    set.seed(inner_seed)
    matched_HC_IPF_sample_df <- matched_HC_df[sample(seq_len(nrow(matched_HC_df)), 
                                                     size = nrow(TA), replace = TRUE, prob = W_matched_HC_ipf),]
    
    matched_HC_IPF_sample <- matched_HC_IPF_sample_df %>% select(MASKID:SEATSYS)
    matched_HC_IPF_sample_weights <- matched_HC_IPF_sample_df$weights
    
    #<- Collect LDM_TA_IPF, LDM_HC_IPF
    LDM_TA_Prop_IPF <- get_LDM_from_count(TA_IPF_sample) %>% select(LDM) %>% rename(LDM_TA_Prop_IPF = LDM)
    LDM_TA_Prop_IPF_list <- c(LDM_TA_Prop_IPF_list, LDM_TA_Prop_IPF)
    
    LDM_HC_IPF <- get_LDM_from_count(unmatched_HC_IPF_sample) %>% select(LDM) %>% rename(LDM_HC_IPF = LDM)
    LDM_HC_IPF_list <- c(LDM_HC_IPF_list, LDM_HC_IPF)
    
    LDM_HC_Prop_IPF <- get_LDM_from_count(matched_HC_IPF_sample) %>% select(LDM) %>% rename(LDM_HC_Prop_IPF = LDM)
    LDM_HC_Prop_IPF_list <- c(LDM_HC_Prop_IPF_list, LDM_HC_Prop_IPF)
    
    inner_seed <- inner_seed + 1
    
  }
  
  #add dataframe
  row <- as.data.frame(list(seed_value, CC_Size, VAR, Level,
                            abs(LDM_TA),
                            abs(rowMeans(as.data.frame(LDM_CC_list))),
                            abs(rowMeans(as.data.frame(LDM_HC_Prop_list))),
                            abs(rowMeans(as.data.frame(LDM_HC_IPF_list))),
                            abs(rowMeans(as.data.frame(LDM_TA_Prop_IPF_list))),
                            abs(rowMeans(as.data.frame(LDM_HC_Prop_IPF_list)))))
  colnames(row) <- colsLDM
  LDM_Summary <- rbind(LDM_Summary, row)
  return(LDM_Summary)
}

######### Generate DataFrame ##########
Final_LDM_Summary <- data.frame(matrix(ncol = 10, nrow = 0))
colsLDM <- c("Seed", "CC_Size", "Var", "Level", "LDM_TA", "LDM_CC", "LDM_HC_Prop", "LDM_HC_IPF", 
             "LDM_TA_Prop_IPF", "LDM_HC_Prop_IPF" )
colnames(Final_LDM_Summary) <- colsLDM

cores=detectCores(logical = TRUE)
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoSNOW(cl)

packages <- c("dplyr", "ipfr", "MatchIt", "marginaleffects")

pb <- txtProgressBar(max = length(total_seeds * length(CC_size_list)), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

Final_LDM_Summary <- foreach(seed_value=seq(1,total_seeds,1), .combine = 'rbind', 
                             .packages = packages, 
                             .options.snow = opts)%:%
  foreach (size=CC_size_list, .combine = "rbind", .packages = packages) %dopar%{
    run_one_scenario(seed_value, size)
    }


close(pb)
stopCluster(cl)
stopImplicitCluster()

directory <- "./Data/Results/PATT_m6/CC_EC_Vary_Size/matchCCEC/"
filename <- "LDM_Summary_Seed20_nested_optimized_gfr_proptrial_v1.csv"
write.csv(Final_LDM_Summary, paste(directory, filename, sep = ""), row.names = FALSE)

#Final_LDM_Summary <- read.csv(paste(directory, filename, sep = ""))

#Final_LDM_Summary <- Final_LDM_Summary %>% filter(!CC_Size==0)

Inf_replace <- max(Final_LDM_Summary %>% mutate_all(function(x) ifelse(is.infinite(x), -1, x)) %>% 
                     select(LDM_CC:LDM_HC_Prop_IPF)) + 0.5
LDM_Summary_Inf <- Final_LDM_Summary
LDM_Summary_Inf[sapply(LDM_Summary_Inf, is.infinite)] <- Inf_replace

mean_summary <- LDM_Summary_Inf %>%
  group_by(CC_Size, Level) %>%
  summarise(median(LDM_TA),
            LDM_TA_Low = confidence_interval(unlist(LDM_TA), 0.95)[1],
            LDM_TA_High = confidence_interval(unlist(LDM_TA), 0.95)[2],
            median(LDM_CC),
            LDM_CC_Low = confidence_interval(unlist(LDM_CC), 0.95)[1],
            LDM_CC_High = confidence_interval(unlist(LDM_CC), 0.95)[2],
            median(LDM_HC_Prop),
            LDM_HC_Prop_Low = confidence_interval(unlist(LDM_HC_Prop), 0.95)[1],
            LDM_HC_Prop_High = confidence_interval(unlist(LDM_HC_Prop), 0.95)[2],
            median(LDM_HC_IPF),
            LDM_HC_IPF_Low = confidence_interval(unlist(LDM_HC_IPF), 0.95)[1],
            LDM_HC_IPF_High = confidence_interval(unlist(LDM_HC_IPF), 0.95)[2],
            median(LDM_TA_Prop_IPF),
            LDM_TA_Prop_IPF_Low = confidence_interval(unlist(LDM_TA_Prop_IPF), 0.95)[1],
            LDM_TA_Prop_IPF_High = confidence_interval(unlist(LDM_TA_Prop_IPF), 0.95)[2],
            median(LDM_HC_Prop_IPF),
            LDM_HC_Prop_IPF_Low = confidence_interval(unlist(LDM_HC_Prop_IPF), 0.95)[1],
            LDM_HC_Prop_IPF_High = confidence_interval(unlist(LDM_HC_Prop_IPF), 0.95)[2],
  ) %>%
  rename(LDM_TA = "median(LDM_TA)",
         LDM_CC = "median(LDM_CC)",
         LDM_HC_Prop = "median(LDM_HC_Prop)",
         LDM_HC_IPF = "median(LDM_HC_IPF)",
         LDM_TA_Prop_IPF = "median(LDM_TA_Prop_IPF)",
         LDM_HC_Prop_IPF = "median(LDM_HC_Prop_IPF)") %>%
  mutate(LDM_CC_Low = ifelse(LDM_CC==Inf_replace, LDM_CC, LDM_CC_Low),
         LDM_CC_High = ifelse(LDM_CC==Inf_replace, LDM_CC, LDM_CC_High)) %>%
  mutate(Plot = ifelse(LDM_CC==Inf_replace, paste(as.character(CC_Size),"*"),
                       as.character(CC_Size))) %>%
  ungroup()


# mean_summary <- LDM_Summary_Inf %>% 
#   group_by(CC_Size, Level) %>% 
#   summarise(median(LDM_TA), 
#             LDM_TA_Low = MedianCI(unlist(LDM_TA), 0.95)[2], 
#             LDM_TA_High = MedianCI(unlist(LDM_TA), 0.95)[3],
#             median(LDM_CC), 
#             LDM_CC_Low = MedianCI(unlist(LDM_CC), 0.95)[2], 
#             LDM_CC_High = MedianCI(unlist(LDM_CC), 0.95)[3],
#             median(LDM_HC_Prop), 
#             LDM_HC_Prop_Low = MedianCI(unlist(LDM_HC_Prop), 0.95)[2], 
#             LDM_HC_Prop_High = MedianCI(unlist(LDM_HC_Prop), 0.95)[3],
#             median(LDM_HC_IPF), 
#             LDM_HC_IPF_Low = MedianCI(unlist(LDM_HC_IPF), 0.95)[2], 
#             LDM_HC_IPF_High = MedianCI(unlist(LDM_HC_IPF), 0.95)[3],
#             median(LDM_TA_Prop_IPF), 
#             LDM_TA_Prop_IPF_Low = MedianCI(unlist(LDM_TA_Prop_IPF), 0.95)[2], 
#             LDM_TA_Prop_IPF_High = MedianCI(unlist(LDM_TA_Prop_IPF), 0.95)[3],
#             median(LDM_HC_Prop_IPF), 
#             LDM_HC_Prop_IPF_Low = MedianCI(unlist(LDM_HC_Prop_IPF), 0.95)[2], 
#             LDM_HC_Prop_IPF_High = MedianCI(unlist(LDM_HC_Prop_IPF), 0.95)[3],
#   ) %>%
#   rename(LDM_TA = "median(LDM_TA)",
#          LDM_CC = "median(LDM_CC)",
#          LDM_HC_Prop = "median(LDM_HC_Prop)",
#          LDM_HC_IPF = "median(LDM_HC_IPF)",
#          LDM_TA_Prop_IPF = "median(LDM_TA_Prop_IPF)",
#          LDM_HC_Prop_IPF = "median(LDM_HC_Prop_IPF)") %>%
#   mutate(LDM_CC_Low = ifelse(LDM_CC==Inf_replace, LDM_CC, LDM_CC_Low),
#          LDM_CC_High = ifelse(LDM_CC==Inf_replace, LDM_CC, LDM_CC_High)) %>%
#   mutate(Plot = ifelse(LDM_CC==Inf_replace, paste(as.character(CC_Size),"*"), 
#                        as.character(CC_Size))) %>%
#   ungroup()



return_plot <- function(var, level){
  ymin_val <- min(ungroup(mean_summary) %>% filter(Level==level) %>% 
                    select(LDM_CC:LDM_HC_Prop_IPF_High))
  ymax_val <- max(ungroup(mean_summary) %>% filter(Level==level) %>% 
                    select(LDM_CC:LDM_HC_Prop_IPF_High))
  cap <- paste(var,":",level)
  p <- mean_summary %>% filter(Level==level) %>%
    #ggplot(aes(x=factor(Plot, levels = mixedsort(Plot)))) +
    ggplot(aes(x=CC_Size)) +
    # geom_line(aes(y=LDM_CC, color="Only CC")) + 
    # geom_errorbar(aes(ymin=LDM_CC_Low, ymax=LDM_CC_High, colour="Only CC")) +
    # geom_point(aes(y=LDM_CC, colour="Only CC")) +
    # geom_line(aes(y=LDM_HC_Prop, color="HC+Propensity")) +
    # geom_errorbar(aes(ymin=LDM_HC_Prop_Low, ymax=LDM_HC_Prop_High, colour="HC+Propensity")) +
    # geom_point(aes(y=LDM_HC_Prop, color="HC+Propensity")) +
    # geom_line(aes(y=LDM_HC_IPF, color="HC+IPF")) +
    # geom_errorbar(aes(ymin=LDM_HC_IPF_Low, ymax=LDM_HC_IPF_High, colour="HC+IPF")) +
    # geom_point(aes(y=LDM_HC_IPF, color="HC+IPF")) +
    geom_line(aes(y=LDM_HC_Prop_IPF, color="HC+Propensity+IPF")) +
    geom_errorbar(aes(ymin=LDM_HC_Prop_IPF_Low, ymax=LDM_HC_Prop_IPF_High, color="HC+Propensity+IPF")) +
    geom_point(aes(y=LDM_HC_Prop_IPF, color="HC+Propensity+IPF")) +
    geom_hline(aes(yintercept= 0.22, color = "Equity Range"), linetype="dashed")+
    #ylim(ymin_val,ymax_val)+
    theme_classic()+
    theme(axis.title =element_text(size=6), axis.text.x = element_text(angle = 90))+
    scale_color_manual(name="Population",
                       values=c("Only CC"="#D55E00", "HC+Propensity"="#336699", "HC+IPF"="#CC66FF",
                                "HC+Propensity+IPF"="#FFCC00", "Equity Range"="#99CC33"))+
    labs(x="CC_Size", y="Median LDM Estimate", caption = cap)
  return (p)
}

g1 <- return_plot("Gender", "Female")
g1
g2 <- return_plot("Gender", "Male")
a1 <- return_plot("Age_Group", "40-59")
a2 <- return_plot("Age_Group", "59+")
r1 <- return_plot("Race_or_Ethnicity", "NH Asian")
r2 <- return_plot("Race_or_Ethnicity", "NH Black")
r3 <- return_plot("Race_or_Ethnicity", "NH White")
r4 <- return_plot("Race_or_Ethnicity", "Hispanic")
r5 <- return_plot("Race_or_Ethnicity", "Other")

ggarrange(g1,g2, nrow=1,
          common.legend = TRUE, legend = "bottom")

ggarrange(a1,a2, nrow=1,
          common.legend = TRUE, legend = "bottom")

ggarrange(r1, r2, r3, r4, r5, nrow=2, ncol=3,
          common.legend = TRUE, legend = "bottom")

