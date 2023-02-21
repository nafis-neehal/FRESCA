setwd("/Users/nafisneehal/ESCA")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

data$Smoker<-ifelse(data$Smoker=='No smoke',0,1)

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1000
total_seeds <- 20
IPF_maxiter <- 1500
randomization_ratio <- 1

PATT_Summary <- data.frame(matrix(ncol = 7, nrow = 0))

colsPATT <- c("Seed", "CC_Size", "TA_CC", "TA_HC", "TA_HC_Prop", "TA_HC_IPF", "TA_HC_Both")
colnames(PATT_Summary) <- colsPATT

get_ATT_with_propensity <- function(population){
  m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                     Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                     GFRESTIMATE, data=population, method = "nearest", distance = "glm",
                   replace = TRUE)
  
  matched_data <- match.data(m.out)
  
  #run linear model
  fit <- lm(SEATSYS ~ RANDASSIGN + Age_Group + Gender + Race_or_Ethnicity , 
            data = matched_data, weights = weights)
  
  #calculate ATT
  comp <- comparisons(fit,
                      variables = "RANDASSIGN",
                      vcov = ~ subclass,
                      newdata = subset(matched_data, RANDASSIGN==1),
                      wts = "weights")
  c <- summary(comp)
  return(c$estimate)
}

get_matched_HC <- function(population){
  m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                     Smoker + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                     GFRESTIMATE, data=population, method = "nearest", distance = "glm",
                   replace = TRUE)
  
  matched_data <- match.data(m.out)
  matched_HC <- matched_data %>% filter(RANDASSIGN == 0)
  
  return (matched_HC)
}

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
  
  #calculate ATT
  # comp <- comparisons(fit,
  #                     variables = "RANDASSIGN",
  #                     newdata = subset(population, RANDASSIGN==1))
  # 
  # c <- summary(comp)
  # return(c$estimate)
}

######### Generate DataFrame ##########
#pb <- txtProgressBar(min=0, max=total_seeds, initial = 0)

for (seed_value in 1:total_seeds){
  
  print(paste(seed_value/total_seeds * 100, "% done", sep=""))
  print(paste("Seed ", seed_value))
  print("----------------------------")
  
  set.seed(seed_value)
  TA <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_size)
  
  set.seed(seed_value)
  CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)
  
  external_population <- setdiff(data, rbind(TA_World, CC_World))
  EC_world <- external_population %>% filter(RANDASSIGN==0)
  
  CC_size_list <- seq(TA_size, 0, -50)
  EC_size_list <- randomization_ratio*TA_size - CC_size_list
  
  W_EC_ipf <- get_ECWorld_Bias_weights(dat = EC_world, maxIter = IPF_maxiter) #maybe bias other vars too?
  
  i <- 1
  for (s in CC_size_list){
    
    #print
    print(paste("CC Size ", s))
    
    ATT_TA_CC_list <- list()
    ATT_TA_HC_list <- list()
    ATT_TA_HC_prop_list <- list()
    ATT_TA_HC_IPF_list <- list()
    ATT_TA_HC_both_list <- list()
    
    flag <- 1
    n_plus_sum <- 1
    success_counter <-1
    total_success <- 10
    inner_seed <- 1
    
    while (success_counter <= total_success){
      
      print(paste("n=", success_counter))
      
      set.seed(inner_seed)
      CC <- sample_n(CC_World %>% filter(RANDASSIGN==0), s)
    
      #<- Collect TA_CC (no prop)
      ATT_TA_CC_list <- c(ATT_TA_CC_list, get_ATT_val(rbind(TA, CC), weights = NULL))
    
      for (n_plus_sum in 0:4){ #for each inner_seed value I'll try 5 times to get a pop without missing subcategory
        set.seed(inner_seed + n_plus_sum)
        EC_ipf_sample <- EC_world[sample(seq_len(nrow(EC_world)), size = randomization_ratio*TA_size - s, 
                                         replace = TRUE, prob = W_EC_ipf),]
        
        HC <- rbind(CC, EC_ipf_sample) #unmatched
        matched_HC_df <- get_matched_HC(rbind(TA, HC))
        matched_HC <- matched_HC_df %>% select(MASKID:SEATSYS)
        
        flag <- check_subgroup_counts(age_sbgrps = 2, gender_sbgrps = 2, race_sbgrps = 5, 
                                      smoke_sbgrps = 2, TA = TA, HC = matched_HC)
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
      
      #<- Collect TA_HC (no prop)
      ATT_TA_HC_list <- c(ATT_TA_HC_list, get_ATT_val(rbind(TA, HC), weights = NULL))
    
      #<- Collect TA_HC_prop
      ATT_TA_HC_prop_list <- c(ATT_TA_HC_prop_list, get_ATT_with_propensity(rbind(TA, HC)))
    
      #Apply IPF on TA and HC + Bootstrap
      W_TA_ipf <- get_IPF_weights(dat = TA, maxIter = IPF_maxiter)
      W_HC_ipf <- get_IPF_weights(dat = HC, maxIter = IPF_maxiter)
      W_matched_HC_ipf <- get_IPF_weights(dat = matched_HC, maxIter = IPF_maxiter)
    
      #<- Collect TA_HC_IPF c[0]
      #<- Collect TA_HC_both c[1]
      set.seed(inner_seed)
      TA_IPF_sample <- TA[sample(seq_len(nrow(TA)), size = nrow(TA), replace = TRUE, prob = W_TA_ipf),]
      
      set.seed(inner_seed)
      HC_IPF_sample <- HC[sample(seq_len(nrow(HC)), size = nrow(TA), replace = TRUE, prob = W_HC_ipf),]
      
      set.seed(inner_seed)
      matched_HC_IPF_sample_df <- matched_HC_df[sample(seq_len(nrow(matched_HC_df)), 
                                                       size = nrow(TA), replace = TRUE, 
                                                       prob = W_matched_HC_ipf),]
      matched_HC_IPF_weights <- matched_HC_IPF_sample_df$weights
      
      matched_HC_IPF_sample <- matched_HC_IPF_sample_df %>% select(MASKID:SEATSYS)
      
      ATT_TA_HC_IPF_list <- c(ATT_TA_HC_IPF_list, get_ATT_val(rbind(TA_IPF_sample, HC_IPF_sample), weights = NULL)) #Unmatched data
      ATT_TA_HC_both_list <- c(ATT_TA_HC_both_list, get_ATT_val(rbind(TA_IPF_sample, matched_HC_IPF_sample), weights = c(rep(1, nrow(TA_IPF_sample)),matched_HC_IPF_weights))) #matched HC
      
      inner_seed <- inner_seed + 1
      
    }
    
    #add dataframe
    row <- as.data.frame(list(seed_value, s, 
                              mean(unlist(ATT_TA_CC_list)), 
                              mean(unlist(ATT_TA_HC_list)), 
                              mean(unlist(ATT_TA_HC_prop_list)), 
                              mean(unlist(ATT_TA_HC_IPF_list)), 
                              mean(unlist(ATT_TA_HC_both_list))))
    colnames(row) <- colsPATT
    PATT_Summary <- rbind(PATT_Summary, row)
    
    i <- i+1
  }
  print("----------------------------")
  #setTxtProgressBar(pb, seed_value)
}

directory <- "./Data/Results/PATT_m6/CC_EC_Vary_Size/matchCCEC/"
filename <- "PATT_Summary_Seed20_withReplace_new3_nosbp.csv"
write.csv(PATT_Summary, paste(directory, filename, sep = ""), row.names = FALSE)

#PATT_Summary <- read.csv(paste(directory, filename, sep = ""))

mean_summary <- PATT_Summary %>% 
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
  geom_line(aes(y=TA_CC)) + geom_point(aes(y=TA_CC)) + 
  geom_errorbar(aes(ymin=TA_CC_Low, ymax=TA_CC_High)) + 
  ylim(min_point, max_point) + 
  geom_hline(yintercept=-15.174, linetype="dashed", color = "red") + 
  theme_classic()+
  labs(x="CC_Size", y="PATT Estimate", caption = "TA + CC")
# p2 <- ggplot(mean_summary, aes(x=CC_Size)) +
#   geom_line(aes(y=TA_HC)) + geom_point(aes(y=TA_HC)) + 
#   geom_errorbar(aes(ymin=TA_HC_Low, ymax=TA_HC_High)) +
#   ylim(min_point, max_point) + 
#   geom_hline(yintercept=-15.174, linetype="dashed", color = "red") + 
#   labs(x="CC_Size", y="PATT Estimate", caption = "TA + HC ")
p2 <- ggplot(mean_summary, aes(x=CC_Size)) +
  geom_line(aes(y=TA_HC_Prop)) + geom_point(aes(y=TA_HC_Prop)) + 
  geom_errorbar(aes(ymin=TA_HC_Prop_Low, ymax=TA_HC_Prop_High)) +
  ylim(min_point, max_point) + 
  geom_hline(yintercept=-15.174, linetype="dashed", color = "red") + 
  theme_classic()+
  labs(x="CC_Size", y="PATT Estimate", caption = "TA + HC + Propensity")
p3 <- ggplot(mean_summary, aes(x=CC_Size)) +
  geom_line(aes(y=TA_HC_IPF)) + geom_point(aes(y=TA_HC_IPF)) + 
  geom_errorbar(aes(ymin=TA_HC_IPF_Low, ymax=TA_HC_IPF_High)) +
  ylim(min_point, max_point) + 
  geom_hline(yintercept=-15.174, linetype="dashed", color = "red") + 
  theme_classic()+
  labs(x="CC_Size", y="PATT Estimate", caption = "TA + HC + IPF")
p4 <- ggplot(mean_summary, aes(x=CC_Size)) +
  geom_line(aes(y=TA_HC_Both)) + geom_point(aes(y=TA_HC_Both)) + 
  geom_errorbar(aes(ymin=TA_HC_Both_Low, ymax=TA_HC_Both_High)) +
  ylim(min_point, max_point) + 
  geom_hline(yintercept=-15.174, linetype="dashed", color = "red") + 
  theme_classic()+
  labs(x="CC_Size", y="PATT Estimate", caption = "TA + HC + Propensity + IPF")

p1 + p2 + p3 + p4 + plot_layout(ncol=2, nrow=2) 




######### Mini Experiment for generating Table1 for TA_World, CC_World, EC_World
tw_copy <- TA_World
tw_copy$group <- "TA_World"
ccw_copy <- CC_World
ccw_copy$group <- "CC_World"
ecw_copy <- EC_world[sample(seq_len(nrow(EC_world)), size = nrow(EC_world), 
                                 replace = TRUE, prob = W_EC_ipf),]
ecw_copy$group <- "EC_World"

population <- do.call("rbind", list(tw_copy, ccw_copy, ecw_copy))
population$Education[population$Education=='<HSG'] <- "Below HSG"
population$CVDHISTORY <- factor(population$CVDHISTORY, levels=c(0,1),
                                labels=c("No", "Yes"))

table1(~ Age_Group + Gender + Race_or_Ethnicity + CVDPOINTS + factor(Smoker) + Education +
        FPG + TC + CVDHISTORY +  SERUMCREAT +
         GFRESTIMATE | group, data=population)




