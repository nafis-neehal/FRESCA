setwd("/Users/nafisneehal/ESCA")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1000
total_seeds <- 5
IPF_maxiter <- 2000
randomization_ratio <- 1

HC_LDM_FINAL <- data.frame(matrix(ncol = 5, nrow = 0))
TA_LDM_FINAL <- data.frame(matrix(ncol = 4, nrow = 0))
colsHC <- c("Seed", "CC_Size", "VAR", "Level", "LDM_HC")
colsTA <- c("Seed", "VAR", "Level", "LDM_TA")
colnames(HC_LDM_FINAL) <- colsHC
colnames(TA_LDM_FINAL) <- colsTA

ATT_FINAL <- data.frame(matrix(ncol = 4, nrow = 0))
colsATT <- c("Seed", "CC_Size", "PATT_Estimate")
colnames(ATT_FINAL) <- colsATT

######### Generate DataFrame ##########
pb <- txtProgressBar(min=1, max=total_seeds, initial = 1)
for (seed_value in 1:total_seeds){
  
  #measure TA LDI Value
  set.seed(seed_value)
  TA <- sample_n(TA_World %>% filter(RANDASSIGN==1), TA_size)
  LDM_TA <- get_LDM_from_count(TA) %>% select(VAR, Level, LDM) %>% rename(LDM_TA = LDM)
  LDM_TA_summary <- as.data.frame(LDM_TA)
  LDM_TA_summary$Seed <- seed_value
  TA_LDM_FINAL <- rbind(TA_LDM_FINAL, LDM_TA_summary %>% select(Seed, VAR, Level, LDM_TA))
  
  #measure CC, EC and HC LDI Value
  set.seed(seed_value)
  CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)
  
  external_population <- setdiff(data, rbind(TA_World, CC_World))
  EC_world <- external_population %>% filter(RANDASSIGN==0)
  
  CC_size_list <- seq(TA_size, 50, -50)
  EC_size_list <- randomization_ratio*TA_size - CC_size_list
  
  W_EC_ipf <- get_ECWorld_Bias_weights(dat = EC_world, maxIter = IPF_maxiter)
  
  i <- 1
  for (s in CC_size_list){
    set.seed(i)
    CC <- sample_n(CC_World %>% filter(RANDASSIGN==0), s)
    
    set.seed(i)
    EC_ipf_sample <- EC_world[sample(seq_len(nrow(EC_world)), size = randomization_ratio*TA_size - s, 
                                     replace = TRUE, prob = W_EC_ipf),]
    
    HC <- rbind(CC, EC_ipf_sample)
    
    #calculate ATT
    merged_sample <- rbind(TA, HC)
    
    m.out <- matchit(RANDASSIGN ~ Age_Group + Gender + Race_or_Ethnicity + Education +
                       Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                       GFRESTIMATE, data=merged_sample, method = "nearest", distance = "glm")
    matched_data <- match.data(m.out)
    
    #run linear model to get the treatment effect
    fit <- lm(SEATSYS ~ RANDASSIGN + Age_Group + Gender + Race_or_Ethnicity + Education + 
                Smoker + SBP + FPG + TC + CVDHISTORY + CVDPOINTS + SERUMCREAT +
                GFRESTIMATE, data = matched_data, weights = weights)
    
    #calculate ATT
    ATT <- coeftest(fit, vcov. = vcovCL, cluster = ~subclass)["RANDASSIGN",,drop = FALSE][1]
    
    #Append ATT
    row <- data.frame(list(seed_value, s, ATT))
    colnames(row) <- colsATT
    ATT_FINAL <- rbind(ATT_FINAL, row)
    
    #calculate LDM
    LDM_HC <- get_LDM_from_count(HC) %>% select(VAR, Level, LDM) %>% rename(LDM_HC = LDM)
    LDM_HC_summary <- as.data.frame(LDM_HC)
    
    LDM_HC_summary$CC_Size <- s
    LDM_HC_summary$Seed <- seed_value
    
    HC_LDM_FINAL <- rbind(HC_LDM_FINAL, LDM_HC_summary%>%select(Seed, CC_Size, VAR, Level, LDM_HC))
    
    i <- i+1
  }
  setTxtProgressBar(pb, seed_value)
}
close(pb)



########### Preprocess for Plot ##########
HC_LDM_FINAL$LDM_HC <- abs(HC_LDM_FINAL$LDM_HC)
HC_LDM_FINAL <- HC_LDM_FINAL %>% mutate(LDM_Bins = cut(LDM_HC, breaks = c(0, 0.21, 0.50, 10000, Inf)))

dat_num_max <- HC_LDM_FINAL %>% filter(LDM_HC != Inf) %>% select(LDM_HC) %>% max()
HC_LDM_FINAL <- HC_LDM_FINAL %>% mutate(LDM_HC_NEW = ifelse(LDM_HC==Inf, dat_num_max+1, LDM_HC))
HC_LDM_FINAL$Seed <- as.factor(HC_LDM_FINAL$Seed)

########### Plotting of Results ###########
## HC Plot ##
group.colors <- c(Missing="red",Present= "black")
axis_label_size <- 7
ylabel <- "Absolute of\n Log Disparate Value"
xlabel <- "CC Sample Size"

# text <- "Examining Absolute Disparate Value with varying CC Size and Seed"
# tgrob <- text_grob(text, size=15)
# plot_0 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,0,0,0,"cm"))

#Gender
var.name <- "Missing_Subgroups"
G1 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='Female') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw()+
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Gender: Female") +
  scale_color_manual(values=group.colors)

G2 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='Male') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Gender: Male") +
  scale_color_manual(values=group.colors)

#Ages
A1 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='40-59') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Age Group: 40-59") +
  scale_color_manual(values=group.colors)

A2 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='59+') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Age Group: 59+") +
  scale_color_manual(values=group.colors)

#Races
R1 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='NH White') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: NH White") +
  scale_color_manual(values=group.colors)

R2 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='NH Black') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: NH Black") +
  scale_color_manual(values=group.colors)

R3 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='NH Asian') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: NH Asian") +
  scale_color_manual(values=group.colors)

R4 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='Hispanic') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: Hispanic") +
  scale_color_manual(values=group.colors)

R5 <- HC_LDM_FINAL %>% mutate(Color = ifelse(LDM_HC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='Other') %>%
  ggplot(aes(x=CC_Size, y=LDM_HC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  theme(axis.title=element_text(size = axis_label_size)) +
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: Other") +
  scale_color_manual(values=group.colors)

ggarrange(G1, G2, A1, A2, R1, R2, R3, R4, R5, nrow=3, ncol=3,
          common.legend = TRUE, legend = "bottom")

### TA Plot ###
TA_LDM_FINAL$LDM_TA <- abs(TA_LDM_FINAL$LDM_TA)
TA_LDM_FINAL <- TA_LDM_FINAL %>% mutate(LDM_Bins = cut(LDM_TA, breaks = c(0, 0.21, 0.50, 10000, Inf)))

dat_num_max <- TA_LDM_FINAL %>% filter(LDM_TA != Inf) %>% select(LDM_TA) %>% max()
TA_LDM_FINAL <- TA_LDM_FINAL %>% mutate(LDM_TA_NEW = ifelse(LDM_TA==Inf, dat_num_max+1, LDM_TA))
TA_LDM_FINAL$Seed <- as.factor(TA_LDM_FINAL$Seed)

TA_LDM_FINAL %>% mutate(Missing_Subgroups = ifelse(LDM_TA==Inf, "Missing", "Present")) %>%
  ggplot(aes(x=Level, y=LDM_TA_NEW, color=Missing_Subgroups, shape=Seed)) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() +
  scale_color_manual(values=group.colors)

### ATT Plot ###
ATT_FINAL$Seed <- as.factor(ATT_FINAL$Seed)
ggplot(ATT_FINAL, aes(x=CC_Size, y=PATT_Estimate, color=Seed)) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept=-15.18, linetype="dashed", color = "black") +
  scale_color_brewer(palette = "Set1")


# ## Test Code
# cc_copy <- rbind(TA_World, CC_World)
# ec_copy <- EC_world[sample(seq_len(nrow(EC_world)), size = nrow(EC_world), 
#                            replace = TRUE, prob = W_EC_ipf),]
# cc_copy$group <- "In Trial"
# ec_copy$group <- "External"
# m <- rbind(cc_copy, ec_copy)
# 
# table1(~ Age_Group + Gender + Race_or_Ethnicity | group, data=m)

### Plot Varying CC and show ratio get's biased
ta_copy <- TA
Count_Final <- data.frame(matrix(ncol = 7, nrow = 0))
cols <- c("CC_Size", "VAR", "Level", "HC_Count", "HC_Ratio", "TA_Count", "TA_Ratio")
vars <- c("Age_Group", "Age_Group", "Gender", "Gender", "Race_or_Ethnicity",
          "Race_or_Ethnicity","Race_or_Ethnicity","Race_or_Ethnicity","Race_or_Ethnicity")
levels <- c("40-59", "59+", "Female", "Male", "Hispanic", "NH Asian", 
            "NH Black", "NH White", "Other")
colnames(Count_Final) <- cols

for (s in CC_size_list){
  cc1 <- sample_n(CC_World %>% filter(RANDASSIGN==0), s)
  ec1 <- EC_world[sample(seq_len(nrow(EC_world)), size = randomization_ratio*TA_size - s, 
                                   replace = TRUE, prob = W_EC_ipf),]
  hc_copy <- rbind(cc1, ec1)
  for (i in 1:length(levels)){
    count_in_ta <- ta_copy %>% filter(!!sym(vars[i])==levels[i]) %>% count()
    count_in_hc <- hc_copy %>% filter(!!sym(vars[i])==levels[i]) %>% count()
    ratio_in_ta <- count_in_ta / nrow(ta_copy)
    ratio_in_hc <- count_in_hc / nrow(hc_copy)
    row <- data.frame(s, vars[i], levels[i], count_in_hc, ratio_in_hc, 
                      count_in_ta, ratio_in_ta)
    colnames(row) <- cols
    Count_Final <- rbind(Count_Final, row)
  }
}

ggplot(Count_Final, aes(x=CC_Size)) + 
  geom_line(aes(y=HC_Ratio, color=Level), linetype="dashed") + 
  geom_line(aes(y=TA_Ratio, color=Level))







