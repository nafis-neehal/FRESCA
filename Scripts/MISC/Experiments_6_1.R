setwd("/Users/nafisneehal/ESCA")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

TA_World <- data %>% filter(RANDASSIGN==1)
TA_size <- 1500
total_seeds <- 5

CC_LDM_FINAL <- data.frame(matrix(ncol = 5, nrow = 0))
TA_LDM_FINAL <- data.frame(matrix(ncol = 4, nrow = 0))
colsCC <- c("Seed", "CC_Size", "VAR", "Level", "LDM_CC")
colsTA <- c("Seed", "VAR", "Level", "LDM_TA")
colnames(CC_LDM_FINAL) <- colsCC
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
  
  #measure CC LDI Value
  set.seed(seed_value)
  CC_World <- sample_n(data %>% filter(RANDASSIGN==0), TA_size)
  
  CC_size_list <- seq(TA_size, 50, -50)
  
  i <- 1
  for (s in CC_size_list){
    set.seed(i)
    CC <- sample_n(CC_World %>% filter(RANDASSIGN==0), s)
    
    #calculate ATT
    merged_sample <- rbind(TA, CC)
    
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
    LDM_CC <- get_LDM_from_count(CC) %>% select(VAR, Level, LDM) %>% rename(LDM_CC = LDM)
    LDM_CC_summary <- as.data.frame(LDM_CC)
    
    LDM_CC_summary$CC_Size <- s
    LDM_CC_summary$Seed <- seed_value
    
    CC_LDM_FINAL <- rbind(CC_LDM_FINAL, LDM_CC_summary%>%select(Seed, CC_Size, VAR, Level, LDM_CC))
        
    i <- i+1
  }
  setTxtProgressBar(pb, seed_value)
}
close(pb)

#preprocess for plot
CC_LDM_FINAL$LDM_CC <- abs(CC_LDM_FINAL$LDM_CC)
CC_LDM_FINAL <- CC_LDM_FINAL %>% mutate(LDM_Bins = cut(LDM_CC, breaks = c(0, 0.21, 0.50, 10000, Inf)))

dat_num_max <- CC_LDM_FINAL %>% filter(LDM_CC != Inf) %>% select(LDM_CC) %>% max()
CC_LDM_FINAL <- CC_LDM_FINAL %>% mutate(LDM_CC_NEW = ifelse(LDM_CC==Inf, dat_num_max+1, LDM_CC))
CC_LDM_FINAL$Seed <- as.factor(CC_LDM_FINAL$Seed)

########### Plotting of Results ###########
## CC Plot ##
group.colors <- c(Missing="red",Present= "black")
axis_label_size <- 7
ylabel <- "Absolute of\n Log Disparate Value"
xlabel <- "CC Sample Size"

# text <- "Examining Absolute Disparate Value with varying CC Size and Seed"
# tgrob <- text_grob(text, size=15)
# plot_0 <- as_ggplot(tgrob) + theme(plot.margin = margin(0,0,0,0,"cm"))

#Gender
var.name <- "Missing_Subgroups"
G1 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='Female') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) + 
  #ylim(0, dat_num_max+1.5) + 
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw()+
  theme(axis.title=element_text(size = axis_label_size)) + 
  labs(x=xlabel, y = ylabel, caption = "Gender: Female") + 
  scale_color_manual(values=group.colors) 

G2 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>% 
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='Male') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() + 
  theme(axis.title=element_text(size = axis_label_size)) + 
  labs(x=xlabel, y = ylabel, caption = "Gender: Male") + 
  scale_color_manual(values=group.colors) 

#Ages
A1 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>% 
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='40-59') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() + 
  theme(axis.title=element_text(size = axis_label_size)) + 
  labs(x=xlabel, y = ylabel, caption = "Age Group: 40-59") + 
  scale_color_manual(values=group.colors)

A2 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>%
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='59+') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) + 
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() + 
  theme(axis.title=element_text(size = axis_label_size)) + 
  labs(x=xlabel, y = ylabel, caption = "Age Group: 59+") + 
  scale_color_manual(values=group.colors)

#Races
R1 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>% 
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='NH White') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) + 
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() + 
  theme(axis.title=element_text(size = axis_label_size)) + 
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: NH White") + 
  scale_color_manual(values=group.colors)

R2 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>% 
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='NH Black') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) + 
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() + 
  theme(axis.title=element_text(size = axis_label_size)) + 
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: NH Black") + 
  scale_color_manual(values=group.colors)

R3 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>% 
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='NH Asian') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) +
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() + 
  theme(axis.title=element_text(size = axis_label_size)) + 
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: NH Asian") + 
  scale_color_manual(values=group.colors)

R4 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>% 
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='Hispanic') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) + 
  #ylim(0, dat_num_max+1.5) +
  geom_point() +
  geom_hline(yintercept=0.22, linetype="dashed", color = "green") +
  theme_bw() + 
  theme(axis.title=element_text(size = axis_label_size)) + 
  labs(x=xlabel, y = ylabel, caption = "Race/Ethnicity: Hispanic") + 
  scale_color_manual(values=group.colors)

R5 <- CC_LDM_FINAL %>% mutate(Color = ifelse(LDM_CC==Inf, "Missing", "Present")) %>% 
  rename("Missing_Subgroups" = Color) %>%
  filter(Level=='Other') %>% 
  ggplot(aes(x=CC_Size, y=LDM_CC_NEW, color=Missing_Subgroups, shape=Seed)) + 
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
