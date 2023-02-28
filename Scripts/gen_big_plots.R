setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")
source("./Scripts/run_patt_scenario.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
directory <- "./Data/Results/M6/"
patt_filename <- "PATT_Summary_extrabias.csv"
ldm_filename <- "LDM_Summary_extrabias.csv"
Final_PATT_Summary <- read.csv(paste(directory, patt_filename, sep = ""))
Final_LDM_Summary <- read.csv(paste(directory, ldm_filename, sep = ""))

######## create mean PATT summary ########
confint <- 0.95
mean_patt_summary <- Final_PATT_Summary %>% 
  group_by(CC_Size) %>% 
  summarise(mean(TA_CC), 
            TA_CC_Low = confidence_interval(TA_CC, confint)[1], 
            TA_CC_High = confidence_interval(TA_CC, confint)[2],
            mean(TA_HC_Prop), 
            TA_HC_Prop_Low = confidence_interval(TA_HC_Prop, confint)[1], 
            TA_HC_Prop_High = confidence_interval(TA_HC_Prop, confint)[2],
            mean(TA_HC_IPF), 
            TA_HC_IPF_Low = confidence_interval(TA_HC_IPF, confint)[1], 
            TA_HC_IPF_High = confidence_interval(TA_HC_IPF, confint)[2], 
            mean(TA_HC_Both), 
            TA_HC_Both_Low = confidence_interval(TA_HC_Both, confint)[1], 
            TA_HC_Both_High = confidence_interval(TA_HC_Both, confint)[2]) %>%
  rename(TA_CC = "mean(TA_CC)",
         TA_HC_Prop = "mean(TA_HC_Prop)",
         TA_HC_IPF = "mean(TA_HC_IPF)",
         TA_HC_Both = "mean(TA_HC_Both)")

######### create mean LDM summary ##########
Inf_replace <- max(Final_LDM_Summary %>% mutate_all(function(x) ifelse(is.infinite(x), -1, x)) %>% 
                     select(LDM_CC:LDM_HC_Prop_IPF)) + 0.5
LDM_Summary_Inf <- Final_LDM_Summary
LDM_Summary_Inf[sapply(LDM_Summary_Inf, is.infinite)] <- Inf_replace

mean_ldm_summary <- LDM_Summary_Inf %>%
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

#### Final Plot ######

##################### PATT #######################

#ylim min and max point generation
min_point <- min(mean_patt_summary %>% 
                   replace(is.na(mean_patt_summary), 0) %>% 
                   select(TA_CC:TA_HC_Both_High)) - 0.25
max_point <- max(mean_patt_summary %>% 
                   replace(is.na(mean_patt_summary), -10000) %>%
                   select(TA_CC:TA_HC_Both_High)) + 0.25

sgt <- -0.314
gt <- -0.239

#plot only CC

get_patt_plot <- function(mean_patt_summary, pop_label, col, col_low, col_high, caption){
  p <- ggplot(mean_patt_summary, aes(y=CC_Size)) + 
    geom_point(aes(x=!!sym(col), colour=pop_label)) + 
    geom_errorbarh(aes(xmin=!!sym(col_low), xmax=!!sym(col_high), colour=pop_label, height=100)) +
    xlim(min_point, max_point) + ylim(-100,1200) + 
    geom_vline(aes(xintercept=sgt, color = "SPRINT"), linetype="dashed") + 
    geom_vline(aes(xintercept=gt, color = "NHANES"), linetype="dashed") + 
    geom_text(aes(x=sgt-0.05, y=600, label = as.character(sgt), colour = "SPRINT"), size = 3) +
    geom_text(aes(x=gt+0.05,  y=600, label = as.character(gt),  colour = "NHANES"), size = 3) +
    theme_classic()+
    labs(y="CC_Size", x="PATT Estimate", caption = caption) +
    scale_color_manual(name='Population', values=c("HC + Propensity + IPF"="black", "SPRINT"="#FF3300", "NHANES"="#339900"))
  return(p)
}

col_label = "TA_CC"
p1 <- get_patt_plot(mean_patt_summary, 
                    pop_label = "Only CC", 
                    col = col_label, 
                    caption = "TA + CC", 
                    col_low = paste(col_label, "_Low", sep=""), 
                    col_high = paste(col_label, "_High", sep=""))
p1

col_label = "TA_HC_Prop"
p2 <- get_patt_plot(mean_patt_summary, 
                    pop_label = "HC + Propensity", 
                    col = col_label, 
                    caption = "TA + HC + Propensity", 
                    col_low = paste(col_label, "_Low", sep=""), 
                    col_high = paste(col_label, "_High", sep=""))
p2

col_label = "TA_HC_IPF"
p3 <- get_patt_plot(mean_patt_summary, 
                    pop_label = "HC + IPF", 
                    col = col_label, 
                    caption = "TA + HC + IPF", 
                    col_low = paste(col_label, "_Low", sep=""), 
                    col_high = paste(col_label, "_High", sep=""))
p3

col_label = "TA_HC_Both"
p4 <- get_patt_plot(mean_patt_summary, 
                    pop_label = "HC + Propensity + IPF", 
                    col = col_label, 
                    caption = "TA + HC + Propensity + IPF", 
                    col_low = paste(col_label, "_Low", sep=""), 
                    col_high = paste(col_label, "_High", sep=""))
p4

##################### Equity #######################

#show LDM just for women
#Final_LDM_Summary <- Final_LDM_Summary %>% filter(CC_Size==cc_size)

get_ldm_plot <- function(mean_ldm_summary, var, level, pop_label, col, col_low, col_high){
  min_point <- min(ungroup(mean_ldm_summary) %>% filter(Level==level) %>% 
                    select(LDM_CC:LDM_HC_Prop_IPF_High))
  max_point <- max(ungroup(mean_ldm_summary) %>% filter(Level==level) %>% 
                    select(LDM_CC:LDM_HC_Prop_IPF_High))
  cap <- paste(var,":",level)
  
  p <- mean_ldm_summary %>% filter(Level==level) %>%
    ggplot(aes(y=CC_Size)) +
    geom_point(aes(x=!!sym(col), colour=pop_label)) +
    geom_errorbarh(aes(xmin=!!sym(col_low), xmax=!!sym(col_high), 
                       colour=pop_label, height=100)) +
    xlim(min_point, max_point) + ylim(-100, 1200) + 
    geom_vline(aes(xintercept=0.22, color = "Equity Range"), linetype="dashed") + 
    geom_text(aes(x=0.1, y=1150, label = as.character(0.22), 
                  colour = "Equity Range"), size = 3) +
    theme_classic()+
    labs(y="CC_Size", x="Median LDM Estimate", caption = cap) +
    scale_color_manual(name="Population",
                       values=c("HC + Propensity + IPF"="black", "Equity Range"="#339900"))
    
  return (p)
}

col_label = "LDM_CC"
l1 <- get_ldm_plot(mean_ldm_summary, var = "Gender", level = "Female", 
                   pop_label = "Only CC", 
                   col = col_label, 
                   col_low = paste(col_label, "_Low", sep=""), 
                   col_high = paste(col_label, "_High", sep=""))
l1

col_label = "LDM_HC_Prop"
l2 <- get_ldm_plot(mean_ldm_summary, var = "Gender", level = "Female", 
                   pop_label = "HC + Propensity", 
                   col = col_label, 
                   col_low = paste(col_label, "_Low", sep=""), 
                   col_high = paste(col_label, "_High", sep=""))
l2

col_label = "LDM_HC_IPF"
l3 <- get_ldm_plot(mean_ldm_summary, var = "Gender", level = "Female", 
                   pop_label = "HC + IPF", 
                   col = col_label, 
                   col_low = paste(col_label, "_Low", sep=""), 
                   col_high = paste(col_label, "_High", sep=""))
l3

col_label = "LDM_HC_Prop_IPF"
l4 <- get_ldm_plot(mean_ldm_summary, var = "Gender", level = "Female", 
                   pop_label = "HC + Propensity + IPF", 
                   col = col_label, 
                   col_low = paste(col_label, "_Low", sep=""), 
                   col_high = paste(col_label, "_High", sep=""))
l4

ggarrange(p1, l1, ncol=2,
          common.legend = TRUE, legend = "bottom")

