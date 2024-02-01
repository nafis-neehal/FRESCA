setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")
source("./Scripts/run_patt_scenario.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

directory <- "./Data/Results/M6/"
patt_filename <- "PATT_Summary_medbias_genpkg.csv"
ldm_filename <- "LDM_Summary_extrabias.csv"
Final_PATT_Summary <- read.csv(paste(directory, patt_filename, sep = ""))
Final_LDM_Summary <- read.csv(paste(directory, ldm_filename, sep = ""))

Final_PATT_Summary <- Final_PATT_Summary %>%
  mutate(TA_CC = exp(TA_CC),
         TA_HC = exp(TA_HC),
         TA_HC_Prop = exp(TA_HC_Prop),
         TA_HC_IPF = exp(TA_HC_IPF),
         TA_HC_Both = exp(TA_HC_Both))


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
# min_point <- min(mean_patt_summary %>% 
#                    replace(is.na(mean_patt_summary), 0) %>% 
#                    select(TA_CC:TA_HC_Both_High)) - 0.25
min_point <- 0
max_point <- max(mean_patt_summary %>% 
                   replace(is.na(mean_patt_summary), -10000) %>%
                   select(TA_CC:TA_HC_Both_High)) + 0.05

sgt <- 0.745 #sample
gt <- 0.798 #population
gt_new <- 0.757 #iptw adjusted population

#plot only CC

# get_patt_plot <- function(mean_patt_summary, pop_label, col, col_low, col_high, caption){
#   p <- ggplot(mean_patt_summary, aes(y=CC_Size)) + 
#     geom_point(aes(x=!!sym(col), colour=pop_label)) + 
#     geom_errorbarh(aes(xmin=!!sym(col_low), xmax=!!sym(col_high), colour=pop_label, height=100)) +
#     xlim(min_point, max_point) + ylim(-100,1200) + 
#     geom_vline(aes(xintercept=sgt, color = "SPRINT"), linetype="dashed") + 
#     geom_vline(aes(xintercept=gt, color = "NHANES"), linetype="dashed") + 
#     geom_text(aes(x=sgt-0.05, y=600, label = as.character(sgt), colour = "SPRINT"), size = 3) +
#     geom_text(aes(x=gt+0.05,  y=600, label = as.character(gt),  colour = "NHANES"), size = 3) +
#     theme_classic()+
#     labs(y="CC_Size", x="PATT Estimate", caption = caption) +
#     scale_color_manual(name='Population', values=c("HC + IPF"="black", "SPRINT"="#FF3300", "NHANES"="#339900"))
#   return(p)
# }

get_patt_plot_2 <- function(mean_patt_summary, pop_label, col, col_low, col_high){
  p <- ggplot(mean_patt_summary, aes(x = CC_Size, y=!!sym(col))) + 
    geom_point(aes(colour=pop_label)) + 
    geom_errorbar(aes(ymin=!!sym(col_low), ymax=!!sym(col_high), colour=pop_label, width=50)) +
    ylim(min_point, max_point) + xlim(-100,1200) + 
    #geom_hline(aes(yintercept=sgt, color = "Sample HR"), linetype="dashed") + 
    geom_hline(aes(yintercept=gt_new, color = "Target Population HR"), linetype="dashed") + 
    theme_classic()+
    labs(y="Hazard Ratio") + #, title=paste("Hazard Ratio for", pop_label, sep = " ")
    scale_x_continuous("CC Size", labels = as.character(mean_patt_summary$CC_Size), breaks = mean_patt_summary$CC_Size) +
    scale_color_manual(name='Method', 
                       values=c("HC + Equity"="black", "Sample HR"="#FF3300", "Target Population HR"="#FF3300")) + 
    theme(axis.text = element_text(size = 25)) + 
    theme(axis.title = element_text(size = 25)) + 
    theme(legend.title = element_text(size = 25)) +
    theme(legend.text = element_text(size = 25)) +
    theme(legend.position = c(0.78, 0.30),
          legend.background = element_rect(fill = "white", color = "black"))+
    # geom_text(aes(y=0.35, x=0, angle=90),
    #           label = "No data for CC Size=0",
    #           size = 6) +
    NULL
  return(p)
}


col_label = "TA_CC"
p1 <- get_patt_plot_2(mean_patt_summary, #%>% filter(CC_Size != 0)
                    pop_label = "CC Only", 
                    col = col_label, 
                    col_low = paste(col_label, "_Low", sep=""), 
                    col_high = paste(col_label, "_High", sep=""))
p1

col_label = "TA_HC_Prop"
p2 <- get_patt_plot_2(mean_patt_summary, 
                    pop_label = "HC + Propensity", 
                    col = col_label, 
                    col_low = paste(col_label, "_Low", sep=""), 
                    col_high = paste(col_label, "_High", sep=""))
p2

col_label = "TA_HC_IPF"
p3 <- get_patt_plot_2(mean_patt_summary, 
                    pop_label = "HC + Equity", 
                    col = col_label, 
                    col_low = paste(col_label, "_Low", sep=""), 
                    col_high = paste(col_label, "_High", sep=""))
p3

col_label = "TA_HC_Both"
p4 <- get_patt_plot_2(mean_patt_summary, 
                    pop_label = "HC + Propensity + Equity", 
                    col = col_label,
                    col_low = paste(col_label, "_Low", sep=""), 
                    col_high = paste(col_label, "_High", sep=""))
p4

##################### Equity #######################

#show LDM just for women
#Final_LDM_Summary <- Final_LDM_Summary %>% filter(CC_Size==cc_size)

# get_ldm_plot <- function(mean_ldm_summary, var, level, pop_label, col, col_low, col_high){
#   min_point <- min(ungroup(mean_ldm_summary) %>% filter(Level==level) %>% 
#                     select(LDM_CC:LDM_HC_Prop_IPF_High))
#   max_point <- max(ungroup(mean_ldm_summary) %>% filter(Level==level) %>% 
#                     select(LDM_CC:LDM_HC_Prop_IPF_High))
#   cap <- paste(var,":",level)
#   
#   p <- mean_ldm_summary %>% filter(Level==level) %>%
#     ggplot(aes(y=CC_Size)) +
#     geom_point(aes(x=!!sym(col), colour=pop_label)) +
#     geom_errorbarh(aes(xmin=!!sym(col_low), xmax=!!sym(col_high), 
#                        colour=pop_label, height=100)) +
#     xlim(min_point, max_point) + ylim(-100, 1200) + 
#     geom_vline(aes(xintercept=0.22, color = "Equity Range"), linetype="dashed") + 
#     geom_text(aes(x=0.1, y=1150, label = as.character(0.22), 
#                   colour = "Equity Range"), size = 3) +
#     theme_classic()+
#     labs(y="CC_Size", x="Median LDM Estimate", caption = cap) +
#     scale_color_manual(name="Population",
#                        values=c("HC + Propensity + IPF"="black", "Equity Range"="#339900"))
#     
#   return (p)
# }

get_ldm_plot_2 <- function(mean_ldm_summary, var, level, pop_label, col, col_low, col_high){
  min_point <- min(ungroup(mean_ldm_summary) %>% filter(Level==level) %>% 
                     select(!!sym(col), !!sym(col_low), !!sym(col_high)))
  max_point <- max(ungroup(mean_ldm_summary) %>% filter(Level==level) %>% 
                     select(!!sym(col), !!sym(col_low), !!sym(col_high)))
  p <- mean_ldm_summary %>% filter(Level==level) %>%
    ggplot(aes(x=CC_Size, y=!!sym(col))) +
    #geom_ribbon(aes(ymin=0, ymax=0.22), alpha=0.2, fill="green", outline.type="upper") +
    geom_rect(aes(xmin=-Inf, xmax = Inf, ymin = 0, ymax = 0.22), fill='#90EE90', alpha = 0.2) +
    geom_point(aes(colour=pop_label)) +
    geom_errorbar(aes(ymin=!!sym(col_low), ymax=!!sym(col_high), 
                      colour=pop_label, width=75)) +
    #ylim(min(0, min_point), max(0.22, max_point)) + 
    ylim(0,1) + 
    geom_hline(aes(yintercept=0.22, color = "Equity Range"), alpha=0.2) +
    theme_classic() +
    labs(y="Log Disparity", x="CC Size") +
    # labs(title=paste("Log Disparity Value for", pop_label, sep = " ")) +
    scale_x_continuous("CC Size", labels = as.character(mean_patt_summary$CC_Size), breaks = mean_patt_summary$CC_Size) +
    scale_color_manual(name="Method", values=c("HC + Propensity + Equity"="black", "Equity Range"="#339900")) +
    theme(axis.text = element_text(size = 25)) +
    theme(axis.title = element_text(size = 25)) +
    theme(legend.title = element_text(size = 25)) +
    theme(legend.text = element_text(size = 25)) +
    theme(legend.position = c(0.78, 0.45),
          legend.background = element_rect(fill = "white", color = "black"))+
    geom_text(aes(x=500, y=0.18),
              label = "Equitable Region",
              size = 6) +
    # geom_text(aes(y=0.60, x=0, angle=90),
    #           label = "No data for CC Size=0",
    #           size = 6) +
    NULL
  return (p)
}



col_label <- "LDM_CC"
pop_label1 <- "CC Only"
l1 <- get_ldm_plot_2(mean_ldm_summary %>% filter(CC_Size!=0), 
                   var = "Gender", level = "Female", 
                   pop_label = pop_label1, 
                   col = col_label, 
                   col_low = paste(col_label, "_Low", sep=""), 
                   col_high = paste(col_label, "_High", sep=""))
l1

col_label <- "LDM_HC_Prop"
pop_label2 <- "HC + Propensity"
l2 <- get_ldm_plot_2(mean_ldm_summary, var = "Gender", level = "Female", 
                   pop_label = pop_label2, 
                   col = col_label, 
                   col_low = paste(col_label, "_Low", sep=""), 
                   col_high = paste(col_label, "_High", sep=""))
l2

col_label <- "LDM_HC_IPF"
pop_label3 <- "HC + Equity"
l3 <- get_ldm_plot_2(mean_ldm_summary, var = "Gender", level = "Female", 
                   pop_label = pop_label3, 
                   col = col_label, 
                   col_low = paste(col_label, "_Low", sep=""), 
                   col_high = paste(col_label, "_High", sep=""))
l3

col_label <- "LDM_HC_Prop_IPF"
pop_label4 <- "HC + Propensity + Equity" 
l4 <- get_ldm_plot_2(mean_ldm_summary, var = "Gender", level = "Female", 
                   pop_label = "HC + Propensity + Equity", 
                   col = col_label, 
                   col_low = paste(col_label, "_Low", sep=""), 
                   col_high = paste(col_label, "_High", sep=""))
l4

# plot1 <- ggarrange(NULL, p1, NULL, l1, NULL, widths = c(0.02, 1, 0.03, 1, 0.02), nrow=1, legend="bottom")
# plot1
# #annotate_figure(plot1, top = text_grob(paste("PATT Estimate and Log Disparity Value for", pop_label1, sep = " "), color = "black", face = "bold", size = 14))
# 
# plot2 <- ggarrange(NULL, p2, NULL, l2, NULL, widths = c(0.02, 1, 0.03, 1, 0.02), nrow=1, legend="bottom")
# plot2
# #annotate_figure(plot2, top = text_grob(paste("PATT Estimate and Log Disparity Value for", pop_label2, sep = " "), color = "black", face = "bold", size = 14))
# 
# plot3 <- ggarrange(NULL, p3, NULL, l3, NULL, widths = c(0.02, 1, 0.03, 1, 0.02), nrow=1, legend="bottom")
# plot3
# #annotate_figure(plot3, top = text_grob(paste("PATT Estimate and Log Disparity Value for", pop_label3, sep = " "), color = "black", face = "bold", size = 14))
# 
# plot4 <- ggarrange(NULL, p4, NULL, l4, NULL, widths = c(0.02, 1, 0.03, 1, 0.02), nrow=1, legend="bottom")
# plot4
# #annotate_figure(plot4, top = text_grob(paste("PATT Estimate and Log Disparity Value for", pop_label4, sep = " "), color = "black", face = "bold", size = 14))



