#### SOURCES ####
##### Note: Export PDF to 15.75x4.00 inches ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")


e <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_summary_4k.csv")
# e <- read.csv("./Data/Results/M6/New_SPRINT_Results/equity/Matchit3/EQUITY_summary.csv")

is.na(e)<-sapply(e, is.infinite)
p_new <- na.replace(e, -1)
p_new <- p_new %>% mutate(CC_Size = as.factor(CC_Size))

# Create a list to store plots
combined_plots <- list()

Combinations <- c("Only CC", "Propensity Matched", "Equity Adjusted", "Both")
Target_Tr1 <- c(0.00, 0.22)


###### Trial 1 : Heart Failure #######

p_new_T1_cc <- p_new %>% filter(combination=='TR1', Population=='Only_CC') %>% select(CC_Size, Median, CI_Lower, CI_Upper)
p_new_T1_prop <- p_new %>% filter(combination=='TR1', Population=='Only_Prop') %>% select(CC_Size, Median, CI_Lower, CI_Upper)
p_new_T1_equity <- p_new %>% filter(combination=='TR1', Population=='Only_Equity') %>% select(CC_Size, Median, CI_Lower, CI_Upper)
p_new_T1_both <- p_new %>% filter(combination=='TR1', Population=='Both') %>% select(CC_Size, Median, CI_Lower, CI_Upper)

forest_data_list <- list(p_new_T1_cc, p_new_T1_prop, p_new_T1_equity, p_new_T1_both)

#+ scale_y_continuous(breaks = forest_data_list[[i]]$CC_Size) + 

# Loop through each forest plot
for (i in 1:4) {
  forest_plot <- ggplot(forest_data_list[[i]], aes(x = Median, y = CC_Size)) +
    geom_rect(aes(xmin=0, xmax = 0.22, ymin = -Inf, ymax = Inf), fill='#90EE90', alpha = 0.2) +
    geom_point() + 
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0.5) +
    geom_vline(xintercept = Target_Tr1[1], colour="black", linetype="dashed") +
    geom_vline(xintercept = Target_Tr1[2], colour="black", linetype="dashed") +
    labs(
      x = "Effect Size",
      y = "CC Size",
      title = paste("Equity (ALLHAT) - ", Combinations[i])
    ) +
    xlim(0.00, 1.00)+
    theme_classic() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none",  # If no legend is needed
      plot.margin = margin(5, 20, 5, 5)  # Adjust plot margins for text positioning
    )+
    geom_text(aes(x=0.15, y=3, angle = 90),
              label = "Equitable Region",
              size = 3) +
    NULL
  if (i==1 || i==3){
    forest_plot <- forest_plot + 
      geom_text(aes(x=0.10, y=1),
                label = "No data",
                size = 3)
    }
  combined_plots[[i]] <- forest_plot
}

# Arrange plots using grid.arrange from gridExtra
combined_layout <- grid.arrange(grobs = combined_plots, ncol = 4)

# Display the combined layout
print(combined_layout)


