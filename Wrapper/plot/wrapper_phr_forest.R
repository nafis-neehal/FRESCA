#### SOURCES ####
##### Note: Export PDF to 15.75x4.00 inches ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")
library(rlang)

p1 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_summary_4000_0.csv")
p2 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_summary_4000_500.csv")
p3 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_summary_4000_1000.csv")
p4 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_summary_4000_2000.csv")
p5 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_summary_4000_4000.csv")

# p1 <- read.csv("./Data/Results/M6/New_SPRINT_Results/phr_extrabias_simplematch_ipfw/Matchit3/PRIMARY_summary_2000_0.csv")
# p2 <- read.csv("./Data/Results/M6/New_SPRINT_Results/phr_extrabias_simplematch_ipfw/Matchit3/PRIMARY_summary_2000_500.csv")
# p3 <- read.csv("./Data/Results/M6/New_SPRINT_Results/phr_extrabias_simplematch_ipfw/Matchit3/PRIMARY_summary_2000_1000.csv")
# p4 <- read.csv("./Data/Results/M6/New_SPRINT_Results/phr_extrabias_simplematch_ipfw/Matchit3/PRIMARY_summary_2000_1500.csv")
# p5 <- read.csv("./Data/Results/M6/New_SPRINT_Results/phr_extrabias_simplematch_ipfw/Matchit3/PRIMARY_summary_2000_2000.csv")

p1$CC_Size <- 0
p2$CC_Size <- 500
p3$CC_Size <- 1000
p4$CC_Size <- 2000
p5$CC_Size <- 4000

# p1$CC_Size <- 0
# p2$CC_Size <- 500
# p3$CC_Size <- 1000
# p4$CC_Size <- 1500
# p5$CC_Size <- 2000

get_df_reshaped <- function(df, colstr){
  df1 <- df %>%
    pivot_wider(
      names_from = values,
      values_from = !!sym(colstr),
      names_prefix = "V_"
    )
  colnames(df1) <- sub("V_", "", colnames(df1))
  return(df1)
}

p <- do.call("rbind", list(p1,p2,p3,p4,p5))
p_new <- na.replace(p, 0)
p_new <- p_new %>% mutate(CC_Size = as.factor(CC_Size))

# Create a list to store plots
combined_plots <- list()

Combinations <- c("Only CC", "Propensity Matched", "Equity Adjusted", "Both")
####Target_Tr1 <- c(1.36, 1.34, 1.39)
Target_Tr1 <- c(1.38, 1.36, 1.41)
# Target_Tr1 <- c(0.79, 0.77, 0.82)


###### Trial 1 : Heart Failure #######

p_new_T1 <- p_new %>% filter(combination=='TR1') %>% select(CC_Size, values, CC_only)
p_new_T1_cc <- get_df_reshaped(p_new_T1, "CC_only")

p_new_T1 <- p_new %>% filter(combination=='TR1') %>% select(CC_Size, values, PROP_only)
p_new_T1_prop <- get_df_reshaped(p_new_T1, "PROP_only")

p_new_T1 <- p_new %>% filter(combination=='TR1') %>% select(CC_Size, values, EQUITY_only)
p_new_T1_equity <- get_df_reshaped(p_new_T1, "EQUITY_only")

p_new_T1 <- p_new %>% filter(combination=='TR1') %>% select(CC_Size, values, BOTH)
p_new_T1_both <- get_df_reshaped(p_new_T1, "BOTH")

forest_data_list <- list(p_new_T1_cc, p_new_T1_prop, p_new_T1_equity, p_new_T1_both)

#+ scale_y_continuous(breaks = forest_data_list[[i]]$CC_Size) + 

# Loop through each forest plot
for (i in 1:4) {
  forest_plot <- ggplot(forest_data_list[[i]], aes(x = mean, y = CC_Size)) +
    geom_point() + 
    geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.5) +
    geom_vline(aes(xintercept = Target_Tr1[1], colour="red"), linetype = "dashed") +
    geom_vline(aes(xintercept = Target_Tr1[2], colour="red")) +
    geom_vline(aes(xintercept = Target_Tr1[3], colour="red")) +
    labs(
      x = "Effect Size",
      y = "CC Size",
      title = paste("PHR (ALLHAT) - ", Combinations[i])
    ) +
    #xlim(0.50, 1.00)+
    xlim(1.00, 1.80)+
    theme_classic() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none",  # If no legend is needed
      plot.margin = margin(5, 20, 5, 5)  # Adjust plot margins for text positioning
    )
  if (i==1 || i==3){
    forest_plot <- forest_plot + 
      geom_text(aes(x=1.05, y=1),
                label = "No data",
                size = 3)
  }
  combined_plots[[i]] <- forest_plot
}

# Arrange plots using grid.arrange from gridExtra
combined_layout <- grid.arrange(grobs = combined_plots, ncol = 4)

# Display the combined layout
print(combined_layout)


