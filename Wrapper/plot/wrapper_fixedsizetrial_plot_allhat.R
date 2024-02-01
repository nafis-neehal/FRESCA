#### SOURCES ####
##### Note: Export PDF to 15.75x4.00 inches ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")


# p1 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_10k/Matchit3/HF_samecost_summary_5000_5000.csv")
# p2 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_10k/Matchit3/HF_samecost_summary_5500_4500.csv")
# p3 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_10k/Matchit3/HF_samecost_summary_6000_4000.csv")
# p4 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_10k/Matchit3/HF_samecost_summary_6500_3500.csv")
# p5 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_10k/Matchit3/HF_samecost_summary_7000_3000.csv")
# p6 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_10k/Matchit3/HF_samecost_summary_7500_2500.csv")
# p7 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_10k/Matchit3/HF_samecost_summary_8000_2000.csv")

# p1 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_8k/HF_samecost_summary_5000_3000.csv")
# p2 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_8k/HF_samecost_summary_5500_2500.csv")
# p3 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_8k/HF_samecost_summary_6000_2000.csv")
# p4 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_8k/HF_samecost_summary_6500_1500.csv")
# p5 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_8k/HF_samecost_summary_7000_1000.csv")
# p6 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_8k/HF_samecost_summary_7500_500.csv")
# p7 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_8k/HF_samecost_summary_8000_0.csv")

p1 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_4k/Matchit3/HF_samecost_summary_3500_500.csv")
p2 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_4k/Matchit3/HF_samecost_summary_3000_1000.csv")
p3 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_4k/Matchit3/HF_samecost_summary_2500_1500.csv")
p4 <- read.csv("./Data/Results/M6/ALLHAT_Results/Samecost/trial_size_4k/Matchit3/HF_samecost_summary_2000_2000.csv")


# p1$TA_Size <- 5000
# p2$TA_Size <- 5500
# p3$TA_Size <- 6000
# p4$TA_Size <- 6500
# p5$TA_Size <- 7000
# p6$TA_Size <- 7500
# p7$TA_Size <- 8000
# p <- do.call("rbind", list(p1,p2,p3,p4,p5,p6,p7))

# p1$TA_Size <- 3500
# p2$TA_Size <- 3000
# p3$TA_Size <- 2500
# p4$TA_Size <- 2000

p1$SC_Size <- 3000
p2$SC_Size <- 2000
p3$SC_Size <- 1000
p4$SC_Size <- 0

p <- do.call("rbind", list(p1,p2,p3,p4))

p_new <- na.replace(p, 0)
p_new <- p_new %>% mutate(SC_Size = as.factor(SC_Size))

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

# Create a list to store plots
combined_plots2 <- list()
Combinations <- c("Both")
Target_Tr2 <- c(1.38, 1.36, 1.41)
p_new_T1 <- p_new %>% filter(combination=='TR1') %>% select(SC_Size, values, BOTH)
p_new_T1_both <- get_df_reshaped(p_new_T1, "BOTH")

forest_data_list <- list(p_new_T1_both)

for (i in 1:1) {
  forest_plot2 <- ggplot(forest_data_list[[i]], aes(x = mean, y = SC_Size)) +
    geom_point() + 
    geom_errorbarh(aes(xmin = lower_CI, xmax = upper_CI), height = 0.5) +
    geom_vline(aes(xintercept = Target_Tr2[1], colour="red"), linetype = "dashed") +
    geom_vline(aes(xintercept = Target_Tr2[2], colour="red")) +
    geom_vline(aes(xintercept = Target_Tr2[3], colour="red")) +
    labs(
      x = "Effect Size",
      y = "SC Size",
      title = paste("PHR (ALLHAT)")
      #title = paste("PHR (ALLHAT) - ", Combinations[i], "[Trial Size 4k]") # \\ CC Size =  (5000, 4500, 4000, ..., 2000) \\ SC Size = (0, 1000, 2000, ..., 6000)]")
    ) +
    #xlim(0.50, 1.00)+
    xlim(1.00, 1.80)+
    theme_classic() +
    theme(
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      aspect.ratio = 1,
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = "none"  # If no legend is needed
      #plot.margin = margin(5, 20, 5, 5)  # Adjust plot margins for text positioning
    )
  
  combined_plots2[[i]] <- forest_plot2
}

combined_layout <- grid.arrange(grobs = c(combined_plots2, combined_plots, NULL, NULL), ncol = 4)



