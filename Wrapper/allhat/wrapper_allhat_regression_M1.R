#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

### Model 1: Run a regression model to predict SqDeviation from target PHR: One model for one outcome ###
# X = TA Size, CC Size, SC Size, Controls (encoded 1: Only CC, 2: Only Prop, 3: Only Equity, 4: Both), SMD Value, CLD Value
###

#1. Load PHR Data
p1 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/HF_full_4000_0.csv")
p2 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/HF_full_4000_500.csv")
p3 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/HF_full_4000_1000.csv")
p4 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/HF_full_4000_2000.csv")
p5 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/HF_full_4000_4000.csv")

# Format the data in DF = <Seed, TA, CC, SC, Controls, SqDeviation>
target_phr <- 1.36

p1$TA_Size <- 4000
p1$CC_Size <- 0
p1$SC_Size <- p1$TA_Size - p1$CC_Size

p2$TA_Size <- 4000
p2$CC_Size <- 500
p2$SC_Size <- p2$TA_Size - p2$CC_Size

p3$TA_Size <- 4000
p3$CC_Size <- 1000
p3$SC_Size <- p3$TA_Size - p3$CC_Size

p4$TA_Size <- 4000
p4$CC_Size <- 2000
p4$SC_Size <- p4$TA_Size - p4$CC_Size

p5$TA_Size <- 4000
p5$CC_Size <- 4000
p5$SC_Size <- p5$TA_Size - p5$CC_Size

p <- do.call("rbind", list(p1,p2, p3, p4, p5))
#p$SYN_only <- NULL
#p <- p%>% filter(CC_Size!=0)

#na.replace(-1.00) %>%

phr_df <- p %>% 
  filter(combination=="TR1" & CC_Size!=0) %>% 
  select(!combination) %>%
  pivot_longer(!c(seed, TA_Size, CC_Size, SC_Size), names_to = "Controls", values_to = "PHR") %>%
  mutate(Deviation = (target_phr-PHR),
         SqDeviation = (target_phr-PHR)**2) %>% #<<<<--------------------------------
  filter(Controls!="SYN_only") %>%
  select(!PHR) %>%
  rename(Seed = seed) %>%
  arrange(Seed, Controls)

remove(p1,p2,p3,p4,p5)

#2. Load Equity Data
e1 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/EQUITY_New_full_4000_0.csv")
e2 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/EQUITY_New_full_4000_500.csv")
e3 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/EQUITY_New_full_4000_1000.csv")
e4 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/EQUITY_New_full_4000_2000.csv")
e5 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/EQUITY_New_full_4000_4000.csv")

e1$TA_Size <- 4000
e1$CC_Size <- 0
e1$SC_Size <- e1$TA_Size - e1$CC_Size

e2$TA_Size <- 4000
e2$CC_Size <- 500
e2$SC_Size <- e2$TA_Size - e2$CC_Size

e3$TA_Size <- 4000
e3$CC_Size <- 1000
e3$SC_Size <- e3$TA_Size - e3$CC_Size

e4$TA_Size <- 4000
e4$CC_Size <- 2000
e4$SC_Size <- e4$TA_Size - e4$CC_Size

e5$TA_Size <- 4000
e5$CC_Size <- 4000
e5$SC_Size <- e5$TA_Size - e5$CC_Size

e <- do.call("rbind", list(e1,e2,e3,e4,e5))
is.na(e)<-sapply(e, is.infinite)
#e <- na.replace(e, -1)

equity_df <- e %>% 
  filter(combination=="TR1" & CC_Size!=0) %>%
  select(!c(combination, Background_Rate,Observed_Rate)) %>%
  group_by(Population, TA_Size, CC_Size, SC_Size, Seed) %>% 
  summarise(CLD = median(LD)) %>%
  select(Seed, TA_Size, CC_Size, SC_Size, Population, CLD) %>%
  rename(Controls=Population) %>%
  filter(Controls!="Adj_TA") %>%
  mutate(Controls = recode(Controls, 'Both'='BOTH', 'Only_CC'='CC_only', 
                           'Only_Prop'='PROP_only', 'Only_Equity'='EQUITY_only')) %>%
  arrange(Seed, Controls)

remove(e1,e2,e3,e4,e5)

#3. Load SMD Data
s1 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/ASMD_new_weighted_4000_0.csv")
s2 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/ASMD_new_weighted_4000_500.csv")
s3 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/ASMD_new_weighted_4000_1000.csv")
s4 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/ASMD_new_weighted_4000_2000.csv")
s5 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/ASMD_new_weighted_4000_4000.csv")

s1$TA_Size <- 4000
s1$CC_Size <- 0
s1$SC_Size <- s1$TA_Size - s1$CC_Size

s2$TA_Size <- 4000
s2$CC_Size <- 500
s2$SC_Size <- s2$TA_Size - s2$CC_Size

s3$TA_Size <- 4000
s3$CC_Size <- 1000
s3$SC_Size <- s3$TA_Size - s3$CC_Size

s4$TA_Size <- 4000
s4$CC_Size <- 2000
s4$SC_Size <- s4$TA_Size - s4$CC_Size

s5$TA_Size <- 4000
s5$CC_Size <- 4000
s5$SC_Size <- s5$TA_Size - s5$CC_Size

s <- do.call("rbind", list(s1,s2,s3,s4,s5))

smd_df <- s %>% 
  filter(combination=="TR1" & CC_Size!=0) %>%
  select(!combination) %>% 
  group_by(Controls, TA_Size, CC_Size, SC_Size, Seed) %>% 
  summarise(SMD = mean(SMD**2)) %>%
  filter(Controls!='only_sc' && Controls!='only_biased_ec') %>%
  mutate(Controls = recode(Controls, 'both'='BOTH', 'only_cc'='CC_only', 
                           'only_prop'='PROP_only', 'only_equity'='EQUITY_only')) %>%
  select(Seed, TA_Size, CC_Size, SC_Size, Controls, SMD) %>%
  arrange(Seed, Controls)
  
remove(s1,s2,s3,s4,s5)

#4. Merge PHR, Equity and SMD dataframes to make the final dataframe
temp_df <- left_join(equity_df, smd_df)
final_df <- left_join(temp_df, phr_df)
copy_final_df <- final_df

#encode Controls to numeric values (Both=1, CC_Only=2, EQUIY_Only=3, PROP_Only=4)
#final_df$Controls <- as.numeric(factor(final_df$Controls, levels = unique(final_df$Controls)))
final_df$Controls <- factor(final_df$Controls)

final_df <- na.omit(final_df)

#write.csv(final_df, file = "./Data/Results/M6/ALLHAT_Results/Regression/final_df.csv", row.names = F)

#5.1. Run a regression model 
# Model1 <- lm(SqDeviation ~ TA_Size + CC_Size + SC_Size + Controls + CLD + SMD, data=final_df)
# Model1 <- lm(SqDeviation ~ CC_Size + SC_Size + Controls + CLD, data=final_df) #<----------
# Model1 <- lm(SqDeviation ~ CC_Size + SC_Size + Controls + SMD, data=final_df) #<----------
# Model1 <- randomForest(SqDeviation ~ TA_Size + CC_Size + SC_Size + Controls + CLD + SMD, data=final_df)
# Model1 <- lm(SqDeviation ~ Controls + CC_Size +  CLD + SMD + CLD*SMD + CC_Size*SMD + CC_Size*CLD + CC_Size*CLD*SMD, data = final_df)
Model1 <- lm(Deviation ~ Controls + CC_Size + I(SMD^2) + CLD, data=final_df)
summary(Model1)
res <- summary(Model1)

############################## Detour #########################

# pl1<- ggplot(final_df, aes(x = CLD, y = SqDeviation)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ poly(x, 1), se = T) +
#   geom_vline(xintercept=0.22, color="red") + 
#   theme_classic() +
#   labs(x = "CLD", y = "Predicted PHR SqDeviation", title = paste("SqDeviation = CLD^1")) + 
#   theme(aspect.ratio = 1)
# 
# 
# pl2<- ggplot(final_df, aes(x = CLD, y = SqDeviation)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = T) +
#   geom_vline(xintercept=0.22, color="red") + 
#   theme_classic() +
#   labs(x = "CLD", y = "Predicted PHR SqDeviation", title = paste("SqDeviation = CLD^2")) + 
#   theme(aspect.ratio = 1) 
# 
# pl3<- ggplot(final_df, aes(x = CLD, y = SqDeviation)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ poly(x, 3), se = T) +
#   geom_vline(xintercept=0.22, color="red") + 
#   theme_classic() +
#   labs(x = "CLD", y = "Predicted PHR SqDeviation", title = paste("SqDeviation = CLD^3")) + 
#   theme(aspect.ratio = 1) 
# 
# pl4<- ggplot(final_df, aes(x = CLD, y = SqDeviation)) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ poly(x, 4), se = T) +
#   geom_vline(xintercept=0.22, color="red") + 
#   theme_classic() +
#   labs(x = "CLD", y = "Predicted PHR SqDeviation", title = paste("SqDeviation = CLD^4")) + 
#   theme(aspect.ratio = 1) 
# 
# 
# combined_layout <- grid.arrange(grobs = list(pl1, pl2, pl3, pl4), ncol=4)
# print(combined_layout)

########################################################


#5.2. Model Performance (Training Error)
predicted_values <- predict(Model1, final_df %>% select(TA_Size, CC_Size, SC_Size, Controls, CLD, SMD))
final_df$Predicted <- predicted_values
mae <- mean(abs(final_df$SqDeviation - final_df$Predicted), na.rm = T)
mse <- mean((final_df$SqDeviation - final_df$Predicted)^2, na.rm = T)
cat("Mean Absolute Error (MAE):", mae, "\n")
cat("Mean Squared Error (MSE):", mse, "\n")




#6. Plot CC Size VS CLD
# generate_lm_predictions_varying_CLD <- function(model, cc_range) {
#   
#   # Check if input_range is a list
#   if (!is.list(cc_range)) {
#     stop("Input range must be a list of numeric vectors")
#   }
#   
#   CLD_Range <- c(0.01, 0.05, 0.22, 0.51, 0.80)
#   
#   plot_df <- as.data.frame(matrix(nrow=0, ncol=3))
#   col_names = c("Input", "CLD", "Predicted")
#   colnames(plot_df) <- col_names
#   
#   for (i in 1:length(CLD_Range)){
#   
#     # Predict values for the input range
#     new_data <- data.frame(TA_Size=4000, CC_Size=cc_range, Controls=1, SMD=0.10, CLD=CLD_Range[[i]])
#     new_data <- new_data %>% rename(CC_Size=x)
#     new_data$SC_Size <- new_data$TA_Size - new_data$CC_Size
#     # new_data <- new_data %>% select(TA_Size, CC_Size, SC_Size, Controls, CLD, SMD)
#     new_data <- new_data %>% select(CC_Size, SC_Size, Controls, CLD)
#     
#     predicted_values <- predict(model, new_data)
#     
#     
#     # Create a data frame for plotting
#     plot_data <- data.frame(Input = unlist(cc_range), CLD=as.character(CLD_Range[[i]]), Predicted = predicted_values)
#     plot_df <- rbind(plot_df, plot_data)
#   }
#   
#   return(plot_df)
#   
# }
# 
# plot_data_CLD <- generate_lm_predictions_varying_CLD(Model1, cc_range = list(x = seq(0, 4000, by = 100)))
# 
# gf_plot_CLD <- ggplot(plot_data_CLD,  aes(x = Input,y = Predicted, color = CLD)) +  
#   geom_line() +
#   labs(x = "CC Size", y = "Predicted PHR SqDeviation") + 
#   ggtitle("LM Model Predictions") + 
#   theme_classic() + 
#   theme(legend.title = element_text(),
#         legend.spacing.y = unit(0, "mm"), 
#         panel.border = element_rect(colour = "black", fill=NA),
#         aspect.ratio = 1, 
#         axis.text = element_text(colour = 1, size = 12),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   NULL
# 
# gf_plot_CLD

#7. Plot CC Size VS SMD

# generate_lm_predictions_varying_SMD <- function(model, cc_range) {
#   
#   # Check if input_range is a list
#   if (!is.list(cc_range)) {
#     stop("Input range must be a list of numeric vectors")
#   }
#   
#   SMD_Range <- c(0.01, 0.10, 0.20, 0.50, 1.00)
#   
#   plot_df <- as.data.frame(matrix(nrow=0, ncol=3))
#   col_names = c("Input", "SMD", "Predicted")
#   colnames(plot_df) <- col_names
#   
#   for (i in 1:length(SMD_Range)){
#     
#     # Predict values for the input range
#     new_data <- data.frame(TA_Size=4000, CC_Size=cc_range, Controls=1, SMD=SMD_Range[[i]], CLD=0.05)
#     new_data <- new_data %>% rename(CC_Size=x)
#     new_data$SC_Size <- new_data$TA_Size - new_data$CC_Size
#     new_data <- new_data %>% select(TA_Size, CC_Size, SC_Size, Controls, CLD, SMD)
#     
#     predicted_values <- predict(model, new_data)
#     
#     # Create a data frame for plotting
#     plot_data <- data.frame(Input = unlist(cc_range), SMD=as.character(SMD_Range[[i]]), Predicted = predicted_values)
#     plot_df <- rbind(plot_df, plot_data)
#   }
#   
#   return(plot_df)
#   
# }
# 
# plot_data_SMD <- generate_lm_predictions_varying_SMD(Model1, cc_range = list(x = seq(0, 4000, by = 100)))
# 
# gf_plot_SMD <- ggplot(plot_data_SMD,  aes(x = Input,y = Predicted, color = SMD)) +  
#   geom_line() +
#   labs(x = "CC Size", y = "Predicted PHR SqDeviation") + 
#   ggtitle("LM Model Predictions") + 
#   theme_classic() + 
#   theme(legend.title = element_text(),
#         legend.spacing.y = unit(0, "mm"), 
#         panel.border = element_rect(colour = "black", fill=NA),
#         aspect.ratio = 1, 
#         axis.text = element_text(colour = 1, size = 12),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   NULL
# 
# gf_plot_SMD

#8. Try 3D Plot
CC_Range <- seq(0, 4000, by = 500)
# SMD_Range <- seq(0, 1, 0.20)
# CLD_Range <- seq(0, 1, 0.20)
# SMD_Range <- c(0, 0.05, 0.10, 0.20)
# CLD_Range <- c(0, 0.10, 0.22, 0.51)
SMD_Range <- c(0, 0.10, 0.20)
CLD_Range <- c(0.03, 0.10, 0.22, 0.75)
# SMD_Range <- c(0.01, 0.10, 0.20, 0.50, 1.00)
# CLD_Range <- c(0.01, 0.05, 0.22, 0.51, 0.80)
df <- expand.grid(CC_Range, SMD_Range, CLD_Range)
colnames(df) <- c("CC_Size", "SMD", "CLD")
df$TA_Size <- 4000
df$Controls <- 1
df$SC_Size <- df$TA_Size - df$CC_Size 
df <- df %>% select(TA_Size, CC_Size, SC_Size, Controls, CLD, SMD)
# df <- df %>% select(CC_Size, SC_Size, Controls, CLD)
# df <- df %>% select(CC_Size, SC_Size, Controls, SMD)

new_predicted_value <- predict(Model1, df)
df$SqDeviation <- new_predicted_value

df$CLD <- as.character(df$CLD)
df$SMD <- as.character(df$SMD)

plot_ly(df, x = ~CC_Size, y = ~SMD, z = ~CLD,  type = "scatter3d", mode="markers",
        marker=list(color = ~SqDeviation, colorscale="RdBu", size=15, showscale=T,
                    colorbar=list(orientation="h",
                                  len=0.5,
                                  title="SqDeviation"))) %>%
  layout(title = "3D Scatter Plot of PHR SqDeviation varying with CC_Size, SMD, and CLD",
         margin = list(r=5, t=-5, b=5, l=5),
         scene = list(xaxis = list(showgrid=T, title = "CC Size"),
                      yaxis = list(showgrid=T, title = "Standardized Mean Difference"),
                      zaxis = list(showgrid=T, title = "Cohort Log Disparity")),
         paper_bgcolor = "#FFFFFF",
         plot_bgcolor = "#FFFFFF")

# plot_ly(df, x = ~CC_Size, y = ~SMD, z = ~SqDeviation,  type = "scatter3d", mode="markers",
#         marker=list(color = ~SqDeviation, colorscale="RdBu", size=15, showscale=T, 
#                     colorbar=list(orientation="h", 
#                                   len=0.5,
#                                   title="SqDeviation"))) %>%
#   layout(title = "3D Scatter Plot of PHR SqDeviation varying with CC_Size, and SMD",
#          margin = list(r=5, t=-5, b=5, l=5),
#          scene = list(xaxis = list(showgrid=T, title = "CC Size"),
#                       yaxis = list(showgrid=T, title = "Standardized Mean Difference"),
#                       zaxis = list(showgrid=T, title = "PHR SqDeviation")),
#          paper_bgcolor = "#FFFFFF",
#          plot_bgcolor = "#FFFFFF") 
  
#tbl <- final_df %>% filter(Controls==1)



