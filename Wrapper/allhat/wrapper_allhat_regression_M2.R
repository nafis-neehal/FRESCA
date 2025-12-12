#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
source("./Wrapper/load_data_formatted.R")

### Model 1: Run a regression model to predict SqDeviation from target PHR: One model for one outcome ###
# X = TA Size, CC Size, SC Size, Controls (encoded 1: Only CC, 2: Only Prop, 3: Only Equity, 4: Both), SMD Value, CLD Value
###

#1. Load PHR Data
p1 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_4000_0.csv")
p2 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_4000_500.csv")
p3 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_4000_1000.csv")
p4 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_4000_2000.csv")
p5 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_4000_4000.csv")
p6 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_5000_0.csv")
p7 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_5000_1000.csv")
p8 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_5000_2000.csv")
p9 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_5000_3000.csv")
p10 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_5000_4000.csv")
p11 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_5000_5000.csv")
p12 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_6000_0.csv")
p13 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_6000_1000.csv")
p14 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_6000_2000.csv")
p15 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_6000_3000.csv")
p16 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_6000_4000.csv")
p17 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_6000_5000.csv")
p18 <- read.csv("./Data/Results/M6/ALLHAT_Results/PHR_ipfw/Matchit3/HF_full_6000_6000.csv")


# Format the data in DF = <Seed, TA, CC, SC, Controls, SqDeviation>
target_phr <- 1.38

p1$TA_Size <- 4000
p1$CC_Size <- 0

p2$TA_Size <- 4000
p2$CC_Size <- 500

p3$TA_Size <- 4000
p3$CC_Size <- 1000

p4$TA_Size <- 4000
p4$CC_Size <- 2000

p5$TA_Size <- 4000
p5$CC_Size <- 4000

p6$TA_Size <- 5000
p6$CC_Size <- 0

p7$TA_Size <- 5000
p7$CC_Size <- 1000

p8$TA_Size <- 5000
p8$CC_Size <- 2000

p9$TA_Size <- 5000
p9$CC_Size <- 3000

p10$TA_Size <- 5000
p10$CC_Size <- 4000

p11$TA_Size <- 5000
p11$CC_Size <- 5000

p12$TA_Size <- 6000
p12$CC_Size <- 0

p13$TA_Size <- 6000
p13$CC_Size <- 1000

p14$TA_Size <- 6000
p14$CC_Size <- 2000

p15$TA_Size <- 6000
p15$CC_Size <- 3000

p16$TA_Size <- 6000
p16$CC_Size <- 4000

p17$TA_Size <- 6000
p17$CC_Size <- 5000

p18$TA_Size <- 6000
p18$CC_Size <- 6000

p <- do.call("rbind", list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18))
p$SC_Size <- p$TA_Size - p$CC_Size

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

remove(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18)

#2. Load Equity Data
e1 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_0.csv")
e2 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_500.csv")
e3 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_1000.csv")
e4 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_2000.csv")
e5 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_4000_4000.csv")

e6 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_5000_0.csv")
e7 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_5000_1000.csv")
e8 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_5000_2000.csv")
e9 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_5000_3000.csv")
e10 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_5000_4000.csv")
e11 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_5000_5000.csv")

e12 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_6000_0.csv")
e13 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_6000_1000.csv")
e14 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_6000_2000.csv")
e15 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_6000_3000.csv")
e16 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_6000_4000.csv")
e17 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_6000_5000.csv")
e18 <- read.csv("./Data/Results/M6/ALLHAT_Results/Equity/Matchit3/EQUITY_New_full_6000_6000.csv")


e1$TA_Size <- 4000
e1$CC_Size <- 0

e2$TA_Size <- 4000
e2$CC_Size <- 500

e3$TA_Size <- 4000
e3$CC_Size <- 1000

e4$TA_Size <- 4000
e4$CC_Size <- 2000

e5$TA_Size <- 4000
e5$CC_Size <- 4000

e6$TA_Size <- 5000
e6$CC_Size <- 0

e7$TA_Size <- 5000
e7$CC_Size <- 1000

e8$TA_Size <- 5000
e8$CC_Size <- 2000

e9$TA_Size <- 5000
e9$CC_Size <- 3000

e10$TA_Size <- 5000
e10$CC_Size <- 4000

e11$TA_Size <- 5000
e11$CC_Size <- 5000

e12$TA_Size <- 6000
e12$CC_Size <- 0

e13$TA_Size <- 6000
e13$CC_Size <- 1000

e14$TA_Size <- 6000
e14$CC_Size <- 2000

e15$TA_Size <- 6000
e15$CC_Size <- 3000

e16$TA_Size <- 6000
e16$CC_Size <- 4000

e17$TA_Size <- 6000
e17$CC_Size <- 5000

e18$TA_Size <- 6000
e18$CC_Size <- 6000


e <- do.call("rbind", list(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18))
is.na(e)<-sapply(e, is.infinite)
e$SC_Size <- e$TA_Size - e$CC_Size

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

remove(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18)

#3. Load SMD Data
s1 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_0.csv")
s2 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_500.csv")
s3 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_1000.csv")
s4 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_2000.csv")
s5 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_4000_4000.csv")

s6 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_5000_0.csv")
s7 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_5000_1000.csv")
s8 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_5000_2000.csv")
s9 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_5000_3000.csv")
s10 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_5000_4000.csv")
s11 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_5000_5000.csv")

s12 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_6000_0.csv")
s13 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_6000_1000.csv")
s14 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_6000_2000.csv")
s15 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_6000_3000.csv")
s16 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_6000_4000.csv")
s17 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_6000_5000.csv")
s18 <- read.csv("./Data/Results/M6/ALLHAT_Results/ASMD_Final/Matchit3/ASMD_new_weighted_6000_6000.csv")

s1$TA_Size <- 4000
s1$CC_Size <- 0

s2$TA_Size <- 4000
s2$CC_Size <- 500

s3$TA_Size <- 4000
s3$CC_Size <- 1000

s4$TA_Size <- 4000
s4$CC_Size <- 2000

s5$TA_Size <- 4000
s5$CC_Size <- 4000

s6$TA_Size <- 5000
s6$CC_Size <- 0

s7$TA_Size <- 5000
s7$CC_Size <- 1000

s8$TA_Size <- 5000
s8$CC_Size <- 2000

s9$TA_Size <- 5000
s9$CC_Size <- 3000

s10$TA_Size <- 5000
s10$CC_Size <- 4000

s11$TA_Size <- 5000
s11$CC_Size <- 5000

s12$TA_Size <- 6000
s12$CC_Size <- 0

s13$TA_Size <- 6000
s13$CC_Size <- 1000

s14$TA_Size <- 6000
s14$CC_Size <- 2000

s15$TA_Size <- 6000
s15$CC_Size <- 3000

s16$TA_Size <- 6000
s16$CC_Size <- 4000

s17$TA_Size <- 6000
s17$CC_Size <- 5000

s18$TA_Size <- 6000
s18$CC_Size <- 6000

s <- do.call("rbind", list(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18))
s$SC_Size <- s$TA_Size - s$CC_Size

smd_df <- s %>% 
  filter(combination=="TR1" & CC_Size!=0) %>%
  select(!combination) %>% 
  group_by(Controls, TA_Size, CC_Size, SC_Size, Seed) %>% 
  summarise(SMD = mean(SMD)) %>%
  filter(Controls!='only_sc' && Controls!='only_biased_ec') %>%
  mutate(Controls = recode(Controls, 'both'='BOTH', 'only_cc'='CC_only', 
                           'only_prop'='PROP_only', 'only_equity'='EQUITY_only')) %>%
  select(Seed, TA_Size, CC_Size, SC_Size, Controls, SMD) %>%
  arrange(Seed, Controls)

remove(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15, s16, s17, s18)

#4. Merge PHR, Equity and SMD dataframes to make the final dataframe
temp_df <- left_join(equity_df, smd_df)
final_df <- left_join(temp_df, phr_df)
copy_final_df <- final_df

#encode Controls to numeric values (Both=1, CC_Only=2, EQUIY_Only=3, PROP_Only=4)
#final_df$Controls <- as.numeric(factor(final_df$Controls, levels = unique(final_df$Controls)))
final_df[final_df=='BOTH'] <- 'FULL'
final_df$Controls <- factor(final_df$Controls)
final_df <- na.omit(final_df)

# write.csv(final_df, file = "./Data/Results/M6/ALLHAT_Results/Regression/final_df2.csv", row.names = F)

#5.1. Run a regression model 
Model1 <- lm(SqDeviation ~ TA_Size + CC_Size + Controls + SMD + CLD, data=final_df)
Model2 <- lm(SqDeviation ~ TA_Size + CC_Size + SMD + CLD, data = final_df)
Model3 <- lm(SqDeviation ~ TA_Size + CC_Size + Controls + CLD, data = final_df)
Model4 <- lm(SqDeviation ~ TA_Size + CC_Size + Controls + SMD, data = final_df)
Model5 <- lm(SqDeviation ~ TA_Size + CC_Size + Controls + SMD * CLD, data=final_df)
Model6 <- lm(SqDeviation ~ TA_Size + CC_Size + Controls * CLD + SMD, data=final_df)
Model7 <- lm(SqDeviation ~ TA_Size + CC_Size * CLD + Controls + SMD, data=final_df)


# Model1 <- lm(Deviation ~ TA_Size + CC_Size + Controls + SMD + CLD, data=final_df)
# Model2 <- lm(Deviation ~ TA_Size + CC_Size + SMD + CLD, data = final_df)
# Model3 <- lm(Deviation ~ TA_Size + CC_Size + Controls + CLD, data = final_df)
# Model4 <- lm(Deviation ~ TA_Size + CC_Size + Controls + SMD, data = final_df)
# Model5 <- lm(Deviation ~ TA_Size + CC_Size + Controls + SMD * CLD, data=final_df)
# Model6 <- lm(Deviation ~ TA_Size + CC_Size + Controls * CLD + SMD, data=final_df)
# Model7 <- lm(Deviation ~ TA_Size + CC_Size * CLD + Controls + SMD, data=final_df)

models <- list(Model1, Model2, Model3, Model4, Model5, Model6, Model7)
mod.names <- c('M1', 'M2', 'M3', 'M4', 'M5', 'M6', 'M7')

summary_glance <- as.data.frame(bind_rows(lapply(models, broom::glance)))
summary_glance$Model <- mod.names
final_summary_glance <- summary_glance %>% select(Model, r.squared, logLik, AIC, BIC)

aictab(cand.set = models, modnames = mod.names)


##### Model Assessment ########
library(jtools)

rescale_dataframe <- function(df, cols_to_rescale, min_val = 0, max_val = 1) {
  # Rescale specified columns between min_val and max_val
  for (col_name in cols_to_rescale) {
    if (is.numeric(df[[col_name]])) {
      min_col <- min(df[[col_name]], na.rm = TRUE)
      max_col <- max(df[[col_name]], na.rm = TRUE)
      
      if (min_col == max_col) {
        df[[col_name]] <- rep(min_val, length(df[[col_name]]))
      } else {
        df[[col_name]] <- ((df[[col_name]] - min_col) / (max_col - min_col)) * (max_val - min_val) + min_val
      }
    }
  }
  
  return(df)
}

final_df2 <- rescale_dataframe(final_df, c('TA_Size', 'CC_Size', 'SMD', 'CLD'), 
                               min_val = 0, max_val = 1)






#final_df <- final_df %>% mutate(TDev = ifelse(abs(Deviation)<=1, (abs(Deviation))**(1/3) , Deviation**2))


Model8 <- lm(SqDeviation ~ TA_Size + CC_Size + CLD + SMD  + Controls, data=final_df2)
summary(Model8)

summ(Model8, confint = T, digits = 5)

plot_summs(Model8, colors = "black", size=3) 

# hist(final_df$CLD)
# library(psych)
# pairs.panels(final_df %>% ungroup() %>% select(TA_Size, CC_Size, Controls, SMD, Deviation, SqDeviation))




########################## Outdated Codes ##########################
#5.2. Search grid with the best Model - Model 7

# TA_Range <- seq(2000, 7000, by = 1000)
# 
# grid_df <- as.data.frame(matrix(nrow = 0, ncol = 4))
# cols <- c("TA_Size",  "CC_Size", "CLD", "SMD")
# colnames(grid_df) <- cols
# 
# for (ta_size in TA_Range){
#   CC_Range <- seq(0, ta_size, by = 500)
#   SMD_Range <- c(0, 0.10, 0.20)
#   CLD_Range <- c(0, 0.10, 0.22)
#   df <- expand.grid(ta_size, CC_Range, CLD_Range, SMD_Range)
#   colnames(df) <- cols
#   grid_df <- rbind(grid_df, df)
# }
# 
# grid_df$Controls <- factor("BOTH")
# grid_df$SC_Size <- grid_df$TA_Size - grid_df$CC_Size
# grid_df <- grid_df %>% select(TA_Size, CC_Size, SC_Size, Controls, CLD, SMD)
# 
# new_predicted_value <- predict(Model7, grid_df)
# 
# #### min max scaling to avoid negative squared deviation prediction ####
# #new_scaled_predicted_value <- rescale(new_predicted_value, to = c(0, 0.5))
# #grid_df$SqDeviation <- abs(new_predicted_value)
# #grid_df$SqDeviation <- abs(new_predicted_value)
# grid_df$SqDeviation <- (new_predicted_value)**2
# 
# 
# 
# #5.3 Plot the best combination?
# # grid_df[which.min(grid_df$SqDeviation),]
# # test_df <- grid_df %>% filter(CLD==0, SMD==0.10)
# # plot_ly() %>%
# #   add_trace(data = test_df,  x=test_df$TA_Size, y=test_df$CC_Size, z=test_df$SqDeviation, type="mesh3d" )
# 
# plot_ly(grid_df, x = ~TA_Size, y = ~CC_Size, z = ~SMD,  type = "scatter3d", mode="markers",
#         marker=list(color = ~SqDeviation, colorscale="RdBu", size=15, showscale=T,
#                     colorbar=list(orientation="h",
#                                   len=0.5,
#                                   title="SqDeviation"))) %>%
#   layout(title = "3D Scatter Plot of PHR Squared Deviation varying with TA Size, CC Size, and SMD",
#          margin = list(r=5, t=-5, b=5, l=5),
#          scene = list(xaxis = list(showgrid=T, title = "TA Size"),
#                       yaxis = list(showgrid=T, title = "CC Size"),
#                       zaxis = list(showgrid=T, title = "SMD")),
#          paper_bgcolor = "#FFFFFF",
#          plot_bgcolor = "#FFFFFF")









