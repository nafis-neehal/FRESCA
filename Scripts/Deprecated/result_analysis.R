setwd("/Users/nafisneehal/ESCA")

suppressPackageStartupMessages({
  library(dplyr)
  library(Metrics)
})

ground <- read.csv("./Data/Results/PATT_ground.csv")
m1 <- read.csv("./Data/Results/PATT_m1.csv")
m2 <- read.csv("./Data/Results/PATT_m2.csv")
m3 <- read.csv("./Data/Results/PATT_m3.csv")

merge <- as.data.frame(do.call("cbind", c(ground %>% select(PATT_ground), m1 %>% select(PATT_m1), m2 %>% select(PATT_m2), 
                   m3 %>% select(PATT_m3))))

long_merge <- melt(merge)

#### Boxplot ####
h = -14.80
ggplot(long_merge, aes(x = variable, y = value)) + theme_classic() +
  geom_boxplot(fill='#A4A4A4', color="black") + 
  geom_hline(yintercept = h, color = "red", linetype = "dashed") +
  geom_text(aes(0, h, label=paste("SATT =", h), vjust = -1, hjust = -0.5))
  #scale_y_continuous(breaks = sort(c(seq(min(long_merge$value), max(long_merge$value), length.out=2), h)))
  
#### Statistical Tests ####
t.test(merge$PATT_ground, merge$PATT_m3)

