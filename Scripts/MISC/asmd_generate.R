setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Scripts/common.R")
source("./Modules/p_val_add.R")
source("./Modules/equity_metrics_miao.R")

directory <- "./Data/Results/M6/"
fextra <- read.csv(paste(directory,"ASMD_extrabias_NC.csv", sep=""))
fmed <- read.csv(paste(directory,"ASMD_medbias_NC.csv", sep=""))
fno <- read.csv(paste(directory,"ASMD_nobias_NC.csv", sep=""))


########## Generate ASMD Plot ##############

## Construct a data frame containing variable name and SMD from all methods
dataPlot <- data.frame(variable = fextra$variable,
                       Unbiased = fno$Unbiased,
                       Sample1  = fmed$Biased,
                       Sample2  = fextra$Biased)

## Create long-format data for ggplot2
dataPlotMelt <- melt(data          = dataPlot,
                     id.vars       = c("variable"),
                     variable.name = "Biased_Degree",
                     value.name    = "SMD")

## Order variable names by magnitude of SMD
varNames <- as.character(dataPlot$variable)[order(dataPlot$Sample1)]

## Order factor levels in the same order
dataPlotMelt$variable <- factor(dataPlotMelt$variable,
                                levels = varNames)

## Plot using ggplot2
ggplot(data = dataPlotMelt,
       mapping = aes(x = variable, y = SMD, group = Biased_Degree, color = Biased_Degree)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.1, color = "black", size = 0.1) +
  coord_flip() +
  theme_bw() +
  theme(legend.key = element_blank())








