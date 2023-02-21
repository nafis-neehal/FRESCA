setwd("/Users/nafisneehal/ESCA")
source("./Scripts/common.R")

directory <- "./Data/Results/PATT_m6/CC_EC_Vary_Size/matchCCEC/summary_seed_results/"

dat_1 <- read.csv(paste(directory, "1500TA_vary_CC_EC_matchCCEC_seed10.csv", sep = ""))
dat_2 <- read.csv(paste(directory, "1500TA_vary_CC_2EC_matchCCEC_seed10.csv", sep = ""))
dat_1_5 <- read.csv(paste(directory, "1500TA_vary_CC_1_5EC_matchCCEC_seed10.csv", sep = ""))

TA_size <- 1500

dat1_slice <- dat_1 %>% select(CC_Sample_Size, PATT, LOWER_CF, UPPER_CF)
dat1_slice$Ratio <- "1"

dat2_slice <- dat_2 %>% select(CC_Sample_Size, PATT, LOWER_CF, UPPER_CF)
dat2_slice$Ratio <- "2"

dat1_5_slice <- dat_1_5 %>% select(CC_Sample_Size, PATT, LOWER_CF, UPPER_CF)
dat1_5_slice$Ratio <- "1.5"

final <- do.call("rbind", list(dat1_slice, dat2_slice, dat1_5_slice))

final1 <- final %>% filter(Ratio=="1")
final2 <- final %>% filter(Ratio=="2")
final1_5 <- final %>% filter(Ratio=="1.5")

# par(mar = c(1, 1, 1, 1))
# p1 <- errbar(final1$CC_Sample_Size, final1$PATT, final1$UPPER_CF, final1$LOWER_CF,
#        xlab = "CC Size", ylab = "PATT Estimation")
# p2 <- errbar(final2$CC_Sample_Size, final2$PATT, final2$UPPER_CF, final2$LOWER_CF,
#              xlab = "CC Size", ylab = "PATT Estimation")
# p1_5 <- errbar(final1_5$CC_Sample_Size, final1_5$PATT, final1_5$UPPER_CF, final1_5$LOWER_CF,
#                xlab = "CC Size", ylab = "PATT Estimation")


#################### GG Plot Error Bar #############
#ggplot(data=final, aes(x=CC_Sample_Size, y=PATT)) + geom_line(aes(colour=Ratio))

p1 <- ggplot(final1, aes(x=CC_Sample_Size, y=PATT, cex.axis=0.75)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0,TA_size,50)) +
  geom_errorbar(aes(ymin=LOWER_CF, ymax=UPPER_CF), width=.2,
                position=position_dodge(0.05)) + labs(caption = "Ratio 1:1") +
  geom_hline(yintercept=-15.174, linetype="dashed", color = "red") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p2 <- ggplot(final2, aes(x=CC_Sample_Size, y=PATT)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0,TA_size,50)) +
  geom_errorbar(aes(ymin=LOWER_CF, ymax=UPPER_CF), width=.2,
                position=position_dodge(0.05)) + labs(caption = "Ratio 1:2") +
  geom_hline(yintercept=-15.174, linetype="dashed", color = "red") +
  theme_classic() + theme(axis.text.x = element_text(angle = 90))
p1_5 <- ggplot(final1_5, aes(x=CC_Sample_Size, y=PATT)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq(0,TA_size,50)) +
  geom_errorbar(aes(ymin=LOWER_CF, ymax=UPPER_CF), width=.2,
                position=position_dodge(0.05)) + labs(caption = "Ratio 1:1.5") +
  geom_hline(yintercept=-15.174, linetype="dashed", color = "red")+
  theme_classic() + theme(axis.text.x = element_text(angle = 90))

p1+
  labs(title = "CC Sample Size VS PATT Estimate (TA=1500)",
       subtitle = "Varying CC Size from 0 to 1500 with 50 stepsize")+
  p1_5+p2+plot_layout(ncol=1)
####################################################

