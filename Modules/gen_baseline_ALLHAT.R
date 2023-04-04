setwd("/Users/nafisneehal/ESCA")

suppressPackageStartupMessages({
  library(sas7bdat)
  library(dplyr)
  library(haven)
  library(survey)
  library(tidyr)
  library(plotly)
  library(ggplot2)
})

allhat_data<-read.sas7bdat("./Data/ALLHAT_data/allhat_key.sas7bdat")

allhat_selected_raw<- allhat_data[c('AGE','SEX','RACE','HISPANIC','EDUCAT','CURSMOKE','BLBMI','BV2SBP','AFGLUC','ACHOL')]

allhat_selected_processed<- mutate(allhat_selected_raw,
                                   # Create age categories for adults aged 18 and over
                                   Age_Group = ifelse(AGE<18,'<18',ifelse(AGE<40,'18-39',ifelse(AGE<60,'40-59','59+'))),
                                   Gender = factor(c(1,2)[SEX],
                                                   labels = c('Male', 'Female')),      
                                   # Create race and Hispanic ethnicity categories for  analysis 
                                   Race_or_Ethnicity = ifelse(HISPANIC==1,'Hispanic',ifelse(HISPANIC==2 & RACE==1,'NH White',ifelse(HISPANIC==2 & RACE==2,'NH Black',ifelse(HISPANIC==2 & RACE==4,'NH Asian','Other')))),
                                   Education = ifelse(EDUCAT< 1,NA,
                                                      ifelse(EDUCAT<11,'<HSG',
                                                             ifelse(EDUCAT==12,'HSG/GED',ifelse(EDUCAT<17,'Some college/TS', ifelse(EDUCAT<26,'>=College grad',NA))))),
                                   Smoker =ifelse(CURSMOKE==1,'Smoke',ifelse(CURSMOKE==2,'No smoke',ifelse(CURSMOKE==3,'No smoke',NA))),
                                   # Create bmi categories for analysis 
                                   BMI = ifelse(BLBMI<18.5,'Underweight',ifelse(BLBMI<25.0,'Normal weight',ifelse(BLBMI<30.0,'Overweight','Obese'))),
                                   SBP = ifelse(BV2SBP<120,'SBP<120',ifelse(BV2SBP<130,'SBP 120-129',ifelse(BV2SBP<140,'SBP 130-139','SBP>=140'))),
                                   
                                   FPG = ifelse(AFGLUC <100,'Glucose<100',ifelse(AFGLUC <126,'Glucose 100-125','Glucose>=126')),
                                   TC = ifelse(ACHOL<200,'Normal TC','High TC'),
                                   
)%>% select(Age_Group:TC)


allhat_selected_processed$Race_or_Ethnicity <- factor(allhat_selected_processed$Race_or_Ethnicity, levels = c("NH White", "NH Black", "NH Asian", "Hispanic", "Other"))
allhat_selected_processed$Education<- factor(allhat_selected_processed$Education, levels = c("<HSG","HSG/GED", "Some college/TS", ">=College grad"))
allhat_selected_processed$BMI <- factor(allhat_selected_processed$BMI, levels = c("Underweight", "Normal weight", "Overweight", "Obese"))
allhat_selected_processed$SBP <- factor(allhat_selected_processed$SBP, levels = c("SBP<120", "SBP 120-129", "SBP 130-139", "SBP>=140"))
allhat_selected_processed$FPG <- factor(allhat_selected_processed$FPG, levels = c("Glucose<100", "Glucose 100-125", "Glucose>=126"))

write.csv(allhat_selected_processed,'ALLHAT_example.csv')