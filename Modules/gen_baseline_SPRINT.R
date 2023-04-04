suppressPackageStartupMessages({
  library(sas7bdat)
  library(dplyr)
  library(haven)
  library(survey)
  library(tidyr)
  library(plotly)
  library(ggplot2)
})

setwd("/Users/nafisneehal/ESCA")
source("modules/equity_metrics_miao.R")

### Helper Functions #####

bmi_calculation_lbs_inch<-function(weight, height){
  return(703 * weight / (height)^2)
}

bmi_calculation_kg_m<-function(weight, height_cm){
  height_m<-0.01*height_cm
  return(weight / (height_m)^2)
}

######
sprint_files <- list.files(path = "./Data/SPRINT_data", pattern = "*.csv", full.names = T)
sprint_data_repeated<- sapply(sprint_files, read.csv, simplify=FALSE)

sprint_data<- sprint_data_repeated[[3]] %>% full_join(sprint_data_repeated[[5]],by= "MASKID")%>% 
  select("MASKID", 'NEWSITEID', 'RANDASSIGN', 'CVDHISTORY', 'GENDER', 'SPANISHNO', 
         'RACE_WHITE', 'RACE_BLACK', 'RACE_ASIAN','SBPAVG','TOTALCHOLEST', 'SERUMCREAT', 'GFRESTIMATE', 
         'CIGSMOKER','RZ_AGE', 'CVDPOINTS')

sprint_labs<-sprint_data_repeated[[4]]%>%select("MASKID",'VISITCODE','FASTING_GLUR','RESULT_GLUR')
sprint_labs<-sprint_labs[which(sprint_labs$VISITCODE=='RZ1' ),] #on baseline
sprint_labs$VISITCODE<-NULL
sprint_labs<-sprint_labs[which(!is.na(sprint_labs$FASTING_GLUR)&!is.na(sprint_labs$RESULT_GLUR) ),]
sprint_edu<-sprint_data_repeated[[1]]%>%select("MASKID",'EDUCATION')
sprint_bmi<-sprint_data_repeated[[2]]%>%select("MASKID",'WEIGHT','HEIGHT')
sprint_bmi$bmi <- bmi_calculation_lbs_inch(sprint_bmi$WEIGHT,sprint_bmi$HEIGHT)

sprint_data<- sprint_data %>% full_join(sprint_edu,by= "MASKID") %>% full_join(sprint_bmi,by= "MASKID")%>% full_join(sprint_labs,by= "MASKID")

sprint_selected_processed<- mutate(sprint_data,
                                   # Create age categories for adults aged 18 and over
                                   Age_Group = ifelse(RZ_AGE<18,'<18',ifelse(RZ_AGE<40,'18-39',ifelse(RZ_AGE<60,'40-59','59+'))),
                                   Gender = factor(c(1,2)[GENDER],
                                                   labels = c('Female','Male')),
                                   # Create race and Hispanic ethnicity categories for  analysis
                                   Race_or_Ethnicity = ifelse(is.na(SPANISHNO),'Hispanic',
                                                              ifelse(SPANISHNO==1 & is.na(RACE_WHITE)& is.na(RACE_BLACK)&is.na(RACE_ASIAN),'Other',
                                                                     ifelse(is.na(SPANISHNO) & is.na(RACE_WHITE)& is.na(RACE_BLACK)&is.na(RACE_ASIAN),'Other',
                                                                            ifelse(SPANISHNO==1 & RACE_WHITE==1 & is.na(RACE_BLACK) & is.na(RACE_ASIAN),'NH White',
                                                                                   ifelse(SPANISHNO==1 & RACE_BLACK==1 & is.na(RACE_WHITE) & is.na(RACE_ASIAN),'NH Black',
                                                                                          ifelse(SPANISHNO==1 & RACE_ASIAN==1 & is.na(RACE_BLACK) & is.na(RACE_WHITE),'NH Asian',
                                                                                                 ifelse(RACE_WHITE==1 & RACE_BLACK==1,'Other',
                                                                                                        ifelse(RACE_BLACK==1 & RACE_ASIAN==1,'Other',
                                                                                                               ifelse(RACE_WHITE==1 & RACE_ASIAN==1,'Other',
                                                                                                                      ifelse(RACE_WHITE==1 & RACE_ASIAN==1 & RACE_BLACK==1,'Other',NA)))))))))),
                                   Education =
                                     ifelse(EDUCATION<4,'<HSG',
                                            ifelse(EDUCATION==5,'HSG/GED',
                                                   ifelse(EDUCATION<9,'Some college/TS',
                                                          ifelse(EDUCATION>8,'>=College grad',NA)))),
                                   Smoker =ifelse(CIGSMOKER==1,'Smoke','No smoke'),
                                   # Create bmi categories for analysis
                                   BMI = ifelse(bmi<18.5,'Underweight',ifelse(bmi<25.0,'Normal weight',ifelse(bmi<30.0,'Overweight','Obese'))),
                                   SBP = ifelse(SBPAVG<120,'SBP<120',ifelse(SBPAVG<130,'SBP 120-129',ifelse(SBPAVG<140,'SBP 130-139','SBP>=140'))),
                                   FPG = ifelse(FASTING_GLUR == "YES" & RESULT_GLUR <100,'Glucose<100',ifelse(FASTING_GLUR == "YES" & RESULT_GLUR <126,'Glucose 100-125',ifelse(FASTING_GLUR == "YES" & RESULT_GLUR >=126, 'Glucose>=126', NA))),
                                   TC = ifelse(TOTALCHOLEST<200,'Normal TC','High TC')
) %>% select(MASKID, NEWSITEID, RANDASSIGN, Age_Group:Smoker, SBP:TC, CVDHISTORY, 
             CVDPOINTS, SERUMCREAT, GFRESTIMATE)


sprint_selected_processed$Race_or_Ethnicity <- factor(sprint_selected_processed$Race_or_Ethnicity, levels = c("NH White", "NH Black", "NH Asian", "Hispanic", "Other"))
sprint_selected_processed$Education <- factor(sprint_selected_processed$Education, levels = c("<HSG", "HSG/GED", "Some college/TS", ">=College grad"))
sprint_selected_processed$SBP <- factor(sprint_selected_processed$SBP, levels = c("SBP<120", "SBP 120-129", "SBP 130-139", "SBP>=140"))
sprint_selected_processed$FPG <- factor(sprint_selected_processed$FPG, levels = c("Glucose<100", "Glucose 100-125", "Glucose>=126"))

### BMI ###
#sprint_selected_processed$BMI <- factor(sprint_selected_processed$BMI, levels = c("Underweight", "Normal weight", "Overweight", "Obese"))
# sprint_all <- count(sprint_selected_processed,Gender,Age_Group,Race_or_Ethnicity,
#                     Education,SBP,BMI,Smoker,TC,FPG) %>%
#   rename(user_n=`n`)
### BMI ###

sprint_all <- count(sprint_selected_processed, Gender,Age_Group,Race_or_Ethnicity) %>% 
  rename(user_n=`n`)

write.csv(sprint_all,"./Data/Processed/SPRINT_count.csv") # contains counts in each subgroup
write.csv(sprint_selected_processed,'./Data/Processed/SPRINT_example.csv') #contains the factorized data for each sprint sample 9361
write.csv(sprint_data,'./Data/Processed/SPRINT_numeric.csv') #contains numeric representations


##### SPRINT Monthly BP Data Load - Treated
sprint_data_1m <- read.csv("./Data/SPRINT_Original_Full_Miao/SPRINT/data/CSV/intbp_manage_1m.csv")
sprint_data_2m <- read.csv("./Data/SPRINT_Original_Full_Miao/SPRINT/data/CSV/intbp_manage_2m.csv")
sprint_data_6m <- read.csv("./Data/SPRINT_Original_Full_Miao/SPRINT/data/CSV/intbp_manage_6m.csv")
sprint_data_18m <- read.csv("./Data/SPRINT_Original_Full_Miao/SPRINT/data/CSV/intbp_manage_18m.csv")
sprint_data_1m <- sprint_data_1m %>% select(MASKID, FORMDAYS, VISITCODE, SEATSYS)
sprint_data_2m <- sprint_data_2m %>% select(MASKID, FORMDAYS, VISITCODE, SEATSYS)
sprint_data_6m <- sprint_data_6m %>% select(MASKID, FORMDAYS, VISITCODE, SEATSYS)
sprint_data_18m <- sprint_data_18m %>% select(MASKID, FORMDAYS, VISITCODE, SEATSYS)

##### SPRINT Monthly BP Data Load - Controls
stbp_anfu <- read.csv("./Data/SPRINT_Original_Full_Miao/SPRINT/data/CSV/stbp_manage_anfu.csv")
stbp_nanfu <- read.csv("./Data/SPRINT_Original_Full_Miao/SPRINT/data/CSV/stbp_manage_nanfu.csv")
stbp_anfu <- stbp_anfu %>% select(MASKID, FORMDAYS, VISITCODE, SEATSYS)
stbp_nanfu <- stbp_nanfu %>% select(MASKID, FORMDAYS, VISITCODE, SEATSYS)

#merge treated and controls
sprint_data_monthly_bp <- do.call("rbind", list(sprint_data_1m, sprint_data_2m, sprint_data_6m,
                                                sprint_data_18m, stbp_anfu, stbp_nanfu))

#write out the data file
write.csv(sprint_data_monthly_bp, "./Data/Processed/SPRINT_monthly_bp.csv") #contains monthly BP examination data
#####