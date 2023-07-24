#### SOURCES ####
setwd("/home/neehan/data/Nafis/ESCA_Primary_Git")
source("./Wrapper/load_packages.R")
source("./Wrapper/wrapper_modules.R")
#source("./Wrapper/load_data_formatted.R")

allhat_data<-read.sas7bdat("./Data/ALLHAT_data/allhat_key.sas7bdat")

allhat_selected_raw<- allhat_data[c('RZGROUP','AGE','SEX','RACE','HISPANIC','EDUCAT', 'BLMEDS2',
                                    'BV2SBP', 'BV2DBP', 'CURSMOKE', 'MISTROKE', 'HXCABG', 
                                    'OASCVD', 'STDEPR', 'DIABETES', 'HDLLT35', 'LVHECG', 'WALL25',
                                    'LCHD', 'BLBMI', 'AFGLUC','ACHOL','EP_CHD', 'DYCHD')]

#1 means Yes, 2 means No --> in all binary cases (unless indicated otherwise)
allhat_selected_processed<- mutate(allhat_selected_raw,
                                   RANDASSIGN = RZGROUP,
                                   # Create age categories for adults aged 18 and over
                                   Age_Group = ifelse(AGE<18,'<18',ifelse(AGE<40,'18-39',ifelse(AGE<60,'40-59','59+'))),
                                   Gender = factor(c(1,2)[SEX],
                                                   labels = c('Male', 'Female')),      
                                   # Create race and Hispanic ethnicity categories for  analysis 
                                   Race_or_Ethnicity = ifelse(HISPANIC==1,'Hispanic',ifelse(HISPANIC==2 & RACE==1,'NH White',ifelse(HISPANIC==2 & RACE==2,'NH Black',ifelse(HISPANIC==2 & RACE==4,'NH Asian','Other')))),
                                   Education = ifelse(EDUCAT< 1,NA,
                                                      ifelse(EDUCAT<11,'<HSG',
                                                             ifelse(EDUCAT==12,'HSG/GED',ifelse(EDUCAT<17,'Some college/TS', 
                                                                                                ifelse(EDUCAT<26,'>=College grad',NA))))),
                                   Prior_HypTreat = ifelse(BLMEDS2==1, 'Yes', ifelse(BLMEDS2==2, 'No', NA)),
                                   SBP = BV2SBP, 
                                   DBP = BV2DBP,

                                   #Risk Categories
                                   Smoker =ifelse(CURSMOKE==1,'Smoke',ifelse(CURSMOKE==2,'No smoke',ifelse(CURSMOKE==3,'No smoke',NA))),
                                   Had_MIS = ifelse(MISTROKE==1, 'Yes', ifelse(MISTROKE==2,'NO',NA)),
                                   Had_CRV = ifelse(HXCABG==1, 'Yes', ifelse(HXCABG==2,'NO',NA)),
                                   Had_ASCVD = ifelse(OASCVD==1, 'Yes', ifelse(OASCVD==2,'NO',NA)),
                                   Had_MajorST = ifelse(STDEPR==1, 'Yes', ifelse(STDEPR==2,'NO',NA)),
                                   Had_T2D = ifelse(DIABETES==1, 'Yes', ifelse(DIABETES==2,'NO',NA)),
                                   Had_HDLC = ifelse(HDLLT35==1, 'Yes', ifelse(HDLLT35==2,'NO',NA)),
                                   Had_LVH_ELCT = ifelse(LVHECG==1, 'Yes', ifelse(LVHECG==2,'NO',NA)),
                                   Had_LVH_ECHO = ifelse(WALL25==1, 'Yes', ifelse(WALL25==2,'NO',NA)),
                                   
                                   #others
                                   Has_CHD = ifelse(LCHD==1, 'Yes', ifelse(LCHD==2,'NO',NA)),
                                   BMI = BLBMI,
                                   
                                   #FPG = ifelse(AFGLUC <100,'Glucose<100',ifelse(AFGLUC <126,'Glucose 100-125','Glucose>=126')),
                                   #TC = ifelse(ACHOL<200,'Normal TC','High TC'),
                                   EVENT = EP_CHD,
                                   EVENTDAYS = DYCHD
                                   
)%>% select(RANDASSIGN:EVENTDAYS)


allhat_selected_processed$Race_or_Ethnicity <- factor(allhat_selected_processed$Race_or_Ethnicity, levels = c("NH White", "NH Black", "NH Asian", "Hispanic", "Other"))
allhat_selected_processed$Education<- factor(allhat_selected_processed$Education, levels = c("<HSG","HSG/GED", "Some college/TS", ">=College grad"))
#allhat_selected_processed$FPG <- factor(allhat_selected_processed$FPG, levels = c("Glucose<100", "Glucose 100-125", "Glucose>=126"))

allhat_final <- na.omit(allhat_selected_processed)

write.csv(allhat_selected_processed,'./Data/ALLHAT_data/ALLHAT_processed_with_NA.csv') 
write.csv(allhat_final,'./Data/ALLHAT_data/ALLHAT_processed.csv') 



