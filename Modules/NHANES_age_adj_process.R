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

# Import the datasets using read_xpt()
nhanes_demo <- read_xpt("./Miao Stuff/preprocess/background_data/demographics_nhanes_15_16.XPT")
nhanes_diabetes <- read_xpt("./Miao Stuff/preprocess/background_data/diabetes_nhanes_15_16.XPT")
nhanes_hyp <- read_xpt("./Miao Stuff/preprocess/background_data/antihypertension_nhanes_15_16.XPT")
nhanes_BP <- read_xpt("./Miao Stuff/preprocess/background_data/blood_pressure_nhanes_15_16.XPT")
nhanes_BM <- read_xpt("./Miao Stuff/preprocess/background_data/body_measure_nhanes_15_16.XPT")
nhanes_smoker <- read_xpt("./Miao Stuff/preprocess/background_data/smoking_nhanes_15_16.XPT")
nhanes_TC <- read_xpt("./Miao Stuff/preprocess/background_data/total_cholesterol_nhanes_15_16.XPT")
nhanes_biochem <- read_xpt("./Miao Stuff/preprocess/background_data/biochem_nhanes_15_16.XPT")
nhanes_hba1c <- read_xpt("./Miao Stuff/preprocess/background_data/hba1c_nhanes_15_16.XPT")
nhanes_FPG <- read_xpt("./Miao Stuff/preprocess/background_data/FPG_nhanes_15_16.XPT")


# Merge the datasets to create nhanes_patients
nhanes_patients <- list(nhanes_demo, nhanes_diabetes,nhanes_hyp,nhanes_BP,nhanes_BM,nhanes_smoker,nhanes_TC,nhanes_biochem,nhanes_hba1c) %>%
  Reduce(function(df1, df2) full_join(df1, df2, by = "SEQN"), .)
# Make sure the dataset is read in as an R data frame
df_nhanes_patients<- as.data.frame(nhanes_patients)
write.csv(df_nhanes_patients,'./Data/NHANES/NHANES_merged/NHANES_background.csv')

# Merge the datasets to create nhanes_patients_lab
nhanes_patients_lab <- list(nhanes_demo, nhanes_diabetes,nhanes_hyp,nhanes_TC,nhanes_hba1c,nhanes_FPG) %>%
  Reduce(function(df1, df2) full_join(df1, df2, by = "SEQN"), .)
# Make sure the dataset is read in as an R data frame
df_nhanes_patients_lab<- as.data.frame(nhanes_patients_lab)
write.csv(df_nhanes_patients_lab,'./Data/NHANES/NHANES_merged/NHANES_background_lab.csv')

nhanes_patients_combined <- list(nhanes_demo, nhanes_diabetes,nhanes_hyp,nhanes_BP,nhanes_BM,nhanes_smoker,nhanes_TC,nhanes_hba1c,nhanes_FPG,nhanes_biochem) %>%
  Reduce(function(df1, df2) full_join(df1, df2, by = "SEQN"), .)
# Make sure the dataset is read in as an R data frame
df_nhanes_patients_combined<- as.data.frame(nhanes_patients_combined)
write.csv(df_nhanes_patients_combined,'./Data/NHANES/NHANES_merged/NHANES_background_combined.csv')


####################### NHANES - Data Count ###################

#upload background datasets
df_nhanes_patients<-read.csv("./Data/NHANES/NHANES_merged/NHANES_background.csv")
df_nhanes_patients$X<-NULL
df_nhanes_patients_lab<-read.csv("./Data/NHANES/NHANES_merged/NHANES_background_lab.csv")
df_nhanes_patients_lab$X<-NULL
df_nhanes_patients_combined<-read.csv("./Data/NHANES/NHANES_merged/NHANES_background_combined.csv")
df_nhanes_patients_combined$X<-NULL

# Background: Hypertension for Demographic characteristics & Risk Factors
df_background_hypertension <- mutate(df_nhanes_patients,
                                     # Create age categories for adults aged 18 and over
                                     Age_Group = ifelse(RIDAGEYR<18,'under 18',ifelse(RIDAGEYR<40,'18-39',ifelse(RIDAGEYR<60,'40-59',"59+"))),
                                     Gender = factor(c(1,2)[RIAGENDR],
                                                     labels = c("Male", "Female")),      
                                     # Create race and Hispanic ethnicity categories for  analysis 
                                     Race_or_Ethnicity = ifelse (RIDRETH3 <3,'Hispanic',ifelse (RIDRETH3 == 3,'NH White', ifelse (RIDRETH3 == 4,'NH Black', ifelse (RIDRETH3 == 6,'NH Asian', ifelse (RIDRETH3 == 7,'Other', NA))))),
                                     # Create age categories for adults aged 18 and over
                                     
                                     Education = factor(c(1,2,3,4,5)[DMDEDUC2],
                                                        labels = c('<HSG','<HSG','HSG/GED','Some college/TS', '>=College grad')), 
                                     Smoker =factor(c(1,2,3)[SMQ040],
                                                    labels = c("Smoke", "Smoke", "No smoke")),
                                     # Create bmi categories for analysis 
                                     BMI = ifelse(BMXBMI <18.5,'Underweight',ifelse(BMXBMI <25.0,'Normal weight',ifelse(BMXBMI <30.0,'Overweight',ifelse(BMXBMI >=30.0,'Obese',NA)))),
                                     SBP = ifelse(BPXSY1 <120,'SBP<120',ifelse(BPXSY1 <130,'SBP 120-129',ifelse(BPXSY1 <140,'SBP 130-139',ifelse(BPXSY1>=140,'SBP>=140', NA)))),
                                     hypertentive = case_when(BPQ020==1 ~ 1, BPQ020!=1~0),
                                     #  Define subpopulation of interest:  people aged 18
                                     inAnalysis= (RIDAGEYR >=50 & (hypertentive == 1))
)

df_background_hypertension$Race_or_Ethnicity <- factor(df_background_hypertension$Race_or_Ethnicity, levels = c("NH White", "NH Black","NH Asian", "Hispanic", "Other"))
df_background_hypertension$Education <- factor(df_background_hypertension$Education, levels = c("<HSG", "HSG/GED", "Some college/TS", ">=College grad"))
df_background_hypertension$BMI <- factor(df_background_hypertension$BMI, levels = c("Underweight", "Normal weight", "Overweight", "Obese"))
df_background_hypertension$SBP <- factor(df_background_hypertension$SBP, levels = c("SBP<120", "SBP 120-129", "SBP 130-139", "SBP>=140"))

NHANES_all_hypertension <- svydesign(data=df_background_hypertension, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)
NHANES_hypertension <- subset(NHANES_all_hypertension, inAnalysis==1)



# Background: Hypertension for Lab Results
df_background_hypertension_lab <- mutate(df_nhanes_patients_lab,
                                         TC = ifelse(LBXTC <200,'Normal TC','High TC'),
                                         HbA1c = ifelse(LBXGH<7,'Normal',ifelse(LBXGH <9,'Elevated','Severaly elevated')) ,
                                         FPG = ifelse(LBXGLU <100,'Glucose<100',ifelse(LBXGLU <126,'Glucose 100-125','Glucose>=126')),
                                         hypertentive = case_when(BPQ020==1 ~ 1, BPQ020!=1~0),
                                         # Define subpopulation of interest:  people aged 18
                                         inAnalysis= (RIDAGEYR >=50 & (hypertentive == 1))
)
df_background_hypertension_lab<-df_background_hypertension_lab[which(!is.na(df_background_hypertension_lab$WTSAF2YR )),]

df_background_hypertension_lab$FPG <- factor(df_background_hypertension_lab$FPG, levels = c("Glucose<100", "Glucose 100-125", "Glucose>=126"))

NHANES_all_hypertension_lab <- svydesign(data=df_background_hypertension_lab, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTSAF2YR, nest=TRUE)
NHANES_hypertension_lab <- subset(NHANES_all_hypertension_lab, inAnalysis==1)



# Background: Hypertension for combined
df_background_hypertension_combined <- mutate(df_nhanes_patients_combined,
                                              # Create age categories for adults aged 18 and over
                                              Age_Group = ifelse(RIDAGEYR<18,'under 18',ifelse(RIDAGEYR<45,'18-44',ifelse(RIDAGEYR<65,'45-64',"64+"))),
                                              Gender = factor(c(1,2)[RIAGENDR],
                                                              labels = c("Male", "Female")),      
                                              # Create race and Hispanic ethnicity categories for  analysis 
                                              Race_or_Ethnicity = ifelse (RIDRETH3 <3,'Hispanic',ifelse (RIDRETH3 == 3,'NH White', ifelse (RIDRETH3 == 4,'NH Black', ifelse (RIDRETH3 == 6,'NH Asian', ifelse (RIDRETH3 == 7,'Other', NA))))),
                                              
                                              # Create age categories for adults aged 18 and over
                                              
                                              Education = factor(c(1,2,3,4,5)[DMDEDUC2],
                                                                 labels = c('<HSG','<HSG','HSG/GED','Some college/TS', '>=College grad')), 
                                              Smoker =factor(c(1,2,3)[SMQ040],
                                                             labels = c("Smoke", "Smoke", "No smoke")),
                                              # Create bmi categories for analysis 
                                              BMI = ifelse(BMXBMI <18.5,'Underweight',ifelse(BMXBMI <25.0,'Normal weight',ifelse(BMXBMI <30.0,'Overweight',ifelse(BMXBMI >=30.0,'Obese',NA)))),
                                              SBP = ifelse(BPXSY1 <120,'SBP<120',ifelse(BPXSY1 <130,'SBP 120-129',ifelse(BPXSY1 <140,'SBP 130-139',ifelse(BPXSY1>=140,'SBP>=140', NA)))),
                                              TC = ifelse(LBXTC <200,'Normal TC','High TC'),
                                              HbA1c = ifelse(LBXGH<7,'Normal',ifelse(LBXGH <9,'Elevated','Severaly elevated')) ,
                                              FPG = ifelse(LBXGLU <100,'Glucose<100',ifelse(LBXGLU <126,'Glucose 100-125','Glucose>=126')),
                                              hypertentive = case_when(BPQ020==1 ~ 1, BPQ020!=1~0),                                 
                                              # Define subpopulation of interest:  people aged 18
                                              inAnalysis= (RIDAGEYR >=50 & (hypertentive == 1))
)

df_background_hypertension_combined<-df_background_hypertension_combined[which(!is.na(df_background_hypertension_combined$WTSAF2YR )),]

df_background_hypertension_combined$Race_or_Ethnicity <- factor(df_background_hypertension_combined$Race_or_Ethnicity, levels = c("NH White", "NH Black","NH Asian", "Hispanic", "Other"))
df_background_hypertension_combined$Education <- factor(df_background_hypertension_combined$Education, levels = c("<HSG", "HSG/GED", "Some college/TS", ">=College grad"))
df_background_hypertension_combined$BMI <- factor(df_background_hypertension_combined$BMI, levels = c("Underweight", "Normal weight", "Overweight", "Obese"))
df_background_hypertension_combined$SBP <- factor(df_background_hypertension_combined$SBP, levels = c("SBP<120", "SBP 120-129", "SBP 130-139", "SBP>=140"))
df_background_hypertension_combined$FPG <- factor(df_background_hypertension_combined$FPG, levels = c("Glucose<100", "Glucose 100-125", "Glucose>=126"))

NHANES_all_hypertension_combined <- svydesign(data=df_background_hypertension_combined, id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTSAF2YR, nest=TRUE)
NHANES_hypertension_combined <- subset(NHANES_all_hypertension_combined, inAnalysis==1)


### Final DF
nhanes_hypertension<-as.data.frame(svytable(~Gender+Age_Group+Race_or_Ethnicity+Education+SBP+BMI+Smoker,NHANES_hypertension,addNA = TRUE,na.action=NULL,round=TRUE))%>%rename(background_n=`Freq`)
nhanes_hypertension_lab<-as.data.frame(svytable(~TC+FPG,NHANES_hypertension_lab,addNA = TRUE,na.action=NULL,round=TRUE))%>%rename(background_n=`Freq`)
nhanes_hypertension_combined<-as.data.frame(svytable(~Gender+Age_Group+Race_or_Ethnicity+Education+SBP+BMI+Smoker+TC+FPG,NHANES_hypertension_combined,addNA = TRUE,na.action=NULL,round=TRUE))%>%rename(background_n=`Freq`)

write.csv(nhanes_hypertension,"./Data/NHANES/NHANES_age_adj/nhanes_hypertension.csv")
write.csv(nhanes_hypertension_lab,"./Data/NHANES/NHANES_age_adj/nhanes_hypertension_lab.csv")
write.csv(nhanes_hypertension_combined,"./Data/NHANES/NHANES_age_adj/nhanes_hypertension_combined.csv")




