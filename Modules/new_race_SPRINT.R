setwd("/Users/nafisneehal")
seatsys_data <- read.csv("./18M_seatsys.csv")
incl_excl_data <- read.csv("./incl_excl.csv")

seatsys_data <- na.omit(seatsys_data)
incl_excl_data <- incl_excl_data %>% replace(is.na(.), 0)

race_dat <- incl_excl_data %>% select(MASKID, SPANISHNO: RACE_OTHER)

#New RACE hierarchy
#Others -> 0, Non-Hispanic Asian -> 1, Non-Hispanic Black -> 2, Non-Hispanic White -> 3, Hispanic -> 4
# race_dat <- race_dat %>% mutate(RACE = ifelse(SPANISHNO==0, 4, ifelse(SPANISHNO==1 & RACE_BLACK==1, 2,
#                                                                                ifelse(SPANISHNO==1 & RACE_ASIAN==1, 1,
#                                                                                       ifelse(SPANISHNO==1 & RACE_WHITE==1, 3, 0)))))
race_dat <- race_dat %>% mutate(# Create race and Hispanic ethnicity categories for  analysis
  RACE = ifelse(is.na(SPANISHNO),'Hispanic',
                             ifelse(SPANISHNO==1 & is.na(RACE_WHITE)& is.na(RACE_BLACK)&is.na(RACE_ASIAN),'Other',
                                    ifelse(is.na(SPANISHNO) & is.na(RACE_WHITE)& is.na(RACE_BLACK)&is.na(RACE_ASIAN),'Other',
                                           ifelse(SPANISHNO==1 & RACE_WHITE==1 & is.na(RACE_BLACK) & is.na(RACE_ASIAN),'NH White',
                                                  ifelse(SPANISHNO==1 & RACE_BLACK==1 & is.na(RACE_WHITE) & is.na(RACE_ASIAN),'NH Black',
                                                         ifelse(SPANISHNO==1 & RACE_ASIAN==1 & is.na(RACE_BLACK) & is.na(RACE_WHITE),'NH Asian',
                                                                ifelse(RACE_WHITE==1 & RACE_BLACK==1,'Other',
                                                                       ifelse(RACE_BLACK==1 & RACE_ASIAN==1,'Other',
                                                                              ifelse(RACE_WHITE==1 & RACE_ASIAN==1,'Other',
                                                                                     ifelse(RACE_WHITE==1 & RACE_ASIAN==1 & RACE_BLACK==1,'Other',NA)))))))))))
table(race_dat %>% select(RACE))
merged <- merge(seatsys_data, race_dat %>% select(MASKID, RACE), by='MASKID')

write.csv(merged, "./18M_seatsys2.csv", row.names = FALSE)
