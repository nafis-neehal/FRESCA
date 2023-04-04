suppressPackageStartupMessages({
  library(MatchIt)
  library(table1)
})

# dat <- read.csv("./Data/XTY/18M_seatsys.csv")
# dat <- na.omit(dat)
# 
# dat$DIV <- 0
# dat[dat$RANDASSIGN==1, ]$DIV <- 1
# 
# DIV_MAP <- data.frame("groupname"=c("Treated","Trial Controls"),
#                       "divvalue"=c(1,0))

p_val_threshold <- 0.05

all_p_controls <- list()

# ####### P-Val Function ########
pvalue <- function(x, ...) {
  
  # Construct vectors of data y, and groups (strata) g
  y <- unlist(x)
  #print(sapply(x, length))
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  
  if (is.numeric(y)) {
    # For numeric variables, perform a standard 2-sample t-test
    p <- t.test(y ~ g)$p.value
  } else {
    # For categorical variables, perform a chi-squared test of independence
    p <- chisq.test(table(y, g))$p.value
  }
  
  # Format the p-value, using an HTML entity for the less-than sign.
  # The initial empty string places the output on the line below the variable label.
  c("", sub("<", "<", format.pval(p, digits=3, eps=p_val_threshold)))
  
}

##### Gen P-Val Function ########
gen_p_val_controls <- function(data, DIV_MAP){

  treated_dat <- data[data$DIV==1, ]
  
  #LOOP through each control group
  for (controls in DIV_MAP$groupname[2:nrow(DIV_MAP)]){
    
    #get the DIV value for that control group
    DIV_val <- DIV_MAP[DIV_MAP$groupname==controls,]$divvalue
    
    #get the slice of data for that control div
    control_dat <- data[data$DIV==DIV_val,]
    
    #merge both groups
    dat <- rbind(treated_dat, control_dat)
    
    #set the variable categories (Don't Bother about the Labels now)
    dat$DIV <- factor(dat$DIV, levels=c(1, DIV_val), 
                      labels = c("Treated", controls))
    dat$SUB_SENIOR <- as.logical((dat$SUB_SENIOR))
    dat$SUB_CKD <- as.logical((dat$SUB_CKD))
    dat$RZ_AGE <- factor(dat$RZ_AGE, levels = 0:3, labels=c("<18", '18-39', '40-59', '59+'))
    dat$CVDHISTORY <- as.logical((dat$CVDHISTORY))
    dat$GENDER <- factor(dat$GENDER, levels=1:2, labels=c("Female", "Male"))
    dat$RACE_ <- factor(dat$RACE_, levels=0:3, labels=c("Asian", "Black", "White", "Others"))
    dat$CIGSMOKER <- as.logical(dat$CIGSMOKER)
    
    
    #plot table
    x<-table1(~ SUB_SENIOR + SUB_CKD +
                CVDHISTORY + RZ_AGE + GENDER + RACE_ + SERUMCREAT + GFRESTIMATE
              + CIGSMOKER + CVDPOINTS + SEATSYS | DIV, data = dat, overall=F,
              extra.col=list(`PVal`=pvalue))

    x<- as.data.frame(x)

    #conditions on p_vals
    x[x$PVal=='<0.05',][[controls]] <- paste(x[x$PVal=='<0.05',][[controls]], "*", sep = "")

    all_p_controls[[controls]] <- x[[controls]]
    
  }
  
  return (all_p_controls)
}

#all_controls <- gen_p_val_controls(dat, DIV_MAP)
