# library(MatchIt)
# library(rdist)
# 
# data <- read.csv("./Data/XTY/18M_seatsys.csv")
# data <- na.omit(data)
# 
# input_threshold <- 100
# input_control_fraction <- 0.7

#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

matching_policy_0 <- function(trial_dat, sca_dat, num_frac_th_controls, num_frac_sca_controls){
  
  #Trial Controls Distance Matrix
  trial_matched.data <- matchit(RANDASSIGN ~ SUB_SENIOR + SUB_CKD + 
                                  CVDHISTORY + RZ_AGE + GENDER + RACE_ + SERUMCREAT + GFRESTIMATE
                                + CIGSMOKER + CVDPOINTS, data = trial_dat, #<-
                                method = 'nearest', distance='glm', 
                                replace=TRUE)
  
  trial_matches <- get_matches(trial_matched.data)
  
  tr_psm <- trial_matches[trial_matches$RANDASSIGN==1,][,7:16]
  cn_psm <- trial_matches[trial_matches$RANDASSIGN==0,][,7:16]
  
  #from each control to all the treated distance in each row
  trial_dist_mat <- cdist(cn_psm, tr_psm, metric="euclidean")
  trial_dist_mat_sum <- rowSums(trial_dist_mat)
  trial_dist_df <- as.data.frame(trial_dist_mat_sum)
  trial_dist_df$cn_ids <- trial_matches[trial_matches$RANDASSIGN==0,]$id
  trial_dist_df$cn_maskIDs <- trial_matches[trial_matches$RANDASSIGN==0,]$MASKID
  
  #SCA Controls Distance Matrix
  sca_matched.data <- matchit(RANDASSIGN ~ SUB_SENIOR + SUB_CKD +
                                CVDHISTORY + RZ_AGE + GENDER + RACE_ + SERUMCREAT + GFRESTIMATE
                              + CIGSMOKER + CVDPOINTS, data = sca_dat, #<-
                              method = 'nearest', distance='glm',
                              replace=TRUE)
  
  sca_matches <- get_matches(sca_matched.data)
  cn_sca_psm <- sca_matches[sca_matches$RANDASSIGN==0,][,7:16]
  sca_dist_mat <- cdist(cn_sca_psm, tr_psm, metric='euclidean')
  sca_dist_mat_sum <- rowSums(sca_dist_mat)
  sca_dist_df <- as.data.frame(sca_dist_mat_sum)
  sca_dist_df$sca_cn_ids <- sca_matches[sca_matches$RANDASSIGN==0,]$id
  sca_dist_df$sca_cn_maskIDs <- sca_matches[sca_matches$RANDASSIGN==0,]$MASKID
  
  #final combined dataframe
  rating <- data.frame(treated_id = trial_matches[trial_matches$RANDASSIGN==1,]$id, 
                       treated_maskID = trial_matches[trial_matches$RANDASSIGN==1,]$MASKID)
  rating <- cbind(rating, trial_dist_df)
  rating <- cbind(rating, sca_dist_df)
  
  rating$trial_dist <- unlist(lapply(as.data.frame(rating$trial_dist_mat_sum), min_max_norm))
  rating$sca_dist <- unlist(lapply(as.data.frame(rating$sca_dist_mat_sum), min_max_norm))
  
  #removing duplicate IDs
  rating_t <- rating[!duplicated(rating$cn_maskIDs),]
  rating_s <- rating[!duplicated(rating$sca_cn_maskIDs),]
  
  num_frac_th_controls <- min(num_frac_th_controls, nrow(rating_t)) #<- 
  num_frac_sca_controls <- min(num_frac_sca_controls, nrow(rating_s)) #<-
  
  #now give top X%  from trial controls and rest from SCA controls
  frac_trial_controls <- rating_t[order(rating_t$trial_dist),][0:num_frac_th_controls,]
  frac_sca_controls <- rating_s[order(rating_s$sca_dist),][0:num_frac_sca_controls,]
  
  #generate final trial and SCA controls
  final_trial_controls <- trial_dat[is.element(trial_dat$MASKID, frac_trial_controls$cn_maskIDs),] 
  final_sca_controls <- sca_dat[is.element(sca_dat$MASKID, frac_sca_controls$sca_cn_maskIDs),]
  
  return(list(final_trial_controls, final_sca_controls))
  
}

gen_fractional_controls <- function(data, input_threshold, input_control_fraction){
  
  counts <- aggregate(x = data$MASKID, by = list(data$NEWSITEID),
                      FUN = function(x) length(unique(x)))
  
  
  #applying thresholds
  site_counts <- counts[order(-counts$x), ]
  th_sites <- site_counts[site_counts$x > input_threshold, ]
  sca_sites <- site_counts[site_counts$x <= input_threshold, ]
  
  #site numbers
  num_all_sites <- length(unique(site_counts$Group.1))
  num_th_sites <- length(unique(th_sites$Group.1))
  num_sca_sites <- num_all_sites - num_th_sites

  #slicing sites by IDs
  th_data <- data[is.element(data$NEWSITEID, th_sites$Group.1), ]
  sca_data <- data[is.element(data$NEWSITEID, sca_sites$Group.1), ]

  #treatment data from threshold sites
  tr_data_from_th <- th_data[th_data$RANDASSIGN == 1, ]
  cn_data_from_th <- th_data[th_data$RANDASSIGN == 0, ]
  cn_data_from_sca <- sca_data[sca_data$RANDASSIGN == 0, ]
  
  # #numbers of treated and controls
  num_treated <- nrow(tr_data_from_th)
  num_th_controls <- nrow(cn_data_from_th)
  num_sca_controls <- nrow(cn_data_from_sca)
  
  #fraction of controls count
  num_frac_th_controls <- round(num_treated * input_control_fraction)
  num_frac_sca_controls <- num_treated - num_frac_th_controls
  
  #make sure each fraction have enough samples,
  #otherwise set size to max samples available
  if (num_frac_th_controls >= num_th_controls) {
    num_frac_th_controls <- num_th_controls
  }
  if (num_frac_sca_controls >= num_sca_controls) {
    num_frac_sca_controls <- num_sca_controls
  }
  
  #Trial Data: Trial treated + Trial Controls
  #SCA Data: Trial treated + SCA Controls
  trial_dat <- rbind(tr_data_from_th, cn_data_from_th)
  sca_dat <- rbind(tr_data_from_th, cn_data_from_sca)
  
  #get final controls: trial controls and SCA controls
  final_controls <- matching_policy_0(trial_dat, sca_dat, num_frac_th_controls, num_frac_sca_controls)
  
  final_trial_controls <- final_controls[[1]] 
  final_sca_controls <- final_controls[[2]]
  
  num_frac_th_controls <- nrow(final_trial_controls)
  num_frac_sca_controls <- nrow(final_sca_controls)
  
  #sample stats count df
  stats <- data.frame(
    Fields = c("Total Sites", "Trial Sites", "SCA Sites", 
               "Total Treated Samples", "Total Trial Controls", 
               "Total SCA Controls", "Fraction Trial Controls", 
               "Fraction SCA Controls"),
    Values = c(num_all_sites, num_th_sites, num_sca_sites,
               num_treated, num_th_controls, num_sca_controls, 
               num_frac_th_controls, num_frac_sca_controls)
  )
  
  filtered_data <- list(tr_data_from_th, final_trial_controls, 
                        final_sca_controls, th_data, sca_data)
  
return(list(stats, filtered_data))
}
  