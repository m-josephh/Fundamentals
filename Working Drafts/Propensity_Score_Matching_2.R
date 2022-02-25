psm_manual <- function(subset.data , quantile_breaks = 0.05){
  # (c) Michael Lahm
  print("Make sure columns are in this order: AGE, GENDER, RAF, ID, HAI FLAG")
  colnames(subset.data) <- c("Age", "Gender", "RAF", "ID", "HAI_Flag")
  treat <- unique(subset.data[subset.data$HAI_Flag_2 == 1, "ID"])
  pm <- glm(HAI_Flag ~ Age + Gender + RAF, data = subset.data, family = binomial)
  
  control_scores <- predict(pm)[subset.data$HAI_Flag == 0]
  trt_scores <- predict(pm)[subset.data$HAI_Flag == 1]
  
  # Histogram Summaries
  par(mfrow = c(1,2))
  hist(control_scores, main = "Control Scores; Before Matching", xlab = "Propensity Scores", col = gray(0.8))
  hist(trt_scores, main = "Treatment Scores; Before Matching", xlab = "Propensity Scores", col = "steel blue")
  
  
  # Creating the buckets for matching
  qq <- quantile(trt_scores, seq(0,1,quantile_breaks))
  control_bucket <- cut(control_scores, breaks = qq)
  qtab <- table(control_bucket)
  control_include <- rep(0, length(control_scores))
  
  # Matching Proper
  for(ii in unique(control_bucket[!is.na(control_bucket)])){
    rn <- which(control_bucket == ii)
    sp <- sample(rn, min(qtab)) # Sampling based on the smallest frequency in the controll group
    control_include[sp] <-1
  }
  
  control.pre <- subset.data[subset.data$HAI_Flag==0,]
  control_id <- control.pre$ID[control_include == 1]
  
  retval_ids <- subset.data[subset.data$ID %in% c(control_id, treat),]
  
  # Histogram after matching
  pm.post_match <- glm(HAI_Flag ~ Age + Gender + RAF, data = retval_ids, family = binomial)
  post_control_scores <- predict(pm.post_match)[retval_ids$HAI_Flag ==0]
  post_treat_scores <- predict(pm.post_match)[retval_ids$HAI_Flag ==1]
  par(mfrow = c(1,2))
  hist(post_control_scores, main = "Control Scores; Post-matching", xlab = "Propensity Scores")
  hist(post_treat_scores, main = "Treatment Scores; Post-matching", xlab = "Propensity Scores")
  
  # Summary of data
  before_matching <- list(
    mean_propensity_score_diff = abs(mean(control_scores) - mean(trt_scores)),
    sd_diff = abs(sd(control_scores) - sd(trt_scores)),
    mean_RAF_diff = abs(mean(subset.data$RAF[subset.data$HAI_Flag == 0]) - mean(subset.data$RAF[subset.data$HAI_Flag==1]))
  )
  after_matching <- list(
    mean_propensity_score_diff = abs(mean(post_control_scores) - mean(post_treat_scores)),
    sd_diff = abs(sd(post_control_scores) - sd(post_treat_scores)),
    mean_RAF_diff = abs(mean(retval_ids$RAF[retval_ids$HAI_Flag == 0]) - mean(retval_ids$RAF[retval_ids$HAI_Flag == 1]))
  )
  summaries <- rbind(before_matching, after_matching)
  
  retval <- list()
  retval$data <- retval_ids
  retval$summaries <- summaries
  
  return(retval)
  
}