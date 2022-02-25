psm_manual <- function(dat, 
                       response, 
                       predictors ,
                       quantile_breaks = 0.08, 
                       positive = 1, 
                       negative = 0 ){
  
  # Renaming ID variable
  names(dat)[1] <- "ID" # The ID is typically in the first column
  
  # renaming the response variable
  names(dat)[names(dat) %in% response] <- "response"
  
  # defining the response and the independent variables as a formula object
  predictors <- paste("response", paste0(predictors, collapse = "+"), sep= "~")
  formula_glm <- as.formula(predictors)
  
  
  treat <- unique(dat[dat$response == positive, "ID"])
  pm <- glm(formula_glm, data = dat, family = binomial)
  
  control_scores <- predict(pm)[dat$response == negative]
  trt_scores <- predict(pm)[dat$response == positive]
  
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
  
  control.pre <- dat[dat$response== negative,]
  control_id <- control.pre$ID[control_include == positive]
  
  retval_ids <- dat[dat$ID %in% c(control_id, treat),]
  
  # Histogram after matching
  pm.post_match <- glm(formula_glm, data = retval_ids, family = binomial)
  post_control_scores <- predict(pm.post_match)[retval_ids$response == negative]
  post_treat_scores <- predict(pm.post_match)[retval_ids$response == positive]
  par(mfrow = c(1,2))
  hist(post_control_scores, main = "Control Scores; Post-matching", xlab = "Propensity Scores")
  hist(post_treat_scores, main = "Treatment Scores; Post-matching", xlab = "Propensity Scores")
  

  # resetting the name of the response variable
  names(retval_ids)[names(retval_ids) %in% "response"] <- response
  
  return(retval_ids)
  
}