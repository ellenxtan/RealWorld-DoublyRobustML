IPTWwrap <- function(mydat, ps, bootdats, mod) {
  
  if (is.null(bootdats)) {
    # https://bolinwu.blog/post/an-overview-of-causal-inference-part-3-inverse-probability-of-treatment-weighting-iptw/
    
    # creat weights from ps
    weight = ifelse(mydat$A == 1, 1/ps, 1/ (1-ps))
    # cat("the maximum weight is", max(weight), "\n",
    #     "the minimum weight is", min(weight), "\n")
    
    # fit a marginal structural model 
    msm = svyglm(Y~A, design = svydesign(~1, weights = ~weight, data = mydat))
    smry <- summary(msm)
    iptw_est <- smry$coefficients[2,1]
    iptw_se <- smry$coefficients[2,2]
    iptw_ci <- c(confint(msm)[2,1], confint(msm)[2,2])
    
  } else {
    iptw_est <- IPTWinner(mydat, ps)
    
    # bootstrap ci
    iptw_se_ci <-IPTWbootCI(bootdats, mod)
    iptw_se <- iptw_se_ci[1]
    iptw_ci <- c(iptw_se_ci[2], iptw_se_ci[3])
  }
  
  res <- c(iptw_est, iptw_se, iptw_ci)
  print(c("IPTW", round(res,4)))
  return(res)
}


IPTWbootCI <- function(bootdats, mod){
  ate_boot <- unlist(lapply(bootdats, function(boot1, mod) {
    IPTWinner(mydat=boot1$mydat, ps=boot1$ps_lst[[mod]])
  }, mod=mod))
  ate_se <- sqrt(var(ate_boot, na.rm=T))
  ate_ci <- as.numeric(quantile(ate_boot, p = c(0.025, 0.975), na.rm=T))
  return(c(ate_se, ate_ci))
}


IPTWinner <- function(mydat, ps){
  A <- mydat$A
  Y <- mydat$Y
  g1hat <- ps
  # print(range(g1hat))
  g1hat <- ifelse(g1hat==1.0, 0.99999999, g1hat)
  g1hat <- ifelse(g1hat==0.0, 0.00000001, g1hat)
  # print(range(g1hat))
  ## Hajek IPW
  # est <- sum(A * Y /(g1hat))/sum(A/g1hat) - sum((1-A)*Y/(1- g1hat))/sum((1-A)/(1-g1hat))
  ## Horvitz-Thompson IPW
  est <- mean(A * Y /g1hat) - mean((1-A)*Y / (1- g1hat))
  return(est)
}



## ate_boot <- replicate(bootT, IPTW1boot(mydat, ps.mod, is_oracle))
# 
# IPTW1boot <- function(mydat, ps.mod, is_oracle){
#   #### boot data
#   resamp_idx <- sample(1:nrow(mydat), replace = TRUE)
#   mydat_boot <- mydat[resamp_idx, , drop = FALSE]
#   
#   #### variable selection
#   sink("NUL")
#   var_select <- VarSelect(mydat_boot, is_oracle)
#   sink()
#   mydat <- var_select$mydat
#   ps.vars <- var_select$ps.vars
#   or.vars <- var_select$or.vars
#   
#   #### PS & OR
#   fits <- FitPSnOR(mydat, ps.vars, or.vars, ps.mod=ps.mod, or.mod=NULL)  # no or.mod for IPTW
#   ps <- fits$myps
#   
#   return(IPTWinner(mydat, ps))
# }
