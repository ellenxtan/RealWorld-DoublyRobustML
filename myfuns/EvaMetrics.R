# metrics: bias (absolute+relative), MSE, coverage, width, variance ratio, power
EvaMetrics <- function(tau, res) {
  res <- as.data.frame(res)
  names(res) <- c("est", "se", "lb", "ub")
  res <- na.omit(res)
  
  bias <- mean(res$est - tau) #mean(abs(res$est - tau))
  # bias.rel <- mean(abs(res$est - tau) / tau * 100)
  mse <- mean((res$est - tau)^2)
  cov <- mean(res$lb<= tau & tau<=res$ub)
  width <- mean(res$ub - res$lb)
  var.ratio <- mean((res$se)^2) / var(res$est) # ratio of model-based variance versus simulation variance
  power <- 1 - mean(res$lb<= 0 & 0<=res$ub)  # 1 - type2error for ATE=3 | type1error for ATE=0
  power.lb <- power - 1.96*sqrt(power*(1-power) / nrow(res))
  power.ub <- power + 1.96*sqrt(power*(1-power) / nrow(res))
  
  mets <- c(bias, mse, cov, width, var.ratio, power, power.lb, power.ub) # bias.abs, bias.rel
  print(mets, digits=3)
  
  return(list(res=res, mets=mets))
}
