IMPwrap <- function(or0, or1, bootdats, mod) {
  est <- mean(or1 - or0)
  
  ate_boot <- unlist(lapply(bootdats, function(boot1, mod) {
    mean(boot1$or1_lst[[mod]] - boot1$or0_lst[[mod]])
  }, mod=mod))
  se <- sqrt(var(ate_boot, na.rm=T))
  ci <- as.numeric(quantile(ate_boot, p = c(0.025, 0.975), na.rm=T))
  
  res <- c(est, se, ci)
  print(c("IMP", round(res,4)))
  return(res)
}
