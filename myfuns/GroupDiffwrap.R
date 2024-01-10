GroupDiffwrap <- function(mydat, bootdats) {
  est <- GroupDiffinner(mydat)
  
  se_ci <- GroupDiffbootCI(bootdats)
  se <- se_ci[1]
  ci <- c(se_ci[2], se_ci[3])
  
  res <- c(est, se, ci)
  print(c("GrpDiff", round(res,4)))
  return(res)
}


GroupDiffbootCI <- function(bootdats){
  ate_boot <- unlist(lapply(bootdats, function(boot1) {
    GroupDiffinner(mydat=boot1$mydat)
  }))
  ate_se <- sqrt(var(ate_boot, na.rm=T))
  ate_ci <- as.numeric(quantile(ate_boot, p = c(0.025, 0.975), na.rm=T))
  return(c(ate_se, ate_ci))
}


GroupDiffinner <- function(mydat) {
  mydat$A <- factor(mydat$A)
  dat0 <- mydat[which(mydat$A==0),]
  dat1 <- mydat[which(mydat$A==1),]
  est <- mean(dat1$Y) - mean(dat0$Y)
  return(est)
}
