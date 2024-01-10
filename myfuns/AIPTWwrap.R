AIPTWwrap <- function(mydat, ps, or0, or1, bootdats, mod, is_print=TRUE) { #, ps.vars, or.vars
  
  # close_form CI (influence function-based)
  point_res <- AIPTWinner(mydat, ps, or0, or1)
  est <- point_res$est
  se <- point_res$se
  ci <- c(est-1.96*se, est+1.96*se)
  
  # bootstrap CI
  res_boot <- NULL
  res_t <- NULL
  if (!is.null(bootdats)) {
    ate_boot <- lapply(bootdats, function(boot1, mod) {
      AIPTWinner(mydat=boot1$mydat, ps=boot1$ps_lst[[mod]], 
                 or0=boot1$or0_lst[[mod]], or1=boot1$or1_lst[[mod]])
    }, mod=mod)
    
    boot_est <- unlist(lapply(ate_boot, "[[", 1))
    boot_se <- unlist(lapply(ate_boot, "[[", 2))
    
    se_boot <- sqrt(var(boot_est, na.rm=T))
    ci_boot <- as.numeric(quantile(boot_est, p = c(0.025, 0.975), na.rm=T))  
    
    t_stat <- (boot_est-est) / boot_se
    c_stat <- as.numeric(quantile(abs(t_stat), probs=0.95))
    ci_t <- c(est-c_stat*se, est+c_stat*se)
    
    res_boot <- c(est, se_boot, ci_boot)
    res_t <- c(est, se, ci_t)
  }
  
  res <- c(est, se, ci)
  if (is_print) {
    print(c("AIPTW", round(res,4)))
    if (!is.null(bootdats)) {
      print(c("AIPTW-boot", round(res_boot,4)))
      print(c("AIPTW-t", round(res_t,4)))
    }
  }
  
  return(list(res=res,res_boot=res_boot,res_t=res_t))
}


AIPTWinner <- function(mydat, ps, or0, or1) {
  Q1W <- or1
  Q0W <- or0
  gW <- ps
  
  # point est
  EY1aiptw <- mean((mydat$A) * (mydat$Y - Q1W) / gW + Q1W)
  EY0aiptw <- mean((1 - mydat$A) * (mydat$Y - Q0W) / (1 - gW) + Q0W)
  est <- EY1aiptw - EY0aiptw
  
  # SE (influence function-based)
  n <- nrow(mydat)
  D1 <- (mydat$A) * (mydat$Y - Q1W) / gW + Q1W - EY1aiptw
  D0 <- (1 - mydat$A) * (mydat$Y - Q0W) / (1 - gW) + Q0W - EY0aiptw 
  varHat_AIPTW <- var(D1 - D0) / n
  se <- sqrt(varHat_AIPTW)
  
  return(list(est=est, se=se))
}
