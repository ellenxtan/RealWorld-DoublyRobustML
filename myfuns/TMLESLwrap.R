TMLESLwrap <- function(mydat, useCovar, bootdats=NULL, is_print=TRUE) { #, mod, ps, or0, or1
  
  covars_all <- grep(paste0("^",useCovar), names(mydat), value=TRUE)
  
  if (length(unique(mydat$Y)) != 2) {  # continuous Y
    sim_fam <- gaussian()
  } else {  # binary Y
    sim_fam <- binomial()
  }
  
  # close_form CI (influence function-based)
  tmle.fit <- tmle(Y = mydat$Y, A = mydat$A, W = mydat[, covars_all], family = sim_fam$family)
  est <- tmle.fit$estimates$ATE$psi
  se <- sqrt(tmle.fit$estimates$ATE$var.psi)
  ci <- tmle.fit$estimates$ATE$CI
  
  # bootstrap CI
  res_boot <- NULL
  res_t <- NULL
  if (!is.null(bootdats)) {
    ate_boot <- lapply(bootdats, function(boot1, covars_all, family) {
      TMLESLinner(mydat=boot1$mydat, covars_all=covars_all, family=family)
    }, covars_all=covars_all, family=sim_fam$family)
    
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
    print(c("TMLESL", round(res,4)))
    if (!is.null(bootdats)) {
      print(c("TMLESL-q", round(res_boot,4)))
      print(c("TMLESL-t", round(res_t,4)))
    }
  }
  
  return(list(res=res,res_boot=res_boot,res_t=res_t))
}


TMLESLinner <- function(mydat, covars_all, family) {#, ps, or0, or1
  tmle.fit <- tmle(Y = mydat$Y, A = mydat$A, W = mydat[, covars_all], family = family)#,
  # Q = cbind(or0, or1), g1W = ps) # Q.SL.library = "SL.glm", g.SL.library = "SL.glm"
  # print(tmle.fit)
  est <- tmle.fit$estimates$ATE$psi
  se <- sqrt(tmle.fit$estimates$ATE$var.psi)
  return(list(est=est, se=se))
}

# tmle_se_ci <-TMLESLbootCI(bootdats, covars_all, family=sim_fam$family) #, mod
# tmle_se <- tmle_se_ci[1]
# tmle_ci <- c(tmle_se_ci[2], tmle_se_ci[3])

# TMLESLbootCI <- function(bootdats, covars_all, family){ #, mod
#   ate_boot <- unlist(lapply(bootdats, function(boot1, covars_all, family) { #, mod
#     TMLESLinner(mydat=boot1$mydat, covars_all=covars_all, family=family) 
#     #ps=boot1$ps_lst[[mod]], or0=boot1$or0_lst[[mod]], or1=boot1$or1_lst[[mod]], 
#   }, covars_all=covars_all, family=family)) #, mod=mod
#   ate_se <- sqrt(var(ate_boot, na.rm=T))
#   ate_ci <- as.numeric(quantile(ate_boot, p = c(0.025, 0.975), na.rm=T))
#   return(c(ate_se, ate_ci))
# }

