TMLEwrap <- function(mydat, ps, or0, or1, useCovar, bootdats, mod, is_print=TRUE) {
  
  covars_all <- grep(paste0("^",useCovar), names(mydat), value=TRUE)
  
  if (length(unique(mydat$Y)) != 2) {  # continuous Y
    sim_fam <- gaussian()
  } else {  # binary Y
    sim_fam <- binomial()
  }
  
  # close_form CI (influence function-based)
  tmle.fit <- tmle(Y = mydat$Y, A = mydat$A, W = mydat[, covars_all], family = sim_fam$family,
                   Q = cbind(or0, or1), g1W = ps)
  est <- tmle.fit$estimates$ATE$psi
  se <- sqrt(tmle.fit$estimates$ATE$var.psi)
  ci <- tmle.fit$estimates$ATE$CI
  
  # bootstrap CI
  res_boot <- NULL
  res_t <- NULL
  if (!is.null(bootdats)) {
    ate_boot <- lapply(bootdats, function(boot1, covars_all, family, mod) {
      TMLEinner(mydat=boot1$mydat, ps=boot1$ps_lst[[mod]], 
                or0=boot1$or0_lst[[mod]], or1=boot1$or1_lst[[mod]], 
                covars_all=covars_all, family=family)
    }, covars_all=covars_all, family=sim_fam$family, mod=mod)
    
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
    print(c("TMLE", round(res,4)))
    if (!is.null(bootdats)) {
      print(c("TMLE-boot", round(res_boot,4)))
      print(c("TMLE-t", round(res_t,4)))
    }
  }
  
  return(list(res=res,res_boot=res_boot,res_t=res_t))
}


TMLEinner <- function(mydat, covars_all, family, ps, or0, or1) {
  tmle.fit <- tmle(Y = mydat$Y, A = mydat$A, W = mydat[, covars_all], family = family,
                   Q = cbind(or0, or1), g1W = ps)
  est <- tmle.fit$estimates$ATE$psi
  se <- sqrt(tmle.fit$estimates$ATE$var.psi)
  return(list(est=est, se=se))
}


# tmle_se_ci <-TMLEbootCI(bootdats, covars_all, family=sim_fam$family, mod)
# tmle_se <- tmle_se_ci[1]
# tmle_ci <- c(tmle_se_ci[2], tmle_se_ci[3])

# TMLEbootCI <- function(bootdats, covars_all, family, mod){
#   ate_boot <- unlist(lapply(bootdats, function(boot1, covars_all, family, mod) {
#     TMLEinner(mydat=boot1$mydat, ps=boot1$ps_lst[[mod]], or0=boot1$or0_lst[[mod]], or1=boot1$or1_lst[[mod]], 
#               covars_all=covars_all, family=family)
#   }, covars_all=covars_all, family=family, mod=mod))
#   ate_se <- sqrt(var(ate_boot, na.rm=T))
#   ate_ci <- as.numeric(quantile(ate_boot, p = c(0.025, 0.975), na.rm=T))
#   return(c(ate_se, ate_ci))
# }

## ate_boot <- replicate(bootT, TMLE1boot(mydat, covars_all, family, ps.mod, or.mod, is_oracle))
# 
# TMLE1boot <- function(mydat, covars_all, family, ps.mod, or.mod, is_oracle){
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
#   fits <- FitPSnOR(mydat, ps.vars, or.vars, ps.mod=ps.mod, or.mod=or.mod)
#   ps <- fits$myps
#   or0 <- fits$myor0
#   or1 <- fits$myor1
# 
#   return(TMLEinner(mydat=mydat, covars_all=covars_all,
#                    family=family, ps=ps, or0=or0, or1=or1))
# }


## TMLE in package
# tmle.fit <- tmle(Y = mydat$Y, A = mydat$A, W = mydat[, covars_all], family = sim_fam$family,
#                  Q = cbind(or0, or1), g1W = ps) # Q.SL.library = "SL.glm", g.SL.library = "SL.glm"
# tmle_est <- tmle.fit$estimates$ATE$psi
# tmle_se <- sqrt(tmle.fit$estimates$ATE$var.psi)
# tmle_ci <- tmle.fit$estimates$ATE$CI
# res <- c(tmle_est, tmle_se, tmle_ci)
# print(c("TMLE", round(res,4)))


## TMLE from scratch (no SuperLearner)
# Y_ori <- mydat$Y
# mydat$Y_ori <- Y_ori
# 
# if (length(unique(Y_ori)) != 2) {  # continuous Y (need to scale to [0,1] and later scale back)
#   Y_min <- min(Y_ori)
#   Y_max <- max(Y_ori)
#   mydat$Y <- (mydat$Y - Y_min) / (Y_max - Y_min)
#   sim_fam <- gaussian()
# } else {  # binary Y
#   sim_fam <- binomial()
# }
#
# or.fml <- as.formula(paste("Y ~ A + ", paste(c(or.vars), collapse = "+")))
# or.fit <- glm(or.fml, family=binomial(), data=mydat)
# QAW <- predict(or.fit, type = "response")
# Q1W <- predict(or.fit, newdata=data.frame(A=1, mydat[, or.vars]), type="response")
# Q0W <- predict(or.fit, newdata=data.frame(A=0, mydat[, or.vars]), type="response")
# 
# ps.fml <- as.formula(paste("A ~ ", paste(c(ps.vars), collapse = "+")))
# ps.fit <- glm(ps.fml, family=binomial(), data=mydat)
# gW <- predict(ps.fit, type="response")
# 
# H1W <- (mydat$A / gW)
# H0W <- (1 - mydat$A) / (1 - gW)
# epsilon <- coef(glm(mydat$Y ~ -1 + H0W + H1W + offset(qlogis(QAW)), family=binomial()))
# 
# Q0W_1 <- plogis(qlogis(Q0W) + epsilon[1] / (1 - gW))
# Q1W_1 <- plogis(qlogis(Q1W) + epsilon[2] / gW)
# tmle_est <- mean(Q1W_1 - Q0W_1)
# EY1tmle1 <- mean(Q1W_1)
# EY0tmle1 <- mean(Q0W_1)
# MORtmle1 <- (EY1tmle1 * (1 - EY0tmle1)) / ((1 - EY1tmle1) * EY0tmle1)
# 
# D1 <- mydat$A/gW*(mydat$Y - Q1W_1) + Q1W_1 - EY1tmle1
# D0 <- (1 - mydat$A)/(1 - gW)*(mydat$Y - Q0W_1) + Q0W_1 - EY0tmle1
# EIC <- D1 - D0
# 
# n <- nrow(mydat)
# varHat.IC <- var(EIC)/n
# 
# tmle_se <- sqrt(varHat.IC)
# tmle_ci <- c(tmle_est - 1.96*tmle_se, tmle_est + 1.96*tmle_se)
# 
# if (length(unique(Y_ori)) != 2) {  # continuous Y (need to scale to [0,1] and later scale back)
#   tmle_est <- (Y_max - Y_min) * tmle_est
#   tmle_se <- (Y_max - Y_min) * tmle_se
#   tmle_ci <- (Y_max - Y_min) * tmle_ci
# }
# 
# res <- c(tmle_est, tmle_se, tmle_ci)
# print(res)
