#### variable selection by lasso
# Lambda.min is more conservative (negligible bias but larger variance)
# lambda.1se is more aggressive (less variance but possibly larger bias)
VarSelect <- function(mydat, useCovar, is_homo, real=FALSE) {
  
  if (real) {
    covars_all <- names(mydat)[which(!names(mydat) %in% c("A","Y"))]
  } else { # simulation
    covars_all <- grep(paste0("^",useCovar), names(mydat), value=TRUE)
  }
  # print(covars_all)
  
  if (length(unique(mydat$Y)) != 2) {  # continuous Y
    sim_fam <- gaussian()
  } else {  # binary Y
    sim_fam <- binomial()
  }
  
  if (is_homo) {
    
    lasso.fit <- cv.glmnet(x=as.matrix(mydat[, covars_all]), y=mydat$Y, family = sim_fam$family)
    coef_min = coef(lasso.fit, s = lasso.fit$lambda.min)
    index_min = which(coef_min != 0) - 1
    index_min = index_min[-1]  # remove intercept
    coef_1se = coef(lasso.fit, s = lasso.fit$lambda.1se)
    index_1se = which(coef_1se != 0) - 1
    index_1se = index_1se[-1]  # remove intercept
    if (real) {
      vars_min0 <- covars_all[index_min]
      vars_1se0 <- covars_all[index_1se]
    } else { # simulation
      if (length(index_min)>0) {
        vars_min0 <- paste0(useCovar, index_min)
      } else {
        vars_min0 <- NULL
      }
      if (length(index_1se)>0) {
        vars_1se0 <- paste0(useCovar, index_1se)
      } else {
        vars_1se0 <- NULL
      }
    }
    cat("lambda.min:", vars_min0, "\n")
    cat("lambda.1se:", vars_1se0, "\n")
    
    return(list(vars_min0=vars_min0, vars_1se0=vars_1se0))
    
    
  } else {
    
    dat0 <- mydat[which(mydat$A==0), ]
    dat1 <- mydat[which(mydat$A==1), ]
    
    ## fit0
    lasso.fit <- cv.glmnet(x=as.matrix(dat0[, covars_all]), y=dat0$Y, family = sim_fam$family)
    coef_min = coef(lasso.fit, s = lasso.fit$lambda.min)
    index_min = which(coef_min != 0) - 1
    index_min = index_min[-1]  # remove intercept
    coef_1se = coef(lasso.fit, s = lasso.fit$lambda.1se)
    index_1se = which(coef_1se != 0) - 1
    index_1se = index_1se[-1]  # remove intercept
    if (real) {
      vars_min0 <- covars_all[index_min]
      vars_1se0 <- covars_all[index_1se]
    } else { # simulation
      if (length(index_min)>0) {
        vars_min0 <- paste0(useCovar, index_min)
      } else {
        vars_min0 <- NULL
      }
      if (length(index_1se)>0) {
        vars_1se0 <- paste0(useCovar, index_1se)
      } else {
        vars_1se0 <- NULL
      }
    }
    cat("lambda.min0:", vars_min0, "\n")
    cat("lambda.1se0:", vars_1se0, "\n")
    
    ## fit1
    lasso.fit <- cv.glmnet(x=as.matrix(dat1[, covars_all]), y=dat1$Y, family = sim_fam$family)
    coef_min = coef(lasso.fit, s = lasso.fit$lambda.min)
    index_min = which(coef_min != 0) - 1
    index_min = index_min[-1]  # remove intercept
    coef_1se = coef(lasso.fit, s = lasso.fit$lambda.1se)
    index_1se = which(coef_1se != 0) - 1
    index_1se = index_1se[-1]  # remove intercept
    if (real) {
      vars_min1 <- covars_all[index_min]
      vars_1se1 <- covars_all[index_1se]
    } else { # simulation
      if (length(index_min)>0) {
        vars_min1 <- paste0(useCovar, index_min)
      } else {
        vars_min1 <- NULL
      }
      if (length(index_1se)>0) {
        vars_1se1 <- paste0(useCovar, index_1se)
      } else {
        vars_1se1 <- NULL
      }
    }
    cat("lambda.min1:", vars_min1, "\n")
    cat("lambda.1se1:", vars_1se1, "\n")
    
    return(list(vars_min0=vars_min0, vars_1se0=vars_1se0,
                vars_min1=vars_min1, vars_1se1=vars_1se1))
    
  }
  
  
}
