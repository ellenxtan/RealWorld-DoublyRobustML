# very slow in high dim (2 versions: CV or not)
OHALwrap <- function(mydat, useCovar) {
  
  covars_all <- grep(paste0("^",useCovar), names(mydat), value=TRUE)
  W <- mydat[, covars_all]
  A <- mydat$A
  Y <- mydat$Y
  
  if (length(unique(Y)) == 2) {  # binary Y
    sim_fam <- binomial()
  } else {  # continuous Y
    sim_fam <- gaussian()
  }
  
  # call ohal_nuisance to get nuisance estimates
  fit <- ohal_nuisance(W = W, A = A, Y = Y, V = 5,
                       outcome_family = sim_fam$family,
                       lambda_seq = exp(seq(-1, -13, length=100)))
  
  # OHAL
  drtmle_ohal_fit <- drtmle(W = W, A = A, Y = Y, 
                            a_0 = c(0,1), tolg = 1e-4,
                            Qn = list(fit$Q0W, fit$Q1W),
                            # Qsteps = 1, # uncomment for one-step targeting
                            gn = list(1 - fit$G1W_OHAL, fit$G1W_OHAL),
                            SL_gr = "SL.hal9001", guard = "g",
                            use_future = TRUE) # use_future: parallel computing
  
  drtmle_ohal_ci <- ci(drtmle_ohal_fit, contrast = c(-1, 1))
  ohal_est <- drtmle_ohal_ci$drtmle[1]
  ohal_lb <- drtmle_ohal_ci$drtmle[2]
  ohal_ub <- drtmle_ohal_ci$drtmle[3]
  ohal_se <- (drtmle_ohal_ci$drtmle[3] - drtmle_ohal_ci$drtmle[2]) / 2 / 1.96
  
  res <- c(ohal_est, ohal_se, ohal_lb, ohal_ub)
  print(c("OHAL", round(res,4)))
  
  # cross-validated OHAL
  drtmle_ohal_fit_cv <- drtmle(W = W, A = A, Y = Y,
                               a_0 = c(0, 1), tolg = 1e-4,
                               Qn = list(fit$cv_Q0W, fit$cv_Q1W),
                               # Qsteps = 1, # uncomment for one-step targeting
                               gn = list(1 - fit$cv_G1W_OHAL, fit$cv_G1W_OHAL),
                               cvFolds = fit$fold_vec,
                               SL_gr = "SL.hal9001",
                               guard = "g", returnModels = FALSE, verbose = FALSE,
                               use_future = TRUE) # use_future: parallel computing
  
  tmp <- ci(drtmle_ohal_fit_cv, contrast = c(-1, 1), est = "aiptw_c")
  # get a confidence interval about drtmle + ohal using cross-validated standard errors
  ohal_se <- (tmp$aiptw_c[3] - tmp$aiptw_c[2]) / 2 / 1.96
  ohal_est <- drtmle_ohal_fit$drtmle$est[2] - drtmle_ohal_fit$drtmle$est[1]
  ohal_lb <- ohal_est - 1.96 * ohal_se
  ohal_ub <- ohal_est + 1.96 * ohal_se
  
  res_cv <- c(ohal_est, ohal_se, ohal_lb, ohal_ub)
  print(c("OHALCV", round(res_cv,4)))
  
  return(list(ohal=res, ohalcv=res_cv))
}



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#' OHAL original codes
#' 
#' Fit a custom version of HAL for outcome adaptive TMLE estimation. The function takes in data 
#' and outputs nuisance estimates evaluated at the observed data for the outcome regression 
#' fit by HAL (and cross-validated versions of that), the propensity score fit by OHAL (and cross-validated versions), 
#' and the propensity score fit by HAL (and cross-validated versions).  
#' 
#' @param W A data.frame of predictors
#' @param A A numeric binary (0/1) treatment vector
#' @param Y A numeric outcome vector
#' @param V The number of cross-validation folds (must be at least 3 or \code{glmnet} will complain)
#' @param outcome_family "gaussian" or "binomial" (implies loss function to use)
#' @param lambda_seq A sequence of lambda values to be passed to \code{glmnet}.
#' In general, the \code{glmnet} defaults do not work well for the purposes of HAL.
#' The default is something that has yielded good results in the past, but needs
#' to be validated for each simulation.
#' @param penalty_threshold Soft threshold for penalty term for OHAL
#' @param penalty_gamma OHAL tuning parameter
#' @param newdata A data.frame of new values of W on which to provide predictions 
#' @return A data.frame. 
#' \describe{
#'   \item{\code{QAW}}{the HAL fit at CV-selected lambda for Qbar evaluated at 
#'                     (A_i, W_i), for i = 1,...,n ; \code{Q1W} = the HAL fit at 
#'                     CV-selected lambda for Qbar evaluated at (1, W_i)}
#'   \item{\code{Q0W}}{the HAL fit at CV-selected lambda for Qbar evaluated at 
#'                     (0, W_i)}
#'   \item{\code{G1W_HAL}}{the HAL fit at CV-selected lambda for G evaluated at 
#'                         W_i}
#'   \item{\code{G1W_OHAL}}{the OHAL fit at CV-selected lambda for G evaluated 
#'                          at W_i}
#'   \item{\code{cv_}}{the cross-validated HAL fit at CV-selected lambda of the 
#'                     above nuisance parameters
#' }
#' @importFrom hal9001 enumerate_basis make_design_matrix make_copy_map
#' }
ohal_nuisance <- function(W, A, Y, V = 10, outcome_family = "binomial",
                          penalty_threshold = 1e-10, penalty_gamma = 1, newdata = NULL, 
                          lambda_seq = exp(seq(0, -13, length=2000)),
                          num_cores=0){
  
  n <- length(A)
  # make a vector of cv-folds
  chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  fold_idx <- chunk2(seq_len(n), V)
  fold_vec <- rep(seq_len(V), unlist(lapply(fold_idx, length)))
  fold_idx <- unlist(fold_idx, use.names = FALSE)
  
  # cross-validation for OR
  cv_out <- sapply(seq_len(V), one_hal_or, 
                   fold_vec = fold_vec, 
                   W = W, A = A, Y = Y, 
                   lambda_seq = lambda_seq, 
                   outcome_family = outcome_family,
                   simplify = FALSE)
  
  # cv risks of OR
  or_risks <- colMeans(Reduce("rbind",lapply(cv_out, "[[", "riskQ")))
  # take smallest lambda that has smallest risk
  lambda_or_idx <- which.min(or_risks)[1]
  # re-fit on full data
  full_hal <- one_hal_or(fold = NULL, fold_vec = NULL, 
                         W = W, A = A, Y = Y, 
                         lambda_seq = lambda_seq, newdata = newdata, 
                         outcome_family = outcome_family)
  
  # QAW at cv selected lambda
  QAW_cvselect <- full_hal$QAW[, lambda_or_idx]
  # Q1W at cv selected lambda
  Q1W_cvselect <- full_hal$Q1W[, lambda_or_idx]
  # Q0W at cv selected lambda
  Q0W_cvselect <- full_hal$Q0W[, lambda_or_idx]
  # cross-validated Q1W at cv selected lambda
  QAW_cv_cvselect <- Reduce("c", lapply(cv_out, function(x){
    x$QAW[, lambda_or_idx]
  })) 
  # cross-validated Q1W at cv selected lambda
  Q1W_cv_cvselect <- Reduce("c", lapply(cv_out, function(x){
    x$Q1W[, lambda_or_idx]
  }))
  # cross-validated Q0W at cv selected lambda
  Q0W_cv_cvselect <- Reduce("c", lapply(cv_out, function(x){
    x$Q0W[, lambda_or_idx]
  }))
  
  # get cv selected OR coefficients for each fold
  or_coef_list <- lapply(cv_out, function(x){
    x$or_coef[,lambda_or_idx]
  })
  
  # PS
  cv_out_ps <- mapply(fold = seq_len(V), or_coef = or_coef_list,
                      FUN = one_hal_ps, 
                      MoreArgs = list(fold_vec = fold_vec, 
                                      W = W, A = A, penalty_gamma = penalty_gamma,
                                      lambda_seq = lambda_seq, 
                                      penalty_threshold = penalty_threshold), 
                      SIMPLIFY = FALSE)
  
  # cv risks of OHAL PS
  ohal_risks <- colMeans(Reduce("rbind",lapply(cv_out_ps, "[[", "ohal_risk")))
  hal_risks <- colMeans(Reduce("rbind",lapply(cv_out_ps, "[[", "hal_risk")))
  
  # take smallest lambda that has smallest risk
  ohal_lambda_ps_idx <- which.min(ohal_risks)[1]
  hal_lambda_ps_idx <- which.min(hal_risks)[1]
  
  # re-fit on full data
  full_hal_ps <- one_hal_ps(fold = NULL, fold_vec = NULL, 
                            W = W, A = A, penalty_gamma = penalty_gamma,
                            lambda_seq = lambda_seq, newdata = newdata, 
                            or_coef = full_hal$or_coef[,lambda_or_idx],
                            penalty_threshold = penalty_threshold)
  # G1W at cv selected lambda
  G1W_OHAL_cvselect <- full_hal_ps$G1W_OHAL[, ohal_lambda_ps_idx]
  G1W_HAL_cvselect <- full_hal_ps$G1W_HAL[, hal_lambda_ps_idx]
  
  # cross-validated G1W at cv selected lambda
  G1W_OHAL_cv_cvselect <- Reduce("c", lapply(cv_out_ps, function(x){
    x$G1W_OHAL[, ohal_lambda_ps_idx]
  }))
  G1W_HAL_cv_cvselect <- Reduce("c", lapply(cv_out_ps, function(x){
    x$G1W_HAL[, hal_lambda_ps_idx]
  }))
  
  # format outcome
  out <- data.frame(QAW = QAW_cvselect,
                    Q1W = Q1W_cvselect,
                    Q0W = Q0W_cvselect,
                    G1W_OHAL = G1W_OHAL_cvselect,
                    G1W_HAL = G1W_HAL_cvselect,
                    cv_QAW = QAW_cv_cvselect, 
                    cv_Q1W = Q1W_cv_cvselect,
                    cv_Q0W = Q0W_cv_cvselect,
                    cv_G1W_OHAL = G1W_OHAL_cv_cvselect,
                    cv_G1W_HAL = G1W_HAL_cv_cvselect,
                    fold_vec = fold_vec)
  
  return(out)
}

#' A helper function for fitting a single fold of HAL. Can be used both for
#' fitting cross-validated HAL and fitting HAL to full data. 
#' 
#' @param fold The validation fold (a numeric 1:V)
#' @param fold_vec A vector of which fold each observation is in. 
#' @param outcome_family "gaussian" or "binomial"
#' @param parametric_fits Leave set to FALSE
#' @param newdata A data.frame on which to return new predictions.
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param lambda_seq See documentation for oat_hal
#' 
#' @return See function itself to understand what it's returning. 
one_hal_or <- function(fold, fold_vec, outcome_family = "binomial",
                       parametric_fits = FALSE, newdata = NULL, 
                       W, A, Y, lambda_seq = exp(seq(-1,-13,length=10000))){
  n <- length(Y)
  n_valid <- sum(fold == fold_vec)
  n_train <- n - n_valid
  
  if(!parametric_fits){
    # if called during cross-validation
    if(!is.null(fold)){
      # only put W in here and we'll cbind A in later 
      x_fit <- as.matrix(W[fold_vec != fold, , drop = FALSE])
      a_fit <- A[fold_vec != fold]
      y_fit <- Y[fold_vec != fold]
      # if called fitting to full data
    }else{
      x_fit <- as.matrix(data.frame(W))
      a_fit <- A
      y_fit <- Y
    }
    
    # for outcome regression
    basis_list <- hal9001::enumerate_basis(x_fit, NULL)
    x_basis <- hal9001:::make_design_matrix(x_fit, basis_list)
    copy_map <- hal9001:::make_copy_map(x_basis)
    unique_columns <- as.numeric(names(copy_map))
    # subset to non-duplicated columns
    x_basis_fit <- x_basis[, unique_columns]
    # add in treatment
    basis_fit <- cbind(a_fit, x_basis_fit)
    
    # fit outcome regression 
    hal_lasso <- glmnet::glmnet(x = basis_fit, y = y_fit, 
                                family = outcome_family, lambda = lambda_seq,
                                standardize = FALSE)
    
    # predictions on validation sample
    if(!is.null(fold)){
      # for selecting lambda
      new_x_fit1 <- as.matrix(data.frame(W)[fold_vec == fold, , drop = FALSE])
      new_a <- A[fold_vec == fold]
      new_a1 <- rep(1, sum(fold_vec == fold))
      new_a0 <- rep(0, sum(fold_vec == fold))
      # for use later in TMLE
      # new_x_fit2 <- as.matrix(rbind(data.frame(A = 1, W)[fold_vec == fold, , drop = FALSE], 
      #                     data.frame(A = 0, W)[fold_vec == fold, , drop = FALSE]))
    }else if(is.null(newdata)){
      # for selecting lambda
      new_x_fit1 <- as.matrix(data.frame(W))
      new_a <- A
      new_a1 <- rep(1, n)
      new_a0 <- rep(1, n)
      # for use later in TMLE
      # new_x_fit2 <- as.matrix(rbind(data.frame(A = 1, W), 
      #                     data.frame(A = 0, W)))
      # new_x_fit3 <- as.matrix(W)
    }else{
      new_x_fit1 <- as.matrix(data.frame(newdata))
      new_a <- rep(1, nrow(newdata))
      new_a1 <- rep(1, nrow(newdata))
      new_a0 <- rep(0, nrow(newdata))
    }
    # make HAL design for getting predictions to select lambda 
    new_x_basis1 <- hal9001:::make_design_matrix(new_x_fit1, basis_list)
    new_x_basis1 <- as.matrix(new_x_basis1[, unique_columns])
    
    new_x_basisa1 <- cbind(new_a1, new_x_basis1)
    new_x_basisa0 <- cbind(new_a0, new_x_basis1)
    new_x_basis2 <- rbind(new_x_basisa1, new_x_basisa0)
    # get outcome regression predictions
    beta_hat <- as.matrix(hal_lasso$beta)
    # intercept 
    alpha_hat <- hal_lasso$a0
    # prediction matrices
    # for (A,W)
    or_pred_matrix1 <- cbind(rep(1, ifelse(is.null(fold), n, n_valid)), 
                             cbind(new_a, new_x_basis1)) %*% rbind(alpha_hat, beta_hat)
    # for rbind((1,W),(0,W))
    or_pred_matrix2 <- cbind(rep(1, ifelse(is.null(fold), 2* n, n_valid*2)), new_x_basis2) %*% rbind(alpha_hat, beta_hat)
    if(outcome_family == "binomial"){
      or_pred_matrix1 <- apply(or_pred_matrix1, 2, plogis)
      or_pred_matrix2 <- apply(or_pred_matrix2, 2, plogis)
    }
    
    # compute MSE/loglik
    if(is.null(newdata)){
      if(outcome_family == "gaussian"){
        or_risk <- apply(or_pred_matrix1, 2, function(x){
          mean((Y[fold_vec == fold] - x)^2)
        })
      }else{
        or_risk <- apply(or_pred_matrix1, 2, function(x){
          mean(ifelse(Y[fold_vec == fold] == 1, -log(x), -log(1 - x)))
        })
      }
    }else{
      or_risk <- NULL
    }
    
    # format output
    out <- list()
    out$QAW <- or_pred_matrix1
    out$riskQ <- NULL
    out$riskG <- NULL
    if(!is.null(fold)){
      out$riskQ <- or_risk
    }
    if(is.null(fold) & is.null(newdata)){
      Q1W_idx <- 1:n
      Q0W_idx <- (n+1):(2*n)
    }else if(!is.null(fold) & is.null(newdata)){
      Q1W_idx <- 1:n_valid
      Q0W_idx <- (n_valid + 1):(2*n_valid)
    }else{
      Q1W_idx <- 1:nrow(newdata)
      Q0W_idx <- (nrow(newdata) + 1):(2*nrow(newdata))
    }
    out$Q1W <- or_pred_matrix2[Q1W_idx, , drop = FALSE]
    out$Q0W <- or_pred_matrix2[Q0W_idx, , drop = FALSE]
    out$or_coef <- hal_lasso$beta
  }
  return(out)
}

#' A helper function for fitting a single fold of HAL. Can be used both for
#' fitting cross-validated HAL and fitting HAL to full data. 
#' 
#' @param fold The validation fold (a numeric 1:V)
#' @param fold_vec A vector of which fold each observation is in. 
#' @param outcome_family "gaussian" or "binomial"
#' @param W Covariates
#' @param A Treatment
#' @param lambda_seq See documentation for oat_hal
#' @param penalty_threshold Threshold for truncating outcome regression
#' coefficients
#' @param penalty_gamma OHAL tuning parameter
#' @param or_coef Vector of coefficients from HAL outcome regression
#' 
#' @return See function itself to understand what it's returning. 
one_hal_ps <- function(fold, fold_vec, W, A, 
                       or_coef, # passed in from HAL OR
                       penalty_gamma = 1,
                       penalty_threshold = 1e-5, newdata = NULL,
                       lambda_seq = exp(seq(-1,-13,length=10000))){
  n <- length(A)
  n_valid <- sum(fold == fold_vec)
  n_train <- n - n_valid
  
  # if called during cross-validation
  if(!is.null(fold)){
    # only put W in here and we'll cbind A in later 
    x_fit <- as.matrix(W[fold_vec != fold, , drop = FALSE])
    a_fit <- A[fold_vec != fold]
    # if called fitting to full data
  }else{
    x_fit <- as.matrix(data.frame(W))
    a_fit <- A
  }
  
  # for PS
  basis_list <- hal9001::enumerate_basis(x_fit, NULL)
  x_basis <- hal9001:::make_design_matrix(x_fit, basis_list)
  copy_map <- hal9001:::make_copy_map(x_basis)
  unique_columns <- as.numeric(names(copy_map))
  # subset to non-duplicated columns
  x_basis_fit <- x_basis[, unique_columns]
  
  # q coefficients
  or_coef_tmp <- abs(or_coef[-1])
  penalty.factor <- or_coef_tmp^(-penalty_gamma)
  penalty.factor[or_coef_tmp < penalty_threshold] <- Inf
  # penalty.factor <- 1 / or_coef_tmps
  # or_coef_tmp[or_coef_tmp < penalty_threshold] <- penalty_threshold
  
  if(!all(or_coef_tmp == 0)){
    # fit ohal
    ohal_lasso <- glmnet::glmnet(x = x_basis_fit, y = a_fit, 
                                 family = 'binomial', lambda = lambda_seq,
                                 standardize = FALSE, 
                                 penalty.factor = penalty.factor)
  }else{
    ohal_lasso <- list(beta = matrix(0, nrow = dim(x_basis_fit)[2], ncol = length(lambda_seq)),
                       a0 = qlogis(mean(a_fit)))
  }
  # hal lasso
  hal_lasso <- glmnet::glmnet(x = x_basis_fit, y = a_fit, 
                              family = 'binomial', lambda = lambda_seq,
                              standardize = FALSE)
  
  # predictions on validation sample
  if(!is.null(fold)){
    # for selecting lambda
    new_x_fit1 <- as.matrix(data.frame(W)[fold_vec == fold, , drop = FALSE])
    # for use later in TMLE
    # new_x_fit2 <- as.matrix(rbind(data.frame(A = 1, W)[fold_vec == fold, , drop = FALSE], 
    #                     data.frame(A = 0, W)[fold_vec == fold, , drop = FALSE]))
  }else if(is.null(newdata)){
    # for selecting lambda
    new_x_fit1 <- as.matrix(data.frame(W))
    # for use later in TMLE
    # new_x_fit2 <- as.matrix(rbind(data.frame(A = 1, W), 
    #                     data.frame(A = 0, W)))
    # new_x_fit3 <- as.matrix(W)
  }else{
    new_x_fit1 <- as.matrix(data.frame(newdata))
  }
  # make HAL design for getting predictions to select lambda 
  new_x_basis1 <- hal9001:::make_design_matrix(new_x_fit1, basis_list)
  new_x_basis1 <- as.matrix(new_x_basis1[, unique_columns])
  
  # get outcome regression predictions
  ohal_beta_hat <- as.matrix(ohal_lasso$beta)
  hal_beta_hat <- as.matrix(hal_lasso$beta)
  # intercept 
  ohal_alpha_hat <- ohal_lasso$a0
  hal_alpha_hat <- hal_lasso$a0
  # prediction matrices
  # for (A,W)
  ps_pred_matrix1 <- cbind(rep(1, ifelse(is.null(fold), n, n_valid)), new_x_basis1) %*% rbind(ohal_alpha_hat, ohal_beta_hat)
  # for rbind((1,W),(0,W))
  ps_pred_matrix2 <- cbind(rep(1, ifelse(is.null(fold), n, n_valid)), new_x_basis1) %*% rbind(hal_alpha_hat, hal_beta_hat)
  
  ps_pred_matrix1 <- apply(ps_pred_matrix1, 2, plogis)
  ps_pred_matrix2 <- apply(ps_pred_matrix2, 2, plogis)
  
  # compute MSE/loglik
  if(is.null(newdata)){
    ohal_ps_risk <- apply(ps_pred_matrix1, 2, function(x){
      mean(ifelse(A[fold_vec == fold] == 1, -log(x), -log(1 - x)))
    })
    hal_ps_risk <- apply(ps_pred_matrix2, 2, function(x){
      mean(ifelse(A[fold_vec == fold] == 1, -log(x), -log(1 - x)))
    })
  }else{
    ohal_ps_risk <- hal_ps_risk <- NULL
  }
  
  
  # format output
  out <- list()
  out$G1W_OHAL <- ps_pred_matrix1
  out$G1W_HAL <- ps_pred_matrix2
  
  out$riskQ <- NULL
  out$riskG <- NULL
  if(!is.null(fold)){
    out$ohal_risk <- ohal_ps_risk
    out$hal_risk <- hal_ps_risk
  }
  return(out)
}

#' Helper function to compute inverse probability weights
#' @param fp Prob(A = 1 | W)
#' @param fA A
#' @param fw Not used
#' @return Vector of inverse probability weights
create_weights = function(fp,fA,fw){
  fw = (fp)^(-1)
  fw[fA==0] = (1 - fp[fA==0])^(-1)
  return(fw)
}

#' Helper function to compute weighted absolute mean difference
#' @param DataM Data set
#' @param varlist Names of variables on which to balance
#' @param trt.var Identifier of treatment variable
#' @param wgt Weight
#' @param beta ?
#' @return See function to understand what it returns
wAMD_function = function(DataM,varlist,trt.var,wgt,beta){
  trt = untrt = diff_vec = rep(NA,length(beta)) 
  names(trt) = names(untrt) = names(diff_vec) = varlist
  for(jj in 1:length(varlist)){ 
    this.var = paste("w",varlist[jj],sep="") 
    DataM[,this.var] = DataM[,varlist[jj]] * DataM[,wgt] 
    trt[jj] = sum( DataM[DataM[,trt.var]==1, this.var ]) / sum(DataM[DataM[,trt.var]==1, wgt]) 
    untrt[jj] = sum(DataM[DataM[,trt.var]==0, this.var]) / sum(DataM[DataM[,trt.var]==0, wgt]) 
    diff_vec[jj] = abs( trt[jj] - untrt[jj] ) 
  } 
  wdiff_vec = diff_vec * abs(beta) 
  wAMD = c( sum(wdiff_vec))
  ret = list( diff_vec = diff_vec, wdiff_vec = wdiff_vec, wAMD = wAMD )
  return(ret) 
}

#' Helper function to estimate stabilized IPTW
#' @param fY Y
#' @param fw IPW
#' @param fA A
ATE_est = function(fY,fw,fA){
  t_ATE = fY*fw
  tt_ATE = ( ( sum(t_ATE[fA==1]) / sum(fw[fA==1]) ) - ( sum(t_ATE[fA==0]) /  sum(fw[fA==0]) ) )
  return(tt_ATE) 
}
#' Function to implement outcome adaptive lasso as proposed by shortreed
#' @param Y The outcome
#' @param A The treatment
#' @param W The covariates
#' @param V Number CV folds
#' @param family The family of the outcome regression for glm
#' @return See function to understand the return value
oal <- function(Y, A, W, V = 5, family){
  library(lqa) 
  
  Data <- as.data.frame(cbind(Y = Y, A = A, W), row.names = NULL)
  
  p = dim(W)[2]
  var.list = colnames(Data)[(1:p)+2]
  
  y.form = formula(paste("Y~."))
  lm.Y = glm(y.form,data=Data, family = family)
  betaXY = coef(lm.Y)[3:(dim(W)[2]+2)] 
  Q1W <- predict(lm.Y, newdata = data.frame(Y = Y, A = 1, W), type = "response")
  Q0W <- predict(lm.Y, newdata = data.frame(Y = Y, A = 0, W), type = "response")
  
  #----------------------------------------------------
  # set vector of possible lambda's to try
  lambda_vec = c(0.49, 0.1, 0.05, seq(0, -10, length.out = 11))
  names(lambda_vec) = as.character(lambda_vec)
  # lambda_n (n)^(gamma/2 - 1) = n^(gamma_convergence_factor)
  gamma_convergence_factor = 2
  # get the gamma value for each value in the lambda vector that corresponds to convergence factor
  gamma_vals = 2*( gamma_convergence_factor - lambda_vec + 1 )
  names(gamma_vals) = names(lambda_vec)
  
  ## Want to save ATE, wAMD and propensity score coefficients for each lambda value
  ATE = wAMD_vec = rep(NA, length(lambda_vec))
  names(ATE) = names(wAMD_vec) = names(lambda_vec)
  coeff_XA = as.data.frame(matrix(NA,nrow=dim(W)[2],ncol=length(lambda_vec)))
  names(coeff_XA) = names(lambda_vec)
  #----------------------------------------------------
  
  ######################################################################################
  #####  Run outcome adaptive lasso for each lambda value 
  ######################################################################################
  # weight model with all possible covariates included, this is passed into lasso function
  w.full.form = formula(paste("A~",paste(var.list,collapse="+")))
  
  gns <- NULL
  gns_cv <- NULL
  
  n <- length(A)
  folds <- by(sample(1:n, n), rep(1:V, length = n), list)
  
  for( lil in names(lambda_vec) ){
    tryCatch({
      il = lambda_vec[lil]
      ig = gamma_vals[lil]
      
      ### create the outcome adaptive lasso penalty with coefficient specific weights determined by outcome model
      oal_pen = adaptive.lasso(lambda=n^(il),al.weights = abs(betaXY)^(-ig) )
      ### run outcome-adaptive lasso model with appropriate penalty
      logit_oal = lqa(w.full.form, data=Data, penalty=oal_pen, family=binomial(logit))
      
      # CV
      gn_cv <- rep(-1, n)
      for(fold in folds){
        tmp <- lqa( w.full.form, data=Data[-fold, c("A", var.list)], 
                    penalty=oal_pen, family=binomial(logit))
        gn_cv[fold] <- predict(tmp, 
                               new.x = cbind(intercept=coef(tmp)[1],
                                             Data[fold,c(var.list)], row.names = NULL)
        )$mu.new
      }
      
      # generate propensity score
      Data[,paste("f.pA",lil,sep="")] = predict(logit_oal)$mu.new
      # save propensity score coefficients
      coeff_XA[var.list,lil] = coef(logit_oal)[var.list]
      # create inverse probability of treatment weights
      Data[,paste("w",lil,sep="")] = create_weights(fp=Data[,paste("f.pA",lil,sep="")],fA=Data$A)
      
      # estimate weighted absolute mean different over all covaraites using this lambda to generate weights
      wAMD_vec[lil] = wAMD_function(DataM=Data,varlist=var.list,trt.var="A",
                                    wgt=paste("w",lil,sep=""),beta=betaXY)$wAMD
      # save ATE estimate for this lambda value
      ATE[lil] = ATE_est(fY=Data$Y,fw=Data[,paste("w",lil,sep="")],fA=Data$A)
      
      gns <- cbind(gns, Data[,paste("f.pA",lil,sep="")])
      gns_cv <- cbind(gns_cv, gn_cv)
      
    }, error=function(e){
      gns <- cbind(gns, rep(NA, n)); gns_cv <- cbind(gns_cv, rep(NA, n)); wAMD_vec[lil] <- NA
    })
  } # close loop through lambda values
  
  return(list(gns[,which.min(wAMD_vec)], gns, gns_cv, Q1W = Q1W, Q0W = Q0W))
}

#' Function to compute Shortreed estimators (IPTW + TMLE)
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param family Family for outcome regression for glm
shortreed_est <- function(W, A, Y, family = binomial()){
  gns <- oal(W = W, A = A, Y = Y, family = family)
  g1hat <- gns[[1]]
  iptw_est <- sum(A * Y /(g1hat))/sum(A/g1hat) - sum((1-A)*Y/(1- g1hat))/sum((1-A)/(1-g1hat))
  Q1W <- gns$Q1W; Q0W <- gns$Q0W
  tmle_fit <- tmle::tmle(Y = Y, A = A, W = W, Q = cbind(Q1W, Q0W), g1W = g1hat, gbound = c(1e-4, 1 - 1e-4))
  tmle_est <- tmle_fit$estimates$ATE$psi
  return(list(iptw = iptw_est, tmle = tmle_est))
}

#' Function to do one bootstrap iteration of Shortreed-based estimators
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param family Family for outcome regression for glm
one_shortreed_boot <- function(W, A, Y, family = binomial()){
  n <- length(Y)
  idx <- sample(1:n, replace = TRUE)
  Yij_vec <- sapply(1:n, function(x,idx){
    sum(idx == x)
  }, idx = idx)
  Wstar <- W[idx, , drop = FALSE]
  Astar <- A[idx]
  Ystar <- Y[idx]
  t_star <- shortreed_est(W = Wstar, A = Astar, Y = Ystar, family = family)
  return(list(Yij = Yij_vec, t_star = t_star))
}

#' Function to do bootstrap confidence interval for Shortreed estimators
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param family Family for outcome regression for glm
shortreed_boot <- function(W, A, Y, nboot = 5e2, family = binomial()){
  rslt <- replicate(nboot, one_shortreed_boot(W = W, A = A, Y = Y, family = family))
  all_Yijs <- Reduce(rbind, rslt[1,])
  all_tstars <- apply(Reduce(rbind, rslt[2,]), 2, unlist, use.names= FALSE)
  all_covs <- apply(all_tstars, 2, function(tstar){
    apply(all_Yijs, 2, function(x){
      cov(x, tstar, use = "complete.obs")
    })
  })
  boot_sds <- apply(all_covs, 2, function(x){ sqrt(sum(x^2)) })
  boot_ests <- colMeans(all_tstars, na.rm = TRUE)
  cis <- rbind(boot_ests - 1.96*boot_sds, boot_ests + 1.96*boot_sds)
  return(list(iptw = cis[,1], aiptw = cis[,2], tmle = cis[,3]))
}

#' Function to compute G-computation estimator
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param formula Regression formula for outcome regression
#' @param family Family for outcome regression for glm
gcomp <- function(W, A, Y, formula = "Y ~ .", family = binomail()){
  out_fit <- glm(as.formula(formula), data = data.frame(W, A = A, Y), family = family)
  Q1hat <- predict(out_fit, newdata = data.frame(W, A = 1, Y), type = "response")
  Q0hat <- predict(out_fit, newdata = data.frame(W, A = 0, Y), type = "response")
  return(mean(Q1hat - Q0hat))
}

#' Function to perform one bootstrap resample of G-computation estimator
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param formula Regression formula for outcome regression
#' @param family Family for outcome regression for glm
gcomp_one_boot <- function(W, A, Y, formula = "Y ~ .", family = binomial()){
  resamp_idx <- sample(1:length(Y), replace = TRUE)
  return(gcomp(W[resamp_idx, , drop = FALSE], A = A[resamp_idx], Y = Y[resamp_idx],
               formula = formula, family = family))
}

#' Function to perform bootstrap confidence interval resample for 
#' G-computation estimator
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param formula Regression formula for outcome regression
#' @param family Family for outcome regression for glm
gcomp_boot_ci <- function(W, A, Y, formula = "Y ~ .", nboot = 5e2, family = binomial()){
  ate_star <- replicate(nboot, gcomp_one_boot(W, A, Y, formula, family = family))
  return(as.numeric(quantile(ate_star, p = c(0.025, 0.975))))
}

#' Function to compute IPTW estimator
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param formula Regression formula for outcome regression
#' @param family Family for outcome regression for glm
iptw <- function(W, A, Y, formula = "A ~ .", family = binomial()){
  ps_fit <- glm(as.formula(formula), data = data.frame(W, A = A), family = family)
  g1hat <- predict(ps_fit, type = "response")
  est <- sum(A * Y /(g1hat))/sum(A/g1hat) - sum((1-A)*Y/(1- g1hat))/sum((1-A)/(1-g1hat))
  return(est)
}

#' Function to perform one bootstrap resample of IPTW estimator
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param formula Regression formula for outcome regression
#' @param family Family for outcome regression for glm
iptw_one_boot <- function(W, A, Y, formula = "A ~ .", family = binomial()){
  resamp_idx <- sample(1:length(A), replace = TRUE)
  return(iptw(W[resamp_idx, , drop = FALSE], A = A[resamp_idx], Y = Y[resamp_idx],
              formula = formula, family = family))
}

#' Function to perform bootstrap confidence interval resample for 
#' IPTW estimator
#' @param W Covariates
#' @param A Treatment
#' @param Y Outcome
#' @param formula Regression formula for outcome regression
#' @param family Family for outcome regression for glm
iptw_boot_ci <- function(W, A, Y, formula = "Y ~ .", nboot = 5e2, family = binomail()){
  ate_star <- replicate(nboot, iptw_one_boot(W, A, Y, formula, family = family))
  return(as.numeric(quantile(ate_star, p = c(0.025, 0.975))))
}

#' Help function to check whether true value in a CI
#' @param truth The true value
#' @param ci A two-length vector CI
truth_in_ci <- function(truth, ci){
  truth > min(ci) & truth < max(ci)
}
