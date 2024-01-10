FitPSnOR <- function(mydat, ps.vars, or0.vars, or1.vars, design, useCovar, 
                     new.ps.dat=NULL, new.or0.dat=NULL, new.or1.dat=NULL,
                     real=FALSE) {
  
  cat("==================\n")
  cat("design", design, "\n")
  print(Sys.time())
  
  if (real) {
    covars_all <- names(mydat)[which(!names(mydat) %in% c("A","Y"))]
  } else { # simulation
    covars_all <- grep(paste0("^",useCovar), names(mydat), value=TRUE)
  }
  # print(covars_all)
  
  if (length(unique(mydat$Y)) == 2) {  # binary Y
    sim_fam <- binomial()
  } else {  # continuous Y
    sim_fam <- gaussian()
  }
  
  myps <- myor0 <- myor1 <- NULL
  ps.fit <- or.fit <- NULL
  
  SL.glmnet2 <- function(...) {
    SL.glmnet(..., useMin = FALSE)
  }
  
  ## default SL libraries in tmle pacakge
  SL.library.ps <- c("SL.glm", "my.dbarts.k.5", "SL.gam")
  SL.library.or <- c("SL.glm", "my.dbarts2", "SL.glmnet", "SL.gam")
  
  ## For maximum accuracy one might try at least the following models: 
  ## glmnet, randomForest, XGBoost, SVM, and bartMachine
  # SL.library.ps <- c("SL.xgboost", "SL.ranger", "SL.gam", "SL.glm", #, "SL.mgcv"
  #                    "my.dbarts.k.5", "SL.glmnet2", "SL.glm.interaction", "SL.glmnet")
  # SL.library.or <- c("SL.xgboost", "SL.ranger", "SL.gam", "SL.glm", #, "SL.mgcv"
  #                    "my.dbarts2", "SL.glmnet2", "SL.glm.interaction", "SL.glmnet")
  
  # data for glm & gam (after variable selection)
  use.ps.dat <- mydat[, ps.vars, drop=FALSE]  # ensure still df
  if (is.null(new.ps.dat)) {
    new.ps.dat <- use.ps.dat
  }
  use.or0.dat <- mydat[, c("A", or0.vars), drop=FALSE]
  if (is.null(new.or0.dat)) {
    new.or0.dat <- data.frame(A=0, mydat[, or0.vars, drop=FALSE])
  }
  use.or1.dat <- mydat[, c("A", or1.vars), drop=FALSE]
  if (is.null(new.or1.dat)) {
    new.or1.dat <- data.frame(A=1, mydat[, or1.vars, drop=FALSE])
  }
  # data for SL (no variable selection beforehead)
  # new.or0.dat.SL <- data.frame(A=0, mydat[, covars_all])
  # new.ps.dat.SL <- mydat[, covars_all]
  # new.or1.dat.SL <- data.frame(A=1, mydat[, covars_all])
  
  # prespecified CV folds
  num_fold <- 5 #2 #
  validRows=list()
  for (v in 1:5) { #5 #2
    validRows[[v]] <- which(mydat$fold_idx==v)
  }
  
  #### PS: A~W
  if (!is.null(ps.vars)) {
    tic()
    
    ps.fit <- suppressWarnings(SuperLearner(Y=mydat$A, X=use.ps.dat, newX=new.ps.dat,
                           family=binomial(), SL.library=SL.library.ps,
                           cvControl = list(V=num_fold, shuffle=FALSE, validRows=validRows)
                           ))
    pred <- predict(ps.fit)
    
    myps <- list()
    myps[["glm"]] <- as.vector(pred$library.predict[, which(colnames(pred$library.predict)=="SL.glm_All")])
    if (design != "A") {  # other models
      myps[["gam"]] <- as.vector(pred$library.predict[, which(colnames(pred$library.predict)=="SL.gam_All")])
      myps[["sl"]] <- as.vector(pred$pred)
    }
    
    print(ps.fit)
    toc()
  }
  
  #### OR: Y~A+W
  if (!is.null(or0.vars) && !is.null(or1.vars)) {
    tic()
    
    ## or.fit0
    or.fit <- suppressWarnings(SuperLearner(Y=mydat$Y, X=use.or0.dat, newX=new.or0.dat,
                           family=sim_fam, SL.library=SL.library.or,
                           cvControl = list(V=num_fold, shuffle=FALSE, validRows=validRows)
                           ))
    pred <- predict(or.fit)
    
    myor0 <- list()
    myor0[["glm"]] <- as.vector(pred$library.predict[, which(colnames(pred$library.predict)=="SL.glm_All")])
    if (design != "A") {  # other models
      myor0[["gam"]] <- as.vector(pred$library.predict[, which(colnames(pred$library.predict)=="SL.gam_All")])
      myor0[["sl"]] <- as.vector(pred$pred)
    }
    
    print(or.fit)
    
    ## or.fit1
    or.fit <- suppressWarnings(SuperLearner(Y=mydat$Y, X=use.or1.dat, newX=new.or1.dat,
                           family=sim_fam, SL.library=SL.library.or,
                           cvControl = list(V=num_fold, shuffle=FALSE, validRows=validRows)
                           ))
    pred <- predict(or.fit)
    
    myor1 <- list()
    myor1[["glm"]] <- as.vector(pred$library.predict[, which(colnames(pred$library.predict)=="SL.glm_All")])
    if (design != "A") {  # other models
      myor1[["gam"]] <- as.vector(pred$library.predict[, which(colnames(pred$library.predict)=="SL.gam_All")])
      myor1[["sl"]] <- as.vector(pred$pred)
    }
    
    print(or.fit)
    toc()
  }
  
  
  return(list(myps=myps, myor1=myor1, myor0=myor0))
  #, ps.fit=ps.fit, or.fit=or.fit))  # not return for memory
}
