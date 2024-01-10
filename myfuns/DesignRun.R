# Main <- function(r, myseed, n, m, tau, set, designs, mods, bootT, is_save, out_path, R) { # is_ohal_cv
# 
#   # print(r)
#   cat(r, "\n")
#   if (r %% 10 == 0) {
#     message(r)
#   }
#   # set.seed(r*100)
#   set.seed(r*myseed, "L'Ecuyer-CMRG")  # Set multicore compatible seed
#   mydat <- GetSimData(n=n, m=m, tau=tau, set=set)
# 
#   #### 4 designs (different covariates are used)
#   for (design in designs) {
#     cat("Design", design, ":\n")
#     DesignRun(mydat, design, mods, bootT, is_save, out_path, R)
#   }
# 
# }



DesignRun <- function (r, mydat, design, mods, bootT, is_homo, 
                       is_save, boot_path, out_path, R, is_boot) {
  
  ## design A (glm only)
  # True covariates (W) and correctly specified models (glm)
  # PS model: logistic regression of A on W1,W2,W3,W7,W8,W9.
  # OR model: linear regression of Y on W1,W2,W3,W7,W8,W9, stratified by A.
  if (design == "A") {
    useCovar <- "W"
    or0.vars <- or1.vars <- c("W1", "W2", "W3", "W7", "W8", "W9")
    mods <- "glm"  # update - glm only
  }
  
  ## design B
  # Covariates (W) are used
  # First apply variable selection, then try different PS and OR models.
  if (design == "B") {
    useCovar <- "W"
  }
  
  ## design C
  # X that are corresponding to true W are used, and try different PS and OR models.
  # PS models: A on X1,X2,X3,X7,X8,X9.
  # OR models: Y on X1,X2,X3,X7,X8,X9, stratified by A.
  if (design == "C") {
    useCovar <- "X"
    or0.vars <- or1.vars <- c("X1", "X2", "X3", "X7", "X8", "X9")
  }
  
  ## design D
  # Covariates (X) are used
  # First apply variable selection, then try different PS and OR models.
  if (design == "D") {
    useCovar <- "X"
  }
  
  ##############################################################################
  ############################ common for 4 designs ############################
  ##############################################################################
  
  if (useCovar == "W") {
    notuseCovar <- "X"
  } else if (useCovar == "X") {
    notuseCovar <- "W"
  }
  # drop W/X for memory
  mydat <- mydat[, !(names(mydat) %in% grep(paste0("^",notuseCovar), names(mydat), value=TRUE))]
  
  # data estimate (VarSelect)
  if (design %in% c("B", "D")) { # need varSelect
    print("Variable selection")
    var_select <- VarSelect(mydat, useCovar, is_homo)
    
    if (is_homo) {
      # or0.vars <- or1.vars <- var_select$vars_min0
      or0.vars <- or1.vars <- var_select$vars_1se0
      
    } else {
      # or0.vars <- var_select$vars_min0
      # or1.vars <- var_select$vars_min1
      or0.vars <- var_select$vars_1se0
      or1.vars <- var_select$vars_1se1
    }
  }
  
  ps.vars <- union(or0.vars, or1.vars)
  
  # data estimate (FitPSnOR)
  print("Fitting PS & OR")
  fits <- FitPSnOR(mydat, ps.vars, or0.vars, or1.vars, design, useCovar)
  ps_lst <- fits$myps
  or0_lst <- fits$myor0
  or1_lst <- fits$myor1
  
  # bootstrap data
  bootdats <- NULL
  if (is_boot) {
    print("Begin boot")
    bootdats <- replicate(bootT, 
                          GetBootDataPSOR(mydat, design, useCovar,
                                          ps.vars, or0.vars, or1.vars),
                          simplify=FALSE)  # return in a list
  }
  
  # save tmp boot & mydat
  tmp_boot <- list(mydat=mydat, design=design, is_homo=is_homo,
                   ps.vars=ps.vars, or0.vars=or0.vars, or1.vars=or1.vars,
                   ps_lst=ps_lst, or0_lst=or0_lst, or1_lst=or1_lst,
                   useCovar=useCovar, bootdats=bootdats, bootT=bootT, mods=mods,
                   is_save=is_save, out_path=out_path)
  mytime <- as.numeric(gsub(":","",strsplit(as.character(Sys.time()), " ")[[1]][2]))
  saveRDS(tmp_boot, file=paste0(boot_path, "/boot_r", r, "_", mytime, ".rds"))
  cat("Save boot data!\n")
  
  # ATE methods on multiple models
  ATEmethods(mydat=mydat, ps_lst=ps_lst, or0_lst=or0_lst, or1_lst=or1_lst, 
             useCovar=useCovar, design=design,
             bootdats=bootdats, bootT=bootT, mods=mods, 
             is_save=is_save, out_path=out_path, r=r, 
             ps.vars=ps.vars, or0.vars=or0.vars, or1.vars=or1.vars)
  
  cat("Finish design", design, "!\n")
  
}


## stop whole program if achieve R replicates
# cnts <- c()
# for (mod in mods) {
#   out_path_full <- paste0(out_path, "_design", design, "_", mod)
#   res_list <- list.files(out_path_full, pattern=".rds", recursive=TRUE)
#   cnts <- c(cnts, length(res_list))
# }
# if (sum(cnts>R) == length(mods)) {
#   cat("Result counts:", cnts, "\n")
#   stop("Reach R replicates! Finished all!")
# }
