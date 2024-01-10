#### boot data & models
GetBootDataPSOR <- function(mydat, design, useCovar, ps.vars, or0.vars, or1.vars, is_fit=FALSE, real=FALSE) {

  #### boot data (stratified by trt)
  # resamp_idx <- sample(1:nrow(mydat), replace = TRUE)
  # bootdat <- mydat[resamp_idx, , drop = FALSE]
  contlID <- which(mydat[, "A"]==0)
  treatID <- which(mydat[, "A"]==1)
  dat0 <- mydat[contlID, ]
  dat1 <- mydat[treatID, ]
  # bootdat <- mydat[c( sample(treatID,replace=T), sample(contlID,replace=T) ),]
  
  # bootstrap stratified by both trt and CVfold
  boot0 <- dat0 %>%
    group_by(fold_idx) %>%
    sample_n(size=length(fold_idx), replace=TRUE)
  boot1 <- dat1 %>%
    group_by(fold_idx) %>%
    sample_n(size=length(fold_idx), replace=TRUE)
  bootdat <- as.data.frame(rbind(boot0,boot1))
  
  print(table(boot0$fold_idx))
  print(table(boot1$fold_idx))
  
  #### variable selection
  # if (design %in% c("B", "D", "real")) { # need varSelect
  #   var_select <- VarSelect(bootdat, useCovar, real)
  #   if (length(var_select$ps.vars) != 0) {
  #     ps.vars <- var_select$ps.vars
  #     or.vars <- var_select$or.vars
  #   }
  # }
  
  #### PS & OR
  if (is_fit) {
    fits <- FitPSnOR(bootdat, ps.vars, or0.vars, or1.vars, design, useCovar, real=real)
    ps_lst <- fits$myps
    or0_lst <- fits$myor0
    or1_lst <- fits$myor1
  } else {
    ps_lst <- NULL
    or0_lst <- NULL
    or1_lst <- NULL
  }
  
  return(list(mydat=bootdat, ps_lst=ps_lst, or0_lst=or0_lst, or1_lst=or1_lst))
  #, ps.fit=ps.fit, or.fit=or.fit))
}

