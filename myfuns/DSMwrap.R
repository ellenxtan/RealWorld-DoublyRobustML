DSMwrap <- function(mydat, ps, or0, or1, useCovar, bootT, real=FALSE) { #, ps.vars, pg.vars
  ## fix NA estimate issue in boot in dsmatch package
  
  if (real) {
    covars_all <- names(mydat)[which(!names(mydat) %in% c("A","Y"))]
  } else { # simulation
    covars_all <- grep(paste0("^",useCovar), names(mydat), value=TRUE)
  }
  
  X <- mydat[, covars_all]
  
  res <- dsmatchATE(mydat$Y, X, mydat$A, method="dsm", 
                    ps=ps, pg=cbind(or0, or1), varest=T, boots=bootT)
  dsm_est <- res$est.ds
  dsm_se <- sqrt(res$bootvar)
  
  # Wald interval
  # dsm_lb <- res$est.ds + qnorm(0.025) * dsm_se
  # dsm_ub <- res$est.ds - qnorm(0.025) * dsm_se
  # quantile interval
  dsm_lb <- as.numeric(res$bootq1)
  dsm_ub <- as.numeric(res$bootq2)
  
  res <- c(dsm_est, dsm_se, dsm_lb, dsm_ub)
  print(c("DSM", round(res,4)))
  return(res)
}
