GetResTab <- function(res, tau) {
  if (is.null(res[[1]])) {  # design 1 could be NULL for mod!="glm"
    return(NULL)
  }
  grp_out <- LstToDF(res, "grp")
  slr_out <- LstToDF(res, "slr")
  ipw_out <- LstToDF(res, "ipw")
  aip_out <- LstToDF(res, "aip")
  tml_out <- LstToDF(res, "tml")
  tms_out <- LstToDF(res, "tms")
  dsm_out <- LstToDF(res, "dsm")
  # pen_out <- LstToDF(res, "pen")
  # hal_out <- LstToDF(res, "hal")
  
  tab <- c()
  tab <- rbind(tab, EvaMetrics(tau, grp_out))
  tab <- rbind(tab, EvaMetrics(tau, slr_out))
  tab <- rbind(tab, EvaMetrics(tau, ipw_out))
  tab <- rbind(tab, EvaMetrics(tau, aip_out))
  tab <- rbind(tab, EvaMetrics(tau, tml_out))
  tab <- rbind(tab, EvaMetrics(tau, tms_out))
  tab <- rbind(tab, EvaMetrics(tau, dsm_out))
  # tab <- rbind(tab, EvaMetrics(tau, pen_out))
  # tab <- rbind(tab, EvaMetrics(tau, hal_out))
  
  tab <- as.data.frame(tab)
  names(tab) <- c("bias.abs", "bias.rel", "mse", "coverage", "width", "var.ratio", "power", "power.lb", "power.ub")
  rownames(tab) <- c("GrpDiff", "SLR", "IPTW", "AIPTW", "TMLE", "TMLESL", "DSM")#, "PENCOMP")#, "OHAL")
  print(tab, digits=4)
  # message(tab)
  
  return(tab)
}
