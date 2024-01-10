# wrapper for all ATE methods
ATEmethods <- function(mydat, ps_lst, or0_lst, or1_lst, useCovar, design,
                       bootdats, bootT, mods, is_save, out_path, r, 
                       ps.vars, or0.vars, or1.vars, pen_tms, real=FALSE) {
  
  if (!is.null(bootdats)) {
    for (mod in mods) {
      cat("model:", mod, "\n")

      # data estimate
      ps <- ps_lst[[mod]]
      or0 <- or0_lst[[mod]]
      or1 <- or1_lst[[mod]]
      
      # ATE estimation
      print(c("method", "est", "se", "lb", "ub"))

      grp <- GroupDiffwrap(mydat, bootdats)
      out_path_full <- paste0(out_path, "_", mod, "/grp")
      suppressWarnings(dir.create(file.path(out_path_full)))
      saveRDS(grp, file=paste0(out_path_full, "/grp_r",r, ".rds"))

      imp <- IMPwrap(or0, or1, bootdats, mod)
      out_path_full <- paste0(out_path, "_", mod, "/imp")
      suppressWarnings(dir.create(file.path(out_path_full)))
      saveRDS(imp, file=paste0(out_path_full, "/imp_r",r, ".rds"))

      ipw <- IPTWwrap(mydat, ps, bootdats, mod)
      out_path_full <- paste0(out_path, "_", mod, "/ipw")
      suppressWarnings(dir.create(file.path(out_path_full)))
      saveRDS(ipw, file=paste0(out_path_full, "/ipw_r",r, ".rds"))

      aip <- AIPTWwrap(mydat, ps, or0, or1, bootdats, mod)
      out_path_full <- paste0(out_path, "_", mod, "/aip")
      suppressWarnings(dir.create(file.path(out_path_full)))
      saveRDS(aip, file=paste0(out_path_full, "/aip_r",r, ".rds"))

      tml <- TMLEwrap(mydat, ps, or0, or1, useCovar, bootdats, mod)
      out_path_full <- paste0(out_path, "_", mod, "/tml")
      suppressWarnings(dir.create(file.path(out_path_full)))
      saveRDS(tml, file=paste0(out_path_full, "/tml_r",r, ".rds"))

      tryCatch ( {
        dsm <- DSMwrap(mydat, ps, or0, or1, useCovar, bootT, real)  # boot within
        out_path_full <- paste0(out_path, "_", mod, "/dsm")
        suppressWarnings(dir.create(file.path(out_path_full)))
        saveRDS(dsm, file=paste0(out_path_full, "/dsm_r",r, ".rds"))
      }, error=function(e) { print("DSM has error!") } )
      
    } # end mod

    print("Finish all mods w/ bootdats")
  } # end !is.null(bootdats)
  

  ##############################################################################
  
  
  # 2. TMLESL & AIPTWSL w/ boot (use all covars)
  if (is.null(bootdats) && pen_tms == "tms") {
    tic()
    
    covars_all <- grep(paste0("^",useCovar), names(mydat), value=TRUE)
    ps.vars <- covars_all
    or0.vars <- covars_all
    or1.vars <- covars_all
    
    print("Fitting PS & OR using covars_all")
    fits <- FitPSnOR(mydat, ps.vars, or0.vars, or1.vars, design, useCovar)
    ps_lst <- fits$myps
    or0_lst <- fits$myor0
    or1_lst <- fits$myor1
    
    print("Begin boot using covars_all")
    bootdats <- replicate(bootT, 
                          GetBootDataPSOR(mydat, design, useCovar,
                                          ps.vars, or0.vars, or1.vars,
                                          is_fit=TRUE),
                          simplify=FALSE)  # return in a list
    tmp_boot <- list(mydat=mydat, design=design, 
                     ps.vars=ps.vars, or0.vars=or0.vars, or1.vars=or1.vars,
                     useCovar=useCovar, mods=mods, pen_tms=pen_tms,
                     bootdats=bootdats, bootT=bootT, 
                     is_fit=TRUE)
    out_path_boot <- paste0(out_path, "_", "sl", "/bootdats_covars_all")
    suppressWarnings(dir.create(file.path(out_path_boot)))
    saveRDS(tmp_boot, file=paste0(out_path_boot, "/boot_covars_all_r", r, ".rds")) #, "_", mytime
    cat("Save boot data for covars_all!\n")
    
    print("TMLESL&ALPTWSL begins")
    for (mod in mods) {
      cat("model:", mod, "\n")
      
      # data estimate
      ps <- ps_lst[[mod]]
      or0 <- or0_lst[[mod]]
      or1 <- or1_lst[[mod]]
      
      # ATE estimation
      print(c("method", "est", "se", "lb", "ub"))
      
      asl <- AIPTWwrap(mydat, ps, or0, or1, bootdats, mod)
      out_path_full <- paste0(out_path, "_", mod, "/asl")
      suppressWarnings(dir.create(file.path(out_path_full)))
      saveRDS(asl, file=paste0(out_path_full, "/asl_r",r, ".rds"))
      
      tsl <- TMLEwrap(mydat, ps, or0, or1, useCovar, bootdats, mod)
      out_path_full <- paste0(out_path, "_", mod, "/tsl")
      suppressWarnings(dir.create(file.path(out_path_full)))
      saveRDS(tsl, file=paste0(out_path_full, "/tsl_r",r, ".rds"))
      
    }
    print("TMLESL&AIPTWSL ends")
    toc()
  }
  
  ##############################################################################
  
  # 3. PENCOMP only run once
  if (is.null(bootdats) && pen_tms == "pen") {
    tic()
    print("PENCOMP begins")  # boot within
    pen <- PENCOMPwrap(mydat, useCovar, bootT, mods, design, ps.vars, or0.vars, or1.vars, real)
    for (mod in mods) {
      tmp <- pen[[mod]]
      out_path_full <- paste0(out_path, "_", mod, "/pen")
      suppressWarnings(dir.create(file.path(out_path_full)))
      saveRDS(tmp, file=paste0(out_path_full,"/pen_r",r, ".rds"))
    }
    print("PENCOMP ends")
    toc()
  }
  

}


# hal <- OHALwrap(mydat, num_cores, is_ohal_cv, x_vars=ps.vars)
# }

# res <- list(grp=grp, slr=slr, ipw=ipw, aip=aip, tml=tml, tms=tms, dsm=dsm) #, pen=pen))  #, hal=hal
# if (is_save) {
#   saveRDS(res, file=paste0(out_path_full, "/all_r",r,"_", mytime, ".rds"))
# }
