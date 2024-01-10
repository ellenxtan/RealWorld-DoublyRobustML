r <- 1 #as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

####tune
set     <- 1 #as.numeric(args[2])     #1,4
is_homo <- FALSE#as.logical(args[3]) #FALSE#TRUE
####

####
mydate  <- "2022-04-05" # "2021-07-30" # Sys.Date()
design  <- "D" #args[1]              #"D"
tau     <- 0
R       <- 1000  #simulated datasets to be used
bootT   <- 500
n       <- 2000
m       <- 9
mods    <- c("glm", "gam", "sl")
is_boot <- TRUE  #as.logical(args[4]) #FALSE#TRUE
pen_tms <- "null"#args[5]             #"pen","tms","null"
is_crossfit <- FALSE#TRUE
####
# sim_sets<- c(4,1)
# taus    <- c(0) #0,1.5,3
# designs <- c("A","B","C","D")
is_sink <- TRUE
is_save <- TRUE
pgmpath <- 'script_ols.R' # program to be executed on each simulated dataset
####

##############################################################################
################################### Setup ####################################
##############################################################################

# load functions
files.sources <- list.files(c("myfuns"), pattern="*.R$", full.names=TRUE)
suppressMessages(invisible(sapply(files.sources, source)))

# main date folder
date_folder <- paste0("results/res_", mydate)
suppressWarnings(dir.create(file.path(date_folder))) 
message(paste0("date_folder: ", date_folder))

# store setting results
set_folder <- paste0(date_folder, "/res_n",n,"_m",m,"_design",design,"_sim",set,"_tau",tau,
                     "_R",R,"_B",bootT,"_homo",is_homo,"_crossfit",is_crossfit)
suppressWarnings(dir.create(file.path(set_folder)))

# store code under each setting
codepath <- paste0(set_folder, "/code")
suppressWarnings(dir.create(file.path(codepath)))
system(paste0("cp ", pgmpath, " ", codepath, "/"))
system(paste0("cp -r myfuns ", codepath, "/"))

# store simulated datasets under each setting
datpath <- paste0(set_folder, "/dats")
suppressWarnings(dir.create(file.path(datpath)))

# store logs under each setting
logpath <- paste0(set_folder, "/logs")
suppressWarnings(dir.create(file.path(logpath)))

# store boot tmp
boot_path <- paste0(set_folder, "/bootdats")
suppressWarnings(dir.create(file.path(boot_path)))

# store result for each mod under each setting
out_path <- paste0(set_folder, "/res")
for (mod in mods) {
    if (design == "A" && mod != "glm") {
        next
    }
    suppressWarnings(dir.create(file.path(paste0(out_path, "_", mod))))
}


## design D
# Covariates (X) are used
# First apply variable selection, then try different PS and OR models.
if (design == "D") {
    useCovar <- "X"
}

if (useCovar == "W") {
    notuseCovar <- "X"
} else if (useCovar == "X") {
    notuseCovar <- "W"
}


##############################################################################
############################## Begin simulation ##############################
##############################################################################

sink_file <- paste0(logpath, "/log_design", design, "_set", set, "_homo",is_homo, 
                    '_boot',is_boot,"_",pen_tms,"_r", r, ".txt") #, "_", mytime

sink(sink_file)

cat("task", r, "\n")
cat("design", design, "\n")
cat("set", set, "\n")
cat("is_boot", is_boot, "\n")
cat("bootT", bootT, "\n")
cat("pen_tms", pen_tms, "\n")
cat("homo", is_homo, "\n")
cat("datpath", datpath, "\n")
cat("logpath", logpath, "\n")
cat("boot_path", boot_path, "\n")
cat("out_path", out_path, "\n")

# mydat <- read.csv(paste0(datpath,'/dat_',r,'.csv'))
set.seed(r*100)
mydat <- GetSimData(n, m, tau, set, is_homo)
write.csv(mydat,paste0(datpath, '/dat_',r,'.csv'), row.names=FALSE)

# drop W/X for memory
mydat <- mydat[, !(names(mydat) %in% grep(paste0("^",notuseCovar), names(mydat), value=TRUE))]

# ATE
tau <- mean(mydat$Y1-mydat$Y0)

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
                                          ps.vars, or0.vars, or1.vars,
                                          is_fit=TRUE),  # use default CI or boot CI
                          simplify=FALSE)  # return in a list
    # save tmp boot & mydat
    tmp_boot <- list(mydat=mydat, design=design, is_homo=is_homo, tau=tau,
                     ps.vars=ps.vars, or0.vars=or0.vars, or1.vars=or1.vars,
                     ps_lst=ps_lst, or0_lst=or0_lst, or1_lst=or1_lst,
                     useCovar=useCovar, mods=mods, pen_tms=pen_tms,
                     bootdats=bootdats, bootT=bootT, 
                     is_save=is_save, out_path=out_path, is_fit=TRUE)  # use default CI or boot CI
    saveRDS(tmp_boot, file=paste0(boot_path, "/boot_r", r, ".rds")) #, "_", mytime
    cat("Save boot data!\n")
}

# estimate ATE with multiple methods
if (!is_crossfit) {
    ATEmethods(mydat=mydat, ps_lst=ps_lst, or0_lst=or0_lst, or1_lst=or1_lst, 
               useCovar=useCovar, design=design,
               bootdats=bootdats, bootT=bootT, mods=mods, 
               is_save=is_save, out_path=out_path, r=r, 
               ps.vars=ps.vars, or0.vars=or0.vars, or1.vars=or1.vars, pen_tms=pen_tms)
    
} else { # cross-fitting for AIPTW & TMLE
    num_fold <- 10
    res_cf <- CrossFit(mydat, ps.vars, or0.vars, or1.vars, design, useCovar, mods, r, num_fold)
    dat_kfold <- res_cf$dat_kfold
    fit_kfold  <- res_cf$fit_kfold
    saveRDS(res_cf, file=paste0(boot_path, "/crossfit_r", r, ".rds")) #, "_", mytime
    cat("Save crossfit data!\n")
    
    for (mod in mods) {
        cat("model:", mod, "\n")
        
        aip_est_lst <- c()
        tml_est_lst <- c()
        aip_psi_lst <- c()
        tml_psi_lst <- c()
        for (k in 1:num_fold) {
            dat <- dat_kfold[[k]]
            fit <- fit_kfold[[k]]
            
            ps <- fit$myps[[mod]]
            or0 <- fit$myor0[[mod]]
            or1 <- fit$myor1[[mod]]
            
            
            aip_k <- AIPTWwrap(dat, ps, or0, or1, bootdats=NULL, mod, is_print=FALSE)
            tml_k <- TMLEwrap(dat, ps, or0, or1, useCovar, bootdats=NULL, mod, is_print=FALSE)
            
            aip_est_lst <- c(aip_est_lst, aip_k)
            tml_est_lst <- c(tml_est_lst, tml_k)
            
            ## 4.2.3. Interactive regression model (IRM)
            # score='ATE' implements the score function
            # https://docs.doubleml.org/stable/guide/scores.html
            # https://towardsdatascience.com/double-machine-learning-for-causal-inference-78e0c6111f9d
            psi <- (or1-or0) + (dat$A*(dat$Y - or1) / ps) - ((1-dat$A)*(dat$Y-or0)/(1-ps))
            aip_psi <- sum((psi - aip_k)^2)
            tml_psi <- sum((psi - tml_k)^2)
            
            aip_psi_lst <- c(aip_psi_lst, aip_psi)
            tml_psi_lst <- c(tml_psi_lst, tml_psi)
        }
        
        aip_cf <- mean(aip_est_lst, na.rm=TRUE)
        tml_cf <- mean(tml_est_lst, na.rm=TRUE)
        aip_sigma2 <- mean(aip_psi_lst, na.rm=TRUE)
        tml_sigma2 <- mean(tml_psi_lst, na.rm=TRUE)
        aip_se <- sqrt(aip_sigma2) / sqrt(nrow(mydat))
        tml_se <- sqrt(tml_sigma2) / sqrt(nrow(mydat))
        aip_ci <- c(aip_cf - 1.96*aip_se, aip_cf + 1.96*aip_se)
        tml_ci <- c(tml_cf - 1.96*tml_se, tml_cf + 1.96*tml_se)
        
        aip_cf_res <- c(aip_cf, aip_se, aip_ci)
        print(c("AIPTW-CF", round(aip_cf_res,4)))
        tml_cf_res <- c(tml_cf, tml_se, tml_ci)
        print(c("TMLE-CF", round(tml_cf_res,4)))
        
        #save
        out_path_full <- paste0(out_path, "_", mod, "/aipcf")
        suppressWarnings(dir.create(file.path(out_path_full)))
        saveRDS(aip_cf_res, file=paste0(out_path_full, "/aipcf_r",r, ".rds"))
        
        out_path_full <- paste0(out_path, "_", mod, "/tmlcf")
        suppressWarnings(dir.create(file.path(out_path_full)))
        saveRDS(tml_cf_res, file=paste0(out_path_full, "/tmlcf_r",r, ".rds"))
        
    }#end mods
    cat("Finish cross-fitting!\n")
}


##############################################################################
############################### End simulation ###############################
##############################################################################

cat("Finish r", r, "!\n")

### end
warnings()
print("Finish one!")

sink()
