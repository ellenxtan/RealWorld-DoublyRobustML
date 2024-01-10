#dataOR=simdat; propenVarList=propenVarList; outcomeVarList0=outcomeVarList0; 
#outcomeVarList1=outcomeVarList1;
#treat.varname=treat.varname; outcome.varname=outcome.varname;
#original=1
### in this paper, we used the simplest basis functions-linear bases
###dataOR=original dataset
###fit models on bootstrap samples
###propenVarList-propensity model
###outcomeVarList0-model for Y0
###outcomeVarList1-model for Y1
###numKnots-number of knots truncated linear basis functions

PENCOMPwrap <- function(mydat, useCovar, bootT, mods, design, 
                        ps.vars, or0.vars, or1.vars, real=FALSE) {
  
  mulResult_all=NULL
  mulResult_all=replicate(bootT, 
                          pencomp(dataOR=mydat, mods=mods, 
                                  useCovar=useCovar, design=design, 
                                  ps.vars=ps.vars, 
                                  or0.vars=or0.vars, or1.vars=or1.vars, real=real) )
  
  print(class(mulResult_all))
  res = processPENCOMP(mulResult_all, mods)
  
  # print(c("PENCOMP", res))
  return(res)
}



#####Zhou&Little2019JASA
pencomp=function(dataOR, mods, useCovar, design, 
                 ps.vars, or0.vars, or1.vars, real,
                 treat.varname="A", outcome.varname="Y", original=0, numKnot=35) {
  
  or.mod <- "glm"
  
  # if (real) {
  #   or.mod <- "glm"
  # } else {
  #   if (useCovar == "W") {
  #     or.mod <- "glm"
  #   } else if (useCovar == "X") {
  #     or.mod <- "gam"
  #   }
  # }

  tryCatch (
    {
      data=NULL
      if (original==1){  # data estimate
        data=dataOR
      } else {  # bootstrap
        
        treatID=which(dataOR[, treat.varname]==1)
        controlID=which(dataOR[, treat.varname]==0)
        data = dataOR[c( sample(treatID,replace=T), sample(controlID,replace=T) ),]  
        ##in the simulation studies, we did stratified bootstraps, results were similar with random bootstraps 
        #data = dataOR[sample(1:dim(dataInput)[1],replace=T),]  ###random bootstraps  ##
      }
      
      treatInd = dataOR[, treat.varname]  ###indicator of treated subjects in the original dataset
      treatBoot = data[, treat.varname]  ###indicator of treated subjects in the bootstrap sample
      Yobs = dataOR[, outcome.varname]
      
      ######## propensity score model ########################################################### 
      # propen.model=NULL
      # propen.model=formulaF(varList=propenVarList, y.name=treat.varname)
      # model2a=glm(propen.model, data=data, family="binomial", control = list(maxit = 50))
      # prob.boot=predict(model2a, newdata=data, type="response")
      
      # var_select <- VarSelect(data, useCovar, real)
      # ps.vars <- var_select$ps.vars
      # or.vars <- var_select$or.vars
      
      propenVarList <- ps.vars
      outcomeVarList0 <- or0.vars
      outcomeVarList1 <- or1.vars
      
      ## fit ps model for all mods (only fit ps model not or) - xtan
      print("Fit PENCOMP PS & OR 1")
      model2a <- FitPSnOR(data, 
                          ps.vars=propenVarList, or0.vars=NULL, or1.vars=NULL, 
                          design=design, useCovar=useCovar, real=real,
                          new.ps.dat=NULL, new.or0.dat=NULL, new.or1.dat=NULL)
      
      newData0=dataOR
      newData0[, treat.varname]=0  # change outcome A
      print("Fit PENCOMP PS & OR 2")
      model0.fit <- FitPSnOR(data, 
                             ps.vars=propenVarList, or0.vars=NULL, or1.vars=NULL, 
                             design=design, useCovar=useCovar, real=real,
                             new.ps.dat=newData0[, propenVarList, drop=FALSE], 
                             new.or0.dat=NULL, new.or1.dat=NULL)
      myps0 <- model0.fit$myps
      
      newData1=dataOR
      newData1[, treat.varname]=1  # change outcome A
      print("Fit PENCOMP PS & OR 3")
      model1.fit <- FitPSnOR(data, 
                             ps.vars=propenVarList, or0.vars=NULL, or1.vars=NULL, 
                             design=design, useCovar=useCovar, real=real,
                             new.ps.dat=newData1[, propenVarList, drop=FALSE], 
                             new.or0.dat=NULL, new.or1.dat=NULL)
      myps1 <- model1.fit$myps
      
      # xtan
      ATEs <- list()
      for (mod in mods) {
        ATEs[[mod]] <- NULL
      }
      
      for (mod in mods) {
        
        ATE <- tryCatch({
          
          ## debug
          # if (mod == "gam") {
          #   stop("test point")
          # }
          
          prob.boot <- model2a$myps[[mod]]  # on boot data
          temp2a=NULL
          temp2a=as.numeric(treatBoot==1)*prob.boot + as.numeric(treatBoot==0)*(1-prob.boot)  ###estimate propensity of observed treatment
          # print(sum(temp2a==1))
          # print(sum(temp2a==0))
          data$pslogit = temp2a#log(temp2a/(1-temp2a))   ####spline on the logit of propensity
          # print(range(data$pslogit))
          
          
          ####fit the prediction model for all the models
          ###for the outcome model separate model y0 and y1
          ##############################################################
          ##############################################################
          ###use equally spaced fixed knots assuming K knots
          data0=data[data[treat.varname]==0,]  # among ctrl
          pspp0=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarList0, 
                           propen.score=data0$pslogit, or.mod=or.mod,  # xtan
                           usedata=data0, num.knot=numKnot)
          
          
          ###imputing the missing potential outcomes (among all subjs) on original data
          # newData0=dataOR
          # newData0[, treat.varname]=0  # change outcome A
          # xtan: SuperLearner form
          # predict2a=1 - as.vector(predict(model2a$ps.fit, newdata=newData0[, propenVarList], onlySL = TRUE)$pred)
          predict2a <- myps0[[mod]]
          
          newData0$pslogit = predict2a#log(predict2a/(1-predict2a))  ####spline on the logit of propensity
          # print(range(newData0$pslogit))
          imputed0 = imputeF(newdata=newData0, model=pspp0, x.varnames=pspp0$x.varnames, 
                             propen.score.new=newData0$pslogit)
          # print(mean(imputed0))
          
          ##########################################################################################################################
          ###prediction model for Y under treatment using treated data
          data1=data[data[treat.varname]==1,]  # among treated
          pspp1=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarList1, 
                           propen.score=data1$pslogit, or.mod=or.mod,  # xtan
                           usedata=data1, num.knot=numKnot)
          
          ###imputing the missing potential outcomes (among all subjs) on original data
          # newData1=dataOR
          # newData1[, treat.varname]=1  # change outcome A
          # xtan: SuperLearner form
          # predict2a=as.vector(predict(model2a$ps.fit, newdata=newData1[, propenVarList], onlySL = TRUE)$pred)
          predict2a <- myps1[[mod]]
          
          newData1$pslogit = predict2a#log(predict2a/(1-predict2a))   ####spline on the logit of propensity
          # print(range(newData1$pslogit))
          imputed1 = imputeF(newdata=newData1, model=pspp1, x.varnames=pspp1$x.varnames, 
                             propen.score.new=newData1$pslogit)
          # print(mean(imputed1))
          
          ##############keep the observed outcome the same#####################
          imputed1[treatInd==1]=Yobs[treatInd==1]
          imputed0[treatInd==0]=Yobs[treatInd==0]
          
          #######include everyone
          diff <- imputed1-imputed0
          nrows = length(diff[!is.na(diff)])  # exclude NA
          ATE=c(mean(diff, na.rm=T), var(diff, na.rm=T)/nrows )
          
          ATE  # return
          
        }, error=function(e) e)
        
        if(inherits(ATE, "error")) next
        
        print(c(mod, ATE))
        ATEs[[mod]] <- ATE
        cat("Finish PENCOMP mod", mod, "\n")
      }
      
      return( list(ATEs=ATEs) ) ###estimate and variance of ATE
      
  }, error=function(e) return(list(ATEs=ATEs)) )
}

# x <- tryCatch({
#   a=1 
#   a
# }, error = function(e) e)


##########################process the bootstrap estimates ##############
processPENCOMP=function(mulResult_all, mods){
  
  res_lst <- list()
  for (mod in mods) {
    res_lst[[mod]] <- NULL
  }
  
  for (mod in mods) {

    pen_res <- tryCatch({
      
      mulResult <- LstToDF(mulResult_all, mod)
      
      if (class(mulResult)[1] == "NULL") {
        stop(NULL)
      }
      
      if (class(mulResult)[1] == "list") {
        # remove NA inner lists
        mulResult <- Filter(function(a) any(!is.na(a)), mulResult)
        mulResult <- as.data.frame(do.call(cbind, mulResult))
      }
      
      x=mulResult[,1] # estimate; xtan
      y=mulResult[,2] # variance; xtan
      x=x[!is.na(x)]
      y=y[!is.na(y)]
      
      numT=length(x)
      theta_d=mean(x)
      Wd=mean(y)  ### within imputation variance
      Bd=(1/(numT-1))*(sum((x-mean(x))^2))###between imputation variance
      Td=Wd+(1+1/numT)*Bd ###total variability associated with mean
      # print(c("Td", Td))
      v=(numT-1)*(1+(1/(1/numT+1))*(Wd/Bd))^2 ##degree of freedom
      
      pen_res <- c(theta_d[1], sqrt(Td[1]), theta_d[1] + c(-1,1)*qt(0.975,v[1])*sqrt(Td[1]))
      
      print(c("PENCOMP", round(pen_res,4), mod))
      
      pen_res  # return
      
    }, error=function(e) e)
    
    if(inherits(pen_res, "error")) next
    
    res_lst[[mod]] <- pen_res
  }
  
  return(res_lst)
      
}


######formula construct
formulaF=function(varList, y.name){
  return ( as.formula(paste(y.name, "~ ", paste(c(varList), collapse = "+"))) )
}


###formula construct for generalized model
formulaGAM=function(varList, y.name, spline){
  return ( as.formula(paste(y.name, " ~ ", paste(c(varList), collapse = "+"), "+", spline)) )
}




#########################################################
### assume a truncated linear basis
### right now it only supports truncated linear basis with equally spaced knot
### will expand on this later
### x.varnames enter as c("L1", "L3"), separate variable names by comma
pencompFit=function(y.varname, x.varnames, propen.score, or.mod, usedata, num.knot) {
  
  # print(sum(is.na(usedata)))

  response=usedata[, y.varname]
  covariateX=NULL
  corr.ps <- c()
  for(i in 1:length(x.varnames)){
    covariateX=cbind(covariateX, usedata[, x.varnames[i] ])
    corr.ps <- c(corr.ps, abs(cor(usedata[, x.varnames[i] ], propen.score, use = "complete.obs")))
  }
  covariateX=cbind(rep(1,nrow(usedata)), covariateX, propen.score)
  
  # drop one covariate if high corr with propen.score (avoid singularity part1)
  if (max(corr.ps)>0.5) {
    max.id <- which.max(corr.ps)
    covariateX <- covariateX[,-(max.id+1)] # add 1 for intercept column
    x.varnames <- x.varnames[-max.id]
  }
  
  space=(max(propen.score)-min(propen.score))/(num.knot+1)
  knots=(min(propen.score)+space*(1:num.knot))
  linearB=NULL
  linearB=outer(propen.score, knots, "-")
  linearB=linearB * (linearB > 0)
  linearB <- ifelse(linearB==0, 0.0000001, linearB)  # (avoid singularity part2)
  
  # mixed model: fixed for linear, random for spline
  all=rep(1, dim(usedata)[1])
  usedata$all <- rep(1, nrow(usedata))  # for random effects
  usedata$intercept <- rep(1, nrow(usedata))
  usedata$propen.score <- propen.score
  
  stopifnot(or.mod %in% c("glm", "gam"))
  
  if (or.mod == "glm") {
    psppM <- nlme::lme(response ~ covariateX-1, random=list(all=nlme::pdIdent(~0+linearB)))
    # Y ~ intercept + W1 + W2 + W3 + propen.score -1
    # or.fml <- as.formula(paste(y.varname, "~ intercept + ", paste(x.varnames, collapse = "+"), "+ propen.score -1"))
    # psppM <- nlme::lme(or.fml, random=list(all=nlme::pdIdent(~0+linearB)), data=usedata)
    fixCoef <- psppM$coefficients$fixed
    names(fixCoef) <- c("intercept", x.varnames, "propen.score")
    randCoef <- psppM$coefficients$random$all
    sigmaRes <- psppM$sigma
    
  } else if (or.mod == "gam") {
    x.fml <- c()
    for (v in 1:length(x.varnames)) {
      if (length(unique(usedata[[x.varnames[v]]]))==2) {  # no smooth for discrete vars
        x.fml <- paste(x.fml, x.varnames[v])
      } else {
        x.fml <- paste(x.fml, "s(", x.varnames[v], ", k=-1)")  # degree - xtan
      }
      if (v != length(x.varnames)) {
        x.fml <- paste(x.fml, "+")
      }
    }  # Y ~ intercept + s(W1) + s(W2) + s(W3) + s(propen.score) -1
    or.fml <- as.formula(paste("Y ~ intercept +", x.fml, "+ propen.score -1"))  # degree - xtan
    # print(or.fml)
    # assign linearB to global envir - o.w. could not find linearB when fit gamm
    assign("linearB", linearB, envir = .GlobalEnv)
    
    ctrl <- nlme::lmeControl(maxIter=50, opt='optim')
    psppM <- mgcv::gamm(or.fml, random=list(all=nlme::pdIdent(~0+linearB)), data=usedata, 
                        control=ctrl)  #method="REML", - REML gives sigular warning; use default ML
    fixCoef <- psppM$lme$coefficients$fixed
    names(fixCoef) <- c("intercept", x.varnames, "propen.score")
    randCoef <- psppM$lme$coefficients$random$all
    sigmaRes <- psppM$lme$sigma
    
    rm(linearB)
    # psppM$coefficients$fixed
    # psppM$coefficients$random$all
    # 
    # psppM1$lme$coefficients$fixed
    # psppM1$gam$coefficients
    # psppM1$lme$coefficients$random$all
    # identical(as.vector(psppM1$lme$coefficients$fixed), as.vector(psppM1$gam$coefficients))
    
  }
  
  return(list(fixed=fixCoef, random=randCoef, knot.loc=knots, 
              sigmaRes=sigmaRes, x.varnames=x.varnames))
}






######################################################################
##########prediction values############################################ 
###assume a truncated linear basis
###model is fitted model from lmeFit
###imputing missing potential outcomes

imputeF=function(newdata, model, x.varnames, propen.score.new) {
  
  knots=model$knot.loc
  
  linearB=NULL
  linearB =outer(propen.score.new, knots, "-")
  linearB =linearB * (linearB > 0)
  
  designM=cbind(rep(1,nrow(newdata)), newdata[, x.varnames], propen.score.new)
  
  designM=as.matrix(designM)
  predicted = designM %*% model$fixed + as.matrix(linearB) %*% as.vector(unlist(model$random)) + rnorm(nrow(newdata), 0, model$sigmaRes)
  
  return(predicted)
  
}



## xtan: 1iter
## FitPSnOR_1iter(refl, ps.vars, or0.vars, or1.vars, design="sim_refl")
pencomp_1iter=function(dataOR, mods, design, #useCovar, 
                 ps_lst, 
                 ps.vars, or0.vars, or1.vars#, #real,
                 # original=1
                 ) {
  
  treat.varname="A"
  outcome.varname="Y"
  numKnot=35
  
  or.mod <- "glm"
  
  # if (real) {
  #   or.mod <- "glm"
  # } else {
  #   if (useCovar == "W") {
  #     or.mod <- "glm"
  #   } else if (useCovar == "X") {
  #     or.mod <- "gam"
  #   }
  # }
  
  tryCatch (
    {
      data=dataOR
      
      treatInd = dataOR[, treat.varname]  ###indicator of treated subjects in the original dataset
      treatBoot = data[, treat.varname]  ###indicator of treated subjects in the bootstrap sample
      Yobs = dataOR[, outcome.varname]
      
      ######## propensity score model ########################################################### 
      # propen.model=NULL
      # propen.model=formulaF(varList=propenVarList, y.name=treat.varname)
      # model2a=glm(propen.model, data=data, family="binomial", control = list(maxit = 50))
      # prob.boot=predict(model2a, newdata=data, type="response")
      
      # var_select <- VarSelect(data, useCovar, real)
      # ps.vars <- var_select$ps.vars
      # or.vars <- var_select$or.vars
      
      propenVarList <- ps.vars
      outcomeVarList0 <- or0.vars
      outcomeVarList1 <- or1.vars
      
      ## fit ps model for all mods (only fit ps model not or) - xtan
      # print("Fit PENCOMP PS & OR 1")
      # model2a <- FitPSnOR_1iter(data, 
      #                     ps.vars=propenVarList, or0.vars=NULL, or1.vars=NULL, 
      #                     design=design, #useCovar=useCovar, real=real,
      #                     new.ps.dat=NULL, new.or0.dat=NULL, new.or1.dat=NULL)
      
      newData0=dataOR
      newData0[, treat.varname]=0  # change outcome A
      print("Fit PENCOMP PS & OR 2")
      model0.fit <- FitPSnOR_1iter(data, 
                             ps.vars=propenVarList, or0.vars=NULL, or1.vars=NULL, 
                             design=design, #useCovar=useCovar, real=real,
                             new.ps.dat=newData0[, propenVarList, drop=FALSE], 
                             new.or0.dat=NULL, new.or1.dat=NULL)
      myps0 <- model0.fit$myps
      
      newData1=dataOR
      newData1[, treat.varname]=1  # change outcome A
      print("Fit PENCOMP PS & OR 3")
      model1.fit <- FitPSnOR_1iter(data, 
                             ps.vars=propenVarList, or0.vars=NULL, or1.vars=NULL, 
                             design=design, #useCovar=useCovar, real=real,
                             new.ps.dat=newData1[, propenVarList, drop=FALSE], 
                             new.or0.dat=NULL, new.or1.dat=NULL)
      myps1 <- model1.fit$myps
      
      # xtan
      ATEs <- list()
      for (mod in mods) {
        ATEs[[mod]] <- NULL
      }
      
      for (mod in mods) {
        
        ATE <- tryCatch({
          
          ## debug
          # if (mod == "gam") {
          #   stop("test point")
          # }
          
          prob.boot <- ps_lst[[mod]]  #model2a$myps[[mod]]  # on boot data
          temp2a=NULL
          temp2a=as.numeric(treatBoot==1)*prob.boot + as.numeric(treatBoot==0)*(1-prob.boot)  ###estimate propensity of observed treatment
          # print(sum(temp2a==1))
          # print(sum(temp2a==0))
          data$pslogit = temp2a#log(temp2a/(1-temp2a))   ####spline on the logit of propensity
          # print(range(data$pslogit))
          
          
          ####fit the prediction model for all the models
          ###for the outcome model separate model y0 and y1
          ##############################################################
          ##############################################################
          ###use equally spaced fixed knots assuming K knots
          data0=data[data[treat.varname]==0,]  # among ctrl
          pspp0=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarList0, 
                           propen.score=data0$pslogit, or.mod=or.mod,  # xtan
                           usedata=data0, num.knot=numKnot)
          
          
          ###imputing the missing potential outcomes (among all subjs) on original data
          # newData0=dataOR
          # newData0[, treat.varname]=0  # change outcome A
          # xtan: SuperLearner form
          # predict2a=1 - as.vector(predict(model2a$ps.fit, newdata=newData0[, propenVarList], onlySL = TRUE)$pred)
          predict2a <- myps0[[mod]]
          
          newData0$pslogit = predict2a#log(predict2a/(1-predict2a))  ####spline on the logit of propensity
          # print(range(newData0$pslogit))
          imputed0 = imputeF(newdata=newData0, model=pspp0, x.varnames=pspp0$x.varnames, 
                             propen.score.new=newData0$pslogit)
          # print(mean(imputed0))
          
          ##########################################################################################################################
          ###prediction model for Y under treatment using treated data
          data1=data[data[treat.varname]==1,]  # among treated
          pspp1=pencompFit(y.varname=outcome.varname, x.varnames=outcomeVarList1, 
                           propen.score=data1$pslogit, or.mod=or.mod,  # xtan
                           usedata=data1, num.knot=numKnot)
          
          ###imputing the missing potential outcomes (among all subjs) on original data
          # newData1=dataOR
          # newData1[, treat.varname]=1  # change outcome A
          # xtan: SuperLearner form
          # predict2a=as.vector(predict(model2a$ps.fit, newdata=newData1[, propenVarList], onlySL = TRUE)$pred)
          predict2a <- myps1[[mod]]
          
          newData1$pslogit = predict2a#log(predict2a/(1-predict2a))   ####spline on the logit of propensity
          # print(range(newData1$pslogit))
          imputed1 = imputeF(newdata=newData1, model=pspp1, x.varnames=pspp1$x.varnames, 
                             propen.score.new=newData1$pslogit)
          # print(mean(imputed1))
          
          ##############keep the observed outcome the same#####################
          imputed1[treatInd==1]=Yobs[treatInd==1]
          imputed0[treatInd==0]=Yobs[treatInd==0]
          
          #######include everyone
          diff <- imputed1-imputed0
          nrows = length(diff[!is.na(diff)])  # exclude NA
          ATE=c(mean(diff, na.rm=T), var(diff, na.rm=T)/nrows )
          
          ATE  # return
          
        }, error=function(e) e)
        
        if(inherits(ATE, "error")) next
        
        print(c(mod, ATE))
        ATEs[[mod]] <- ATE
        cat("Finish PENCOMP mod", mod, "\n")
      }
      
      return( list(ATEs=ATEs) ) ###estimate and variance of ATE
      
    }, error=function(e) return(list(ATEs=ATEs)) )
}





