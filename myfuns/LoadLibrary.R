# .libPaths(c(.Library,'/home/l014856/R/x86_64-pc-linux-gnu-library/4.0'))

##download repo from github and unzip
# devtools::install("dsmatch-master")  # fix NA estimate issue in boot in dsmatch package
# devtools::install("lqa-master")
# devtools::install_github("benkeser/drtmle")

# /lrlhps/apps/R/qualified/R-4.0.3/bin/R -f "LoadLibrary.R" >& "LoadLibrary.R.txt" &
# /lrlhps/apps/R/qualified/R-3.6.0/bin/R -f "LoadLibrary.R" >& "LoadLibrary.R.txt" &

# .libPaths(c(.Library,'/home/l014856/R/x86_64-pc-linux-gnu-library/4.0'))

library(dplyr)

library(survey)

# library(DoubleML)
# library(mlr3)
# library(mlr3learners)
# # remotes::install_github("mlr-org/mlr3extralearners")
# library(mlr3extralearners)
# # a=list_mlr3learners(select = c("id", "mlr3_package", "required_packages"))
# # View(a)
# library(mlr3pipelines)

# ps & or models
library(glmnet)
library(ranger)
# library(polspline)  # "SL.polymars"
library(xgboost)
# library(e1071)  # svm
# library(ck37r)
# library(gbm)
# library(gam)
# options(java.parameters = "-Xmx2500m")
# library(bartMachine)

# dsm
library(dsmatch)#, lib.loc="/home/l014856/R/x86_64-pc-linux-gnu-library/4.0")

# pencomp
library(nlme)
library(mgcv)
library(MASS)

# tmle
library(SuperLearner)#, lib.loc="/home/l014856/R/x86_64-pc-linux-gnu-library/4.0")
library(tmle)#, lib.loc="/home/l014856/R/x86_64-pc-linux-gnu-library/4.0")
library(dbarts)   # need comment out all .libPath & use R-4.0.3;   "tmle.SL.dbarts2"

# tmle-sl
# library(dbarts)

# ohal
# library(gam)
# library(lqa)
# library(drtmle)
# library(hal9001)

library(parallel)#, lib.loc="/home/l014856/R/x86_64-pc-linux-gnu-library/4.0")
library(doParallel)
# library(future.batchtools)
library(tictoc)#, lib.loc="/home/l014856/R/x86_64-pc-linux-gnu-library/4.0")

print("Finished load lib!")
