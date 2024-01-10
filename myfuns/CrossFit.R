CrossFit <- function(mydat, ps.vars, or0.vars, or1.vars, design, useCovar, mods, r,
                     num_fold=5) {
  yourData <- mydat
  
  #Randomly shuffle the data
  yourData<-yourData[sample(nrow(yourData)),]
  
  #Create k equally size folds
  folds <- cut(seq(1,nrow(yourData)),breaks=num_fold,labels=FALSE)
  
  #Perform k fold cross validation
  dat_kfold <- list()
  fit_kfold <- list()
  for(k in 1:num_fold){
    #Segement your data by fold using the which() function 
    testIndexes <- which(folds==k, arr.ind=TRUE)
    testData <- yourData[testIndexes, ]
    trainData <- yourData[-testIndexes, ]
    #Use the test and train data partitions however you desire...
    
    new.ps.dat <- testData[, ps.vars, drop=FALSE]
    new.or0.dat <- data.frame(A=0, testData[, or0.vars, drop=FALSE])
    new.or1.dat <- data.frame(A=1, testData[, or1.vars, drop=FALSE])

    fits <- FitPSnOR(trainData, ps.vars, or0.vars, or1.vars, design, useCovar,
                     new.ps.dat=new.ps.dat,
                     new.or0.dat=new.or0.dat, new.or1.dat=new.or1.dat)
    
    #save
    dat_kfold[[k]] <- testData
    fit_kfold[[k]] <- fits
  }
  
  return(list(dat_kfold=dat_kfold, fit_kfold=fit_kfold))
}


# library(caret)
# 
# data(iris)
# flds <- createFolds(iris$Sepal.Length, k = 5, list = TRUE, returnTrain = TRUE)
# # names(flds)[1] <- "train"
# 
# for (k in 1:5) {
#   print(nrow(iris[ flds[[k]], ]))
# }
