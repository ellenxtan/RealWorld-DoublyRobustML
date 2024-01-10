## tmle.SL.dbarts2 is in the default library for estimating Q. It uses the default setting in the dbarts package, k=2. tmle.SL.dbarts.k.5 is used to estimate the components of g. It sets k=0.5, to avoid shrinking predicted values too far from (0,1). See bart documentation for more information.
# https://rdrr.io/cran/tmle/man/tmle.Sl.dbarts2.html

## add prediction based on:
# https://github.com/ecpolley/SuperLearner/issues/102

## example - not run
# data(Boston, package = "MASS")
# Y = Boston$medv
# # Remove outcome from covariate dataframe.
# X = Boston[, -14]
# set.seed(1)
# # Sample rows to speed up example.
# row_subset = sample(nrow(X), 30)
# 
# library(dbarts)
# library(SuperLearner)
# sl = SuperLearner(Y[row_subset], X[row_subset, ], newX=X[1, ], family = gaussian(),
#                   cvControl = list(V = 2), SL.library = c("SL.mean", "my.dbarts.k.5", "my.dbarts2"))
# print(sl)
# predict(sl)

my.dbarts.k.5 = function(Y, X, newX, family, obsWeights, id,
                         sigest = NA,
                         sigdf = 3,
                         sigquant = 0.90,
                         k = 0.5,  # xtan
                         power = 2.0,
                         base = 0.95,
                         binaryOffset = 0.0,
                         ntree = 200,
                         ndpost = 1000,
                         nskip = 100,
                         printevery = 100,
                         keepevery = 1,
                         keeptrainfits = T,
                         usequants = F,
                         numcut = 100,
                         printcutoffs = 0,
                         nthread = 1,
                         keepcall = T,
                         verbose = F,
                         ...) {
  
  # .SL.require("dbarts")
  
  model =
    dbarts::bart(x.train = X,
                 y.train = Y,
                 # We need to pass newX in directly due to lack of prediction.
                 x.test = newX,
                 sigest = sigest,
                 sigdf = sigdf,
                 sigquant = sigquant,
                 k = k,
                 power = power,
                 base = base,
                 binaryOffset = binaryOffset,
                 weights = obsWeights,
                 ntree = ntree,
                 ndpost = ndpost,
                 nskip = nskip,
                 printevery = printevery,
                 keepevery = keepevery,
                 keeptrainfits = keeptrainfits,
                 usequants = usequants,
                 numcut = numcut,
                 printcutoffs = printcutoffs,
                 nthread = nthread,
                 keepcall = keepcall,
                 verbose = verbose)
  
  # TODO: there is no predict!
  #pred = predict(model, newdata = newX)
  if (family$family == "gaussian") {
    pred = model$yhat.test.mean
  } else {
    # No mean is provided for binary Y :/
    pred = colMeans(pnorm(model$yhat.test))
  }
  
  fit = list(object = model)
  class(fit) = c("SL.dbarts")
  out = list(pred = pred, fit = fit)
  return(out)
}



my.dbarts2 = function(Y, X, newX, family, obsWeights, id,
                      sigest = NA,
                      sigdf = 3,
                      sigquant = 0.90,
                      k = 2.0,  # xtan
                      power = 2.0,
                      base = 0.95,
                      binaryOffset = 0.0,
                      ntree = 200,
                      ndpost = 1000,
                      nskip = 100,
                      printevery = 100,
                      keepevery = 1,
                      keeptrainfits = T,
                      usequants = F,
                      numcut = 100,
                      printcutoffs = 0,
                      nthread = 1,
                      keepcall = T,
                      verbose = F,
                      ...) {
  
  # .SL.require("dbarts")
  
  model =
    dbarts::bart(x.train = X,
                 y.train = Y,
                 # We need to pass newX in directly due to lack of prediction.
                 x.test = newX,
                 sigest = sigest,
                 sigdf = sigdf,
                 sigquant = sigquant,
                 k = k,
                 power = power,
                 base = base,
                 binaryOffset = binaryOffset,
                 weights = obsWeights,
                 ntree = ntree,
                 ndpost = ndpost,
                 nskip = nskip,
                 printevery = printevery,
                 keepevery = keepevery,
                 keeptrainfits = keeptrainfits,
                 usequants = usequants,
                 numcut = numcut,
                 printcutoffs = printcutoffs,
                 nthread = nthread,
                 keepcall = keepcall,
                 verbose = verbose)
  
  # TODO: there is no predict!
  #pred = predict(model, newdata = newX)
  if (family$family == "gaussian") {
    pred = model$yhat.test.mean
  } else {
    # No mean is provided for binary Y :/
    pred = colMeans(pnorm(model$yhat.test))
  }
  
  fit = list(object = model)
  class(fit) = c("SL.dbarts")
  out = list(pred = pred, fit = fit)
  return(out)
}
