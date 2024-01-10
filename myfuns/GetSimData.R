GetSimData <- function(n, m, tau, set, is_homo=TRUE) {
  X = matrix(runif(n * m), n, m)
  X = as.data.frame(X)
  colnames(X) = paste0("X", seq(m))
  
  W1 <- exp(X$X1 / 2)
  W2 <- exp(X$X2 / 3)
  W3 <- X$X3^2
  W4 <- X$X4^2
  W5 <- X$X5
  W6 <- X$X6
  W7 <- X$X7 + X$X8
  W8 <- X$X7^2 + X$X8^2
  W9 <- X$X9^3

  
  W <- cbind(W1,W2,W3,W4,W5,W6,W7,W8,W9)
  W <- scale(W)  # standardize to have mean 0 and sd 1
  W <- as.data.frame(W)
  rm(W1,W2,W3,W4,W5,W6,W7,W8,W9)
  # colMeans(W)
  
  if (set == 1) {  # normal ps
    ps_logit <- (-3 + -1*W$W1 + 2*W$W2 - 3*W$W3 + 3*W$W4 + 2*W$W5 + W$W6) / 15
  } else if (set == 2) {  # moderate ps
    ps_logit <- (-3 + -1*W$W1 + 2*W$W2 - 3*W$W3 + 3*W$W4 + 2*W$W5 + W$W6) / 10
  } else if (set == 3) {  # extreme ps
    ps_logit <- (-3 + -1*W$W1 + 2*W$W2 - 3*W$W3 + 3*W$W4 + 2*W$W5 + W$W6) / 5
  } else if (set == 4) {  # more extreme
    ps_logit <- (-8*W$W1 + 1.5*W$W2 + 0.5*W$W3 - 0.5*W$W4 + 2.5*W$W5 - 0.5*W$W6) / 5
  }
  ps <- 1 / (1 + exp(-ps_logit))
  A <- rbinom(n, 1, ps)
  Y0 <- -2 + 1.5*W$W1 - 2*W$W2 + 1.5*W$W3 + 2.5*W$W7 - W$W8 + W$W9 + rnorm(n)
  
  if (is_homo) {
    Y1 <- Y0 + tau
  } else {  # hetero
    Y1 <- Y0 + tau + 5*W$W1 + 3*W$W3 + 2*W$W1*W$W3#2
  }
  
  Y <- Y1 * A + Y0 * (1 - A)
  
  # hist(ps)
  # hist(Y0)
  
  mydat <- cbind(Y, A, Y0, Y1, ps, ps_logit, X, W)
  mydat <- as.data.frame(mydat)
  
  #cut in 5 folds for CV in SuperLearner (add a column `fold_idx` in mydat)
  num_fold <- 5
  contlID <- which(mydat[, "A"]==0)
  treatID <- which(mydat[, "A"]==1)
  dat0 <- mydat[contlID, ]
  dat1 <- mydat[treatID, ]
  
  fold0_idx <- sample(1:num_fold,nrow(dat0),replace=TRUE)
  fold1_idx <- sample(1:num_fold,nrow(dat1),replace=TRUE)
  
  dat0$fold_idx <- fold0_idx
  dat1$fold_idx <- fold1_idx
  
  mydat <- as.data.frame(rbind(dat0,dat1))
  mydat <- mydat[ order(as.numeric(rownames(mydat))) , ]
  
  print(table(dat0$fold_idx))
  print(table(dat1$fold_idx))
  
  return(mydat)
}
