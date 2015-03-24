# ARL0 = 1/alpha
# ARL1 = 1/(1-beta)
ControlLimitT2 <- function(dataset) {
  kIteration <- 10000
  cov <- cov(dataset)  
  inv <- solve(cov)
  T2s <- vector(mode="numeric", length=kIteration)
  for (i in 1:kIteration) {
    T2s[i] <- dataset[i,] %*% inv %*% dataset[i,] # Compute each chart stastics
  }
  T2s <- sort(T2s)
  kARL0 = 200
  flag <- kIteration - kIteration / kARL0 
  return(T2s[flag])  
}

ARLT2 <- function(dataset, control, sig) {
  cov <- cov(dataset)  
  inv <- solve(cov)
  size <- nrow(dataset)
  err.num <- 0
  ARL <- 0
  if (sig == 0) {
    for (i in 1:size) {
      T2 <- dat[i,] %*% inv %*% dat[i,]
      if (T2 > control) err.num <- err.num + 1
    }
    ARL <- size / err.num   
  } else if (sig > 0) {
    for (i in 1:size) {
      T2 <- dat[i,] %*% inv %*% dat[i,]
      if (T2 <= control) err.num <- err.num + 1
    }
    ARL <- size / (size - err.num)
  } else {
    stop("Error Sig!")
  }
  return(ARL)
}

ControlLimitBNL1 <- function(dataset, df, wei) {
  size <- nrow(dataset)  
  df.len <- length(df)
  res.set <- matrix(data=NA, nrow=size, ncol=df.len)
  for (i in 1:size) {
    T2 <- dataset[i,] %*% inv %*% dataset[i,]
    # solve S2
    common <- diag(ncol(dataset)) # (x^Tx)^-1x^T
    least <- common %*% dataset[i,]
    lambdas <- sort(abs(least), decreasing = T)$x[df+1]



  }
  for (j in 1:df.len) {
    res.set[,j] <- sort(res.set[,j])
  }
  kARL0 = 200
  flag <- size - size / kARL0 
  return(res.set[flag,]) 
}

ARLBNL1 <- function(dataset, control, sig, df, dag) {
  size <- nrow(dataset)  
  cov <- trueCov(dag) # true covariance matrix 
  inv <- solve(cov)
  wei.mat <- t(chol(cov)) # true weight matrix
  df.len <- length(df)
  err.num <- rep(0, df.len)
  ARL <- rep(0, df.len) 
  if (sig == 0) {
    for (i in 1:size) {
      T2 <- dat[i,] %*% inv %*% dat[i,]
      S2 <- solver(wei.mat, dataset[i,], type = "LASSO-BN")
      for (j in 1:df.len) {
        res <- T2 - S2$obj.val[df[j] + 1]  
        if (res > control[j]) err.num[j] <- err.num[j] + 1
      }
    }
    ARL <- size / err.num   
  } else if (sig > 0) {
    for (i in 1:size) {
      T2 <- dat[i,] %*% inv %*% dat[i,]
      S2 <- solver(wei.mat, dataset[i,], type = "LASSO-BN")
      for (j in 1:df.len) {
        res <- T2 - S2$obj.val[df[j] + 1] # because first row is df 0  
        if (res <= control[j]) err.num[j] <- err.num[j] + 1
      }
    }
    ARL <- size / (size - err.num)
  } else {
    stop("Error Sig!")
  }
  return(ARL)
}