rm(list = ls())

library(pcalg)
library(genlasso)
set.seed(2015)

node.num.set <- c(30, 50, 100) # Amount of nodes in BN
var <- 1 
kIteration <- 10000 # if ARL = 200, type I error number is 10000/200 = 50
sig.set <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5) # Mean shift magnitude
var.df <- c(2, 3, 4, 5) # guessed amount of mean shift vars

load(paste("./dat/simu/weighM", node.num.set[var], sep="_"))

df.len <- length(var.df)

## control limit, ARL0
load(paste("./dat/simu/dat", node.num.set[var], 0, sep="_"))
inv <- solve(cov(dat))
size <- nrow(dat)  
res.set <- matrix(data=NA, nrow=size, ncol=df.len)
X <- solve(W)
for (i in 1:size) {
  T2 <- dat[i,] %*% inv %*% dat[i,]  
  y <- X %*% dat[i,]
  fit <- genlasso(y, X, diag(node.num.set[var])/2)
  for (j in 1:length(var.df)) {
    beta <- round(fit$beta[,var.df[j]+1], digits = 3)
    S2 <- sum((y - X %*% beta)^2) #+ fit$lambda[var.df[j]+1]*sum(abs(beta))
    res.set[i,j] <- T2 -S2
  }
}
for (j in 1:df.len) {
  res.set[,j] <- sort(res.set[,j])
}
kARL0 = 200
flag <- size - size / kARL0 
cl <- res.set[flag,] 
rm(dat, res.set, T2, y, fit, S2, kARL0, flag)

## out of control to ARL1
ARLs <- matrix(data=NA, nrow=length(sig.set), ncol=length(var.df))
for (i in 2:length(sig.set)) {
  load(paste("./dat/simu/dat", node.num.set[var], sig.set[i], sep="_"))
  inv <- solve(cov(dat))
  size <- nrow(dat)
  err.num <- rep(0, df.len)
  X <- solve(W)
  for (j in 1:size) {
    T2 <- dat[j,] %*% inv %*% dat[j,]
    # solve S2
    y <- X %*% dat[j,]
    fit <- genlasso(y, X, diag(node.num.set[var])/2)
    for (k in 1:length(var.df)) {
      beta <- round(fit$beta[,var.df[k]+1], digits = 3)
      S2 <- sum((y - X %*% beta)^2) #+ fit$lambda[var.df[k]+1]*sum(abs(beta))
      if ((T2-S2)<cl[k]) err.num[k] <- err.num[k] + 1
    }
  }
  ARLs[i,] = size/(size-err.num)
  print(i)
}
ARLs[1,] <- rep(200, df.len)
write.csv(ARLs, file = "./dat/ARL.csv")
