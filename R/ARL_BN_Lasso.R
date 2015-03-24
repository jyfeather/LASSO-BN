rm(list = ls())

library(pcalg)
set.seed(2015)

node.num.set <- c(30, 50, 100) # Amount of nodes in BN
kIteration <- 10000 # if ARL = 200, type I error number is 10000/200 = 50
sig.set <- c(0.5, 1, 2, 4) # Mean shift magnitude
var.df <- c(2, 3, 4, 5) # guessed amount of mean shift vars

## control limit, ARL0
load(paste("./dat/simu/dat", node.num.set[1], 0, sep="_"))
load(paste("./dat/simu/bn", node.num.set[1], sep=""))
load(paste("./dat/simu/weighM", node.num.set[1], sep="_"))

size <- nrow(dat)  
df.len <- length(var.df)
res.set <- matrix(data=NA, nrow=size, ncol=df.len)
inv <- solve(trueCov(bn.dag))
for (i in 1:size) {
  T2 <- dat[i,] %*% inv %*% dat[i,]
  # solve S2
  common <- diag(ncol(dat)) # (x^Tx)^-1x^T
  least <- common %*% dat[i,]
  lambdas <- sort(abs(least), decreasing = T)[var.df+1]
  for (j in 1:length(var.df)) {
    beta <- sign(least)*((abs(least)-lambdas[j]) * ((abs(least)-lambdas[1])>0))
    S2 <- sum((solve(W) %*% dat[i,] - beta)^2) + lambdas[j]*sum(abs(beta))
    res.set[i,j] <- T2 -S2
  }
}
for (j in 1:df.len) {
  res.set[,j] <- sort(res.set[,j])
}
kARL0 = 200
flag <- size - size / kARL0 
res.set[flag,] 
rm(dat)

## out of control to ARL1
ARLs <- matrix(data=NA, nrow=length(sig.set), ncol=length(var.df))
for (i in 1:length(sig.set)) {
  filename <- paste("./data/Simu/dat/dat", node.num.set[var.net], sig.set[i], sep="_")
  load(filename)
  ARLs[i,] = ARLBNL1(dat, control.limit, sig.set[i], var.df, bn.dag)
  print(i)
  rm(dat) 
}
rm(bn.dag)
write.csv(t(ARLs), file = "./data/ARL.csv")
