rm(list = ls())

library(pcalg)
set.seed(2015)

node.num.set <- c(30, 50, 100) # Amount of nodes in BN
var <- 1
kIteration <- 10000 # if ARL = 200, type I error number is 10000/200 = 50
sig.set <- c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5) # Mean shift magnitude
var.df <- c(2, 3, 4, 5) # guessed amount of mean shift vars

load(paste("./dat/simu/bn", node.num.set[var], sep=""))
load(paste("./dat/simu/weighM", node.num.set[var], sep="_"))

df.len <- length(var.df)

## control limit, ARL0
load(paste("./dat/simu/dat", node.num.set[var], 0, sep="_"))
inv <- solve(cov(dat))
size <- nrow(dat)  
res.set <- vector(mode = "numeric", length = size)
for (i in 1:size) {
  T2 <- dat[i,] %*% inv %*% dat[i,]
  res.set[i] <- T2
}
res.set <- sort(res.set)
kARL0 = 200
flag <- size - size / kARL0 
cl <- res.set[flag] 
rm(dat, res.set, T2, kARL0, flag, inv)

## out of control to ARL1
ARLs <- vector(mode = "numeric", length = length(sig.set))
for (i in 2:length(sig.set)) {
  load(paste("./dat/simu/dat", node.num.set[var], sig.set[i], sep="_"))
  inv <- solve(cov(dat))
  size <- nrow(dat)
  err.num <- 0 
  for (j in 1:size) {
    T2 <- dat[j,] %*% inv %*% dat[j,]
    if (T2 <= cl) err.num <- err.num + 1
  }
  ARLs[i] = size/(size-err.num)
}
ARLs[1] <- 200
write.csv(ARLs, file = "./dat/ARL.csv")
