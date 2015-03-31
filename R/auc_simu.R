rm(list = ls())

require("ROCR")       # ROC
require("Matrix")
require("genlasso")

set.seed(2015)

node.num.set <- c(30, 50, 100) # Amount of nodes in BN
kIteration <- 200 
sig.set <- c(0.1, 0.3, 0.5, 0.7, 1, 1.5) # Mean shift magnitude
var.df <- 3 # guessed amount of mean shift vars
ns <- 2 # 1, 2, 5, 10

## LASSO BN
auc = array(data=NA, dim=c(length(node.num.set), length(sig.set)))
for (i in 1:length(node.num.set)) {
  load(paste("./dat/simu/shift", node.num.set[i], sep="_"))
  load(paste("./dat/simu/weighM", node.num.set[i], sep="_"))
  for (j in 1:length(sig.set)) {
    load(paste("./dat/simu/dat", node.num.set[i], sig.set[j], sep="_"))
    dat <- dat[1:(kIteration*ns),]
    size <- nrow(dat)
    tmp.coef <- matrix(data=0, nrow=node.num.set[i], ncol=size)
    tmp.x <- diag(node.num.set[i])
    tmp.least <- solve(t(tmp.x) %*% tmp.x) %*% tmp.x
    for (m in 1:size) {
      tmp.y <- solve(W) %*% dat[m,]
      tmp.least2 <- tmp.least %*% tmp.y        
      tmp.coef[sort(abs(tmp.least2), decreasing=T, index.return=T)$ix[1:var.df], m] = 1
    }
    shift.real <- rep(0, node.num.set[i])
    shift.real[shift.pos] = 1 # 1 is true, 0 is false 
    shift.l1 <- rep(0, node.num.set[i])
    for (m in 1:node.num.set[i]) {
      shift.l1[m] <- nnzero(tmp.coef[m,]) / size
    }    
    roc.pred <- prediction(shift.l1, shift.real)
    roc.perf <- performance(roc.pred, "auc")
    auc[i,j] = as.numeric(roc.perf@y.values)
  }
  print(i)
}
print(t(auc))
write.csv(t(auc), file = "./dat/auc.csv")

## VS-MSPC
auc = array(data=NA, dim=c(length(node.num.set), length(sig.set)))
for (i in 1:length(node.num.set)) {
  load(paste("./dat/simu/shift", node.num.set[i], sep="_"))
  for (j in 1:length(sig.set)) {
    load(paste("./dat/simu/dat", node.num.set[i], sig.set[j], sep="_"))
    dat <- dat[1:(kIteration*ns),]
    W <- t(chol(cov(dat)))
    size <- nrow(dat)
    tmp.coef <- matrix(data=0, nrow=node.num.set[i], ncol=size)
    tmp.x <- solve(W)
    for (m in 1:size) {
      tmp.y <- solve(W) %*% dat[m,]
      fit <- genlasso(tmp.y, tmp.x, diag(node.num.set[i])/2)
      tmp.coef[sort(fit$beta[,var.df+1], decreasing = T, index.return = T)$ix[1:var.df], m] <- 1
    }
    shift.real <- rep(0, node.num.set[i])
    shift.real[shift.pos] = 1 # 1 is true, 0 is false 
    shift.l1 <- rep(0, node.num.set[i])
    for (m in 1:node.num.set[i]) {
      shift.l1[m] <- nnzero(tmp.coef[m,]) / size  
    }    
    roc.pred <- prediction(shift.l1, shift.real)
    roc.perf <- performance(roc.pred, "auc")
    auc[i,j] = as.numeric(roc.perf@y.values)
  }
  print(i)
}
print(auc)
write.csv(t(auc), file = "./dat/auc.csv")


## test
#library(genlasso)
#tmp <- genlasso(y = solve(W)%*%dat[m,], X = diag(100), D = diag(100))
#sort(abs(coef.genlasso(tmp, df = 5)$beta), decreasing = T, index.return=T)$ix[1:5]
#tmp <- solve(t(tmp.x) %*% tmp.x) %*% tmp.x %*% (solve(W)%*%dat[m,])
#sort(abs(tmp), decreasing = T, index.return=T)$ix[1:5]
