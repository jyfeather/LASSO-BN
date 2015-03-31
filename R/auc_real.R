rm(list = ls())

require("genlasso")   # Lasso solver
require("ROCR")       # ROC
require("Matrix")

set.seed(2015)

kIteration <- 200
node.num <- 22
sig.set <- c(0.1, 0.3, 0.5, 0.7, 1, 1.5) # Mean shift magnitude
var.df <- 1 # 1, 2, 3, 4, 5 guessed amount of mean shift vars
ns <- 1 # 1, 2, 5, 10

load("./dat/real/weighM")
load("./dat/real/shifts")

## LASSO-BN
auc <- matrix(data = NA, nrow = length(sig.set), ncol = length(shifts.pos))
for (pos in 1:length(shifts.pos)) {
  shift.real <- rep(0, node.num)
  shift.real[shifts.pos[pos]] = 1
  for (i in 1:length(sig.set)) {
    load(paste("./dat/real/dat",sig.set[i], shifts.pos[pos], sep = "_"))
    dat <- dat[1:(kIteration*ns),]
    size <- nrow(dat)
    tmp.coef <- matrix(data=0, nrow=node.num, ncol=size)
    tmp.x <- diag(node.num)
    tmp.least <- solve(t(tmp.x) %*% tmp.x) %*% tmp.x
    for (m in 1:size) {
      tmp.y <- solve(W) %*% dat[m,]
      tmp.least2 <- tmp.least %*% tmp.y        
      tmp.coef[sort(abs(tmp.least2), decreasing=T, index.return=T)$ix[1:var.df], m] = 1
    }    
    shift.l1 <- rep(0, node.num)
    for (k in 1:node.num) {
      shift.l1[k] <- nnzero(tmp.coef[k,]) / size
    }    
    roc.pred <- prediction(shift.l1, shift.real)
    roc.perf <- performance(roc.pred, "auc")
    auc[i, pos] = as.numeric(roc.perf@y.values)  
  }
}
print(auc)
write.csv(auc, file="./dat/auc.csv")

## VS-MSPC
auc <- matrix(data = NA, nrow = length(sig.set), ncol = length(shifts.pos))
for (pos in 1:length(shifts.pos)) {
  shift.real <- rep(0, node.num)
  shift.real[shifts.pos[pos]] = 1
  for (i in 1:length(sig.set)) {
    load(paste("./dat/real/dat",sig.set[i], shifts.pos[pos], sep = "_"))
    dat <- dat[1:(kIteration*ns),]
    size <- nrow(dat)
    W <- t(chol(cov(dat)))
    tmp.coef <- matrix(data=0, nrow=node.num, ncol=size)
    tmp.x <- solve(W)
    for (m in 1:size) {
      tmp.y <- tmp.x %*% dat[m,]
      fit <- genlasso(tmp.y, tmp.x, diag(node.num))
      tmp.coef[sort(fit$beta[,var.df+1], decreasing = T, index.return = T)$ix[1:var.df], m] <- 1
    }    
    shift.l1 <- rep(0, node.num)
    for (k in 1:node.num) {
      shift.l1[k] <- nnzero(tmp.coef[k,]) / size  
    }    
    roc.pred <- prediction(shift.l1, shift.real)
    roc.perf <- performance(roc.pred, "auc")
    auc[i, pos] = as.numeric(roc.perf@y.values)  
  }
}
print(auc)
write.csv(auc, file="./dat/auc.csv")
