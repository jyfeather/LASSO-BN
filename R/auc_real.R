rm(list = ls())

require("genlasso")   # Lasso solver
require("ROCR")       # ROC
require("Matrix")

set.seed(2015)

kIteration <- 500
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
    tmp.coef <- matrix(data=0, nrow=node.num, ncol=kIteration)
    tmp.x <- diag(node.num)
    tmp.least <- solve(t(tmp.x) %*% tmp.x) %*% tmp.x
    for (m in 1:kIteration) {
      tmp.coef2 <- matrix(0, nrow = node.num, ncol = ns)
      for (n in 1:ns) {
        l <- sample(nrow(dat), 1)
        tmp.y <- solve(W) %*% dat[l,]
        tmp.least2 <- tmp.least %*% tmp.y        
        tmp.coef2[, n] = abs(tmp.least2)
      }
      tmp.coef[sort(rowSums(tmp.coef2), decreasing = T, index.return = T)$ix[1:var.df], m] <- 1
    }
    
    shift.l1 <- rep(0, node.num)
    for (k in 1:node.num) {
      shift.l1[k] <- nnzero(tmp.coef[k,]) / kIteration  
    }    
    roc.pred <- prediction(shift.l1, shift.real)
    roc.perf <- performance(roc.pred, "auc")
    auc[i, pos] = as.numeric(roc.perf@y.values)  
  }
  print(auc)
}
write.csv(auc, file="./dat/auc.csv")


## VS-MSPC
auc <- matrix(data = NA, nrow = length(sig.set), ncol = length(shifts.pos))
for (pos in 1:length(shifts.pos)) {
  shift.real <- rep(0, node.num)
  shift.real[shifts.pos[pos]] = 1
  for (i in 1:length(sig.set)) {
    load(paste("./dat/real/dat",sig.set[i], shifts.pos[pos], sep = "_"))
    tmp.coef <- matrix(data=0, nrow=node.num, ncol=kIteration)
    tmp.x <- solve(W)
    for (m in 1:kIteration) {
      tmp.coef2 <- matrix(0, nrow = node.num, ncol = ns)
      for (n in 1:ns) {
        l <- sample(nrow(dat), 1)
        tmp.y <- solve(W) %*% dat[l,]
        fit <- genlasso(tmp.y, tmp.x, diag(node.num))
        tmp.coef2[,n] <- fit$beta[,var.df+1]  
      }
      tmp.coef[sort(rowSums(tmp.coef2), decreasing = T, index.return = T)$ix[1:var.df], m] <- 1
    }
    
    shift.l1 <- rep(0, node.num)
    for (k in 1:node.num) {
      shift.l1[k] <- nnzero(tmp.coef[k,]) / kIteration  
    }    
    roc.pred <- prediction(shift.l1, shift.real)
    roc.perf <- performance(roc.pred, "auc")
    auc[i, pos] = as.numeric(roc.perf@y.values)  
  }
  print(auc)
}
write.csv(auc, file="./dat/auc.csv")
