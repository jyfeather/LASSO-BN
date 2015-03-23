rm(list = ls())

require("ROCR")       # ROC
require("Matrix")

set.seed(2015)

node.num.set <- c(30, 50, 100) # Amount of nodes in BN
kIteration <- 500 # if ARL = 200, type I error number is 10000/200 = 50
sig.set <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.5) # Mean shift magnitude
var.df <- c(1, 2, 5, 10) # guessed amount of mean shift vars
#type = "LASSO-BN" # VS-MSPC or LASSO-BN
type = "VS-MSPC" # VS-MSPC or LASSO-BN

auc = array(data=NA, dim=c(length(node.num.set), length(sig.set), length(var.df)))
for (i in 1:length(node.num.set)) {
  load(paste("./dat/simu/shift", node.num.set[i], sep="_"))
  load(paste("./dat/simu/weighM", node.num.set[i], sep="_"))
  
  for (j in 1:length(sig.set)) {
    load(paste("./dat/simu/dat", node.num.set[i], sig.set[j], sep="_"))
    for (k in 1:length(var.df)) {
      tmp.coef <- matrix(data=0, nrow=node.num.set[i], ncol=kIteration)
      if (type == "LASSO-BN") {
        tmp.x <- diag(node.num.set[i])
      } else if (type == "VS-MSPC") {
        tmp.x <- solve(W)
      } else {
        stop("Wrong Type!")
      }
      tmp.least <- solve(t(tmp.x) %*% tmp.x) %*% tmp.x
      for (m in 1:kIteration) {
        tmp.y <- solve(W) %*% dat[m,]
        tmp.least2 <- tmp.least %*% tmp.y        
        tmp.coef[sort(abs(tmp.least2), decreasing = T, index.return=T)$ix[1:k], m] = 1
      }
      shift.real <- rep(0, node.num.set[i])
      shift.real[shift.pos] = 1 # 1 is true, 0 is false 
      shift.l1 <- rep(0, node.num.set[i])
      for (m in 1:node.num.set[i]) {
        shift.l1[m] <- nnzero(tmp.coef[m,]) / kIteration  
      }    
#      print(shift.l1)
      roc.pred <- prediction(shift.l1, shift.real)
      roc.perf <- performance(roc.pred, "auc")
      auc[i,j,k] = as.numeric(roc.perf@y.values)
    }
  }
  print(i)
}
print(auc)
write.csv(auc[1,,], file = "./dat/simu/auc.csv")

## test
#library(genlasso)
#tmp <- genlasso(y = solve(W)%*%dat[m,], X = diag(100), D = diag(100))
#sort(abs(coef.genlasso(tmp, df = 5)$beta), decreasing = T, index.return=T)$ix[1:5]
#tmp <- solve(t(tmp.x) %*% tmp.x) %*% tmp.x %*% (solve(W)%*%dat[m,])
#sort(abs(tmp), decreasing = T, index.return=T)$ix[1:5]
