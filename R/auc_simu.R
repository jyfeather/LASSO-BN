rm(list = ls())

setwd("C:\\Users\\JIN\\Dropbox\\Research\\Monitoring\\code")
#setwd("C:\\Users\\jyfea_000\\Dropbox\\Research\\Monitoring\\code")

## Import packages
require("pcalg")      # Causal network
require("genlasso")   # Lasso solver
require("ROCR")       # ROC
require("Matrix")

## Import functions
source("./R/chart.R")

set.seed(2014)

node.num.set <- c(30, 60, 100) # Amount of nodes in BN
kIteration <- 500 # if ARL = 200, type I error number is 10000/200 = 50
sig.set <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.5) # Mean shift magnitude
var.df <- c(2, 3, 4, 5, 6) # 2, 3, 4, 5, 6 guessed amount of mean shift vars
#type = "LASSO-BN" # VS-MSPC or LASSO-BN
type = "VS-MSPC" # VS-MSPC or LASSO-BN

auc = array(data=NA, dim=c(length(node.num.set), length(sig.set), length(var.df)))
for (i in 1:length(node.num.set)) {
  for (j in 1:length(sig.set)) {
    load(paste("./data/4ROC/dat/dat", node.num.set[i], sig.set[j], sep="_"))
    load(paste("./data/4ROC/dat/shift", node.num.set[i], sig.set[j], sep="_"))
    load(paste("./data/4ROC/dat/weight", node.num.set[i], sig.set[j], sep="_"))
    #cov.est <- cov(dat)
    #W <- chol(cov.est)
    for (k in 1:length(var.df)) {
      tmp.coef <- matrix(data=NA, nrow=node.num.set[i], ncol=kIteration)
      for (m in 1:kIteration) {
        tmp.fit <- solver(W, dat[m,], var.df[k], type)
        tmp.coef[,m] = tmp.fit$coef.set  
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
}
print(auc)