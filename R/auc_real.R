rm(list = ls())

setwd("C:\\Users\\JIN\\Dropbox\\Research\\Monitoring\\code")
#setwd("C:\\Users\\jyfea_000\\Dropbox\\Research\\Monitoring\\code")

## Import packages
#require("lars")       # Lasso solver
require("genlasso")   # Lasso solver
require("ROCR")       # ROC
require("Matrix")

## Import functions
source("./R/chart.R")

set.seed(2014)

kIteration <- 500
node.num <- 22
sig.set <- c(0.1, 0.3, 0.5, 0.7, 0.9, 1.5) # Mean shift magnitude
#sig.set <- c(0.1) # Mean shift magnitude
var.df <- c(1,2,3,4,5) # 1, 2, 3, 4, 5 guessed amount of mean shift vars
type = "LASSO-BN" # VS-MSPC or LASSO-BN
#type = "VS-MSPC" # VS-MSPC or LASSO-BN

## Solve model
auc <- matrix(data = NA, nrow = length(sig.set), ncol = length(var.df))
for (i in 1:length(sig.set)) {
  dat <- read.csv(paste("./data/RealExample/dat",sig.set[i], ".csv", sep = ""))
  dat <- dat[, -1] # 1st col is index
  load("./data/RealExample/weighM")
  #cov.est <- cov(dat)
  #W <- chol(cov.est)
  for (j in 1:length(var.df)) {
    tmp.coef <- matrix(data=NA, nrow=node.num, ncol=kIteration)
    for (k in 1:kIteration) {
      tmp.fit <- solver(W, t(dat[k,]), var.df[j], type)
      tmp.coef[,k] = tmp.fit$coef.set  
    }
    shift.real <- c(1, rep(0, 21))
    shift.l1 <- rep(0, node.num)
    for (k in 1:node.num) {
      shift.l1[k] <- nnzero(tmp.coef[k,]) / kIteration  
    }    
    print(shift.l1)
    roc.pred <- prediction(shift.l1, shift.real)
    roc.perf <- performance(roc.pred, "auc")
    auc[i, j] = as.numeric(roc.perf@y.values)  
  }  
}
print(auc)