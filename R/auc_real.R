rm(list = ls())

setwd("C:/Users/JIN/Documents/GitHub/LASSO-BN")

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
var.df <- c(1,2,3,4,5) # 1, 2, 3, 4, 5 guessed amount of mean shift vars
type = "LASSO-BN" # VS-MSPC or LASSO-BN
#type = "VS-MSPC" # VS-MSPC or LASSO-BN

## Solve model
auc <- matrix(data = NA, nrow = length(sig.set), ncol = length(var.df))

load("./data/TEP/weighM")
load("./data/TEP/dat/shifts")
shift.real <- rep(0, node.num)
shift.real[shifts.pos] = 1

for (i in 1:length(sig.set)) {
  load(paste("./data/TEP/dat/dat",sig.set[i], sep = ""))
  for (j in 1:length(var.df)) {
    tmp.coef <- matrix(data=NA, nrow=node.num, ncol=kIteration)
    for (k in 1:kIteration) {
      tmp.fit <- solver(W, dat[k,], var.df[j], type)
      tmp.coef[,k] = tmp.fit$coef.set  
    }
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
write.csv(auc, file="./data/auc.csv")
