rm(list = ls())

setwd("C:/Users/JIN/Documents/GitHub/LASSO-BN")

## Import packages
require("pcalg")      # Causal network
require("genlasso")   # Lasso solver

source("./R/chart.R")
source("./R/ARL.R")

set.seed(2014)

node.num.set <- c(30, 60, 100) # Amount of nodes in BN
var.net <- 2 # 1, 2, 3
kIteration <- 10000 # if ARL = 200, type I error number is 10000/200 = 50
sig.set <- c(0.5, 1, 2, 4) # Mean shift magnitude
var.df <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, node.num.set[var.net]/2, node.num.set[var.net]-1) # guessed amount of mean shift vars

load(paste("./data/Simu/dat/dat", node.num.set[var.net], 0, sep="_"))
load(paste("./data/Simu/bn/bn", node.num.set[var.net], sep=""))
#control.limit <- ControlLimitBNL1(dat, var.df, bn.dag)
control.limit <- c(10.40094, 13.78475, 16.62168, 19.16862, 21.78490, 24.54558, 26.74896, 28.92770, 31.69569, 66.83572, 90.57255)
rm(dat)

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
