rm(list = ls())

setwd("C:\\Users\\Administrator\\Dropbox\\Research\\Monitoring\\code")
#setwd("C:\\Users\\JIN\\Dropbox\\Research\\Monitoring\\code")

## Import packages
require("pcalg")      # Causal network
require("genlasso")   # Lasso solver

## Import functions
source("./R/chart.R")
source("./R/ARL.R")

set.seed(2014)

node.num.set <- c(30, 60, 100) # Amount of nodes in BN
var.net <- 1 # 1, 2, 3
load(paste("./data/4ARL/bn/bn", node.num.set[var.net], sep="")) # load Bayesian Network
#plot(bn.dag)

kIteration <- 10000 # if ARL = 200, type I error number is 10000/200 = 50
sig.set <- c(0, 0.5, 1, 2, 4) # Mean shift magnitude
var.df <- 2 # 2, 3, 4, 5, 6 guessed amount of mean shift vars
 
#################### Hotelling T squared chart ###########################
load(paste("./data/4ARL/dat/dat", node.num.set[var.net], 0, sep="_"))
control.limit <- ControlLimitT2(dat)
rm(dat)

ARLs <- vector(mode="numeric", length=length(sig.set))
for (i in 1:length(sig.set)) {
  filename <- paste("./data/4ARL/dat/dat", node.num.set[var.net], sig.set[i], sep="_")
  load(filename)
  ARLs[i] = ARLT2(dat, control.limit, sig.set[i])
  rm(dat)
}

#################### Lasso Based (VS-MSPC) ###########################
load(paste("./data/4ARL/dat/dat", node.num.set[var.net], 0, sep="_"))
control.limit <- ControlLimitL1(dat, var.df)
rm(dat)

ARLs <- vector(mode="numeric", length=length(sig.set))
for (i in 1:length(sig.set)) {
  filename <- paste("./data/4ARL/dat/dat", node.num.set[var.net], sig.set[i], sep="_")
  load(filename)
  ARLs[i] = ARLL1(dat, control.limit, sig.set[i], var.df)
  print(i)
  rm(dat) 
}

#################### Bayesian Network Lasso Based (BNL1) ###########################
load(paste("./data/4ARL/dat/dat", node.num.set[var.net], 0, sep="_"))
load(paste("./data/4ARL/bn//bn", node.num.set[var.net], sep=""))
cov.true <- trueCov(bn.dag)
weight <- WeightConvert(node.num.set[var.net], cov.true, bn.dag)
control.limit <- ControlLimitBNL1(dat, var.df, weight, cov.true)
rm(dat)

ARLs <- vector(mode="numeric", length=length(sig.set))
for (i in 1:length(sig.set)) {
  filename <- paste("./data/4ARL/dat/dat", node.num.set[var.net], sig.set[i], sep="_")
  load(filename)
  ARLs[i] = ARLBNL1(dat, control.limit, sig.set[i], var.df, weight, cov.true)
  print(i)
  rm(dat) 
}
rm(bn.dag)