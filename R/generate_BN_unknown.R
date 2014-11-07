rm(list = ls())

setwd("C:/Users/JIN/Documents/GitHub/LASSO-BN")

## Import packages
require("pcalg")      # Causal network
require("MASS")

set.seed(2014)

kIteration <- 10000 
#sig.set <- c(0, 0.5, 1, 2, 4) # Mean shift magnitude
sig.set <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.5, 2, 4) # Mean shift magnitude
kNumShifts <- 3 # Amount of mean shift varibles

## Generate and save a node number specified random Bayesian Network
node.num.set <- c(30, 60, 100) # Amount of nodes in BN
edge.num.set <- c(0.15, 0.07, 0.04) # Related to number of arcs For ARL
#edge.num.set <- c(0.3, 0.15, 0.05) # Related to number of arcs For AUC 
for (i in 1:length(node.num.set)) {
  bn.dag <- randomDAG(n = node.num.set[i], prob = edge.num.set[i])  
  save(bn.dag, file=paste("./data/Simu/bn/bn", node.num.set[i], sep=""))
  
  W <- t(chol(trueCov(bn.dag))) # find causal effect matrix
  save(W, file=paste("./data/Simu/dat/weighM", node.num.set[i], sep="_"))

  shifts <- rep(0, node.num.set[i])
  shift.pos <- sample(1:node.num.set[i], kNumShifts) # for ARL
  save(shift.pos, file=paste("./data/Simu/dat/shift", node.num.set[i], sep="_"))
    
  ## Generate and save the dataset based on BN
  for (j in 1:length(sig.set)) {
    shifts[shift.pos] <- sig.set[j]
    e.mat <- mvrnorm(kIteration, mu=shifts, Sigma=diag(node.num.set[i])) # simulate true meam shift
    dat <- W %*% t(e.mat)
    dat <- t(dat)
    save(dat, file=paste("./data/Simu/dat/dat", node.num.set[i], sig.set[j], sep="_"))
  }
}
