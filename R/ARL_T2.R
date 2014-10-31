rm(list = ls())

#setwd("C:\\Users\\Administrator\\Dropbox\\Research\\Monitoring\\code")
setwd("C:\\Users\\JIN\\Dropbox\\Research\\Monitoring\\code")

## Import functions
source("./R/ARL.R")

set.seed(2014)

node.num.set <- c(30, 60, 100) # Amount of nodes in BN
var.net <- 2 # 1, 2, 3
kIteration <- 10000 # if ARL = 200, type I error number is 10000/200 = 50
sig.set <- c(0, 0.5, 1, 2, 4) # Mean shift magnitude
 
## Hotelling T squared chart
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