rm(list = ls())


library(bnlearn)
library(igraph)
library(MASS)
library(graph)

# Create the Bayesian network structure
# BN structure based on Mississippi Eastman Process
bn.struct = model2network("[F1][F2][F3][F4][F19][F10|F1:F2:F3:F4][J20|F10][F5|J20][F6|F1:F2:F3:F4:F5]
                   [P7|F6][L8|F6:P7][T9|P7][L12|L8][T21|T9][T11|T9][T22|T11][P13|T11:L12]
                   [F14|P13][F17|F14][L15|F17][P16|L15:F19][T18|P16:F19]")

#plot(bn, highlight = c("F5", mb(bn, "F5")))
g <- igraph.from.graphNEL(as.graphNEL(bn.struct))
plot(g)

# Estimate the Bayesian network parameters
inputMat <- read.table("./dat/real/incontrol.csv", header = T)
inputMat <- as.data.frame(scale(inputMat))
nodeName <- colnames(inputMat)
nodeSymb <- c("F1", "F2", "F3", "F4", "F5", "F6",
	"P7", "L8", "T9", "F10", "T11", "L12", "P13",
	"F14", "L15", "P16", "F17", "T18", "F19", "J20",
	"T21", "T22")
colnames(inputMat) <- nodeSymb
fit = bn.fit(bn.struct, inputMat)

#rbn(fit, n = 10, inputMat)

########################################################
# The idea is like mu = W %*% b, mu is the exhibited mean shift, W is the causal effect matrix, b is the true mean shift
# To generate the data, as following steps
#    - step 1: compute the causal effect matrix from the learned bayesian network, W = (I - omega)^-1
#    - step 2: generate l according to multivariate normal distribution, with some of them has mean shift
#    - step 3: get the e according to the equation e = W %*% b

# Step 1
numNode <- length(nodeSymb)
omega <- matrix(data = 0, nrow = numNode, ncol = numNode)
for (i in 1:numNode) {
  tmp.coef <- coef(fit)[[nodeSymb[i]]]
  tmp.coef2 <- tmp.coef[-1]
  tmp.pos <- which(nodeSymb %in% names(tmp.coef2))
  omega[i, tmp.pos] = tmp.coef2
}
omega[is.na(omega)] <- 0
W <- solve(diag(numNode) - omega) # causal effect matrix
save(W, file = "./dat/real/weighM")

# shift positions
shifts.pos <- c(1,3,5)
save(shifts.pos, file = "./dat/real/shifts") 

# Step 2
sig.set <- c(0.1, 0.3, 0.5, 0.7, 1, 1.5)
for (i in 1:length(sig.set)) {
  for (j in 1:length(shifts.pos)) {
    shifts <- rep(0, numNode)
    shifts[shifts.pos[j]] <- sig.set[i]
    shifts.sigma <- diag(numNode) # identitiy covariance matrix
    numSample <- 10000 # sample size
    shifts.sample <- mvrnorm(n = numSample, mu = shifts, Sigma = shifts.sigma)
    
    # Step 3
    dat <- W %*% t(shifts.sample)
    dat <- t(dat)
    save(dat, file = paste("./dat/real/dat",sig.set[i], shifts.pos[j], sep = "_"))   
  }
}