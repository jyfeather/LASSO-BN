rm(list = ls())

####################### create 5 node BN ###########################
library(bnlearn)
library(igraph)
library(MASS)
library(graph)

# build BN structure
bn.struct = model2network("[Z1][Z4][Z2|Z1][Z3|Z1:Z4][Z5|Z2:Z3]")
g <- igraph.from.graphNEL(as.graphNEL(bn.struct))
plot(g)

# known BN path coefficient matrix, omega
numNode <- 5
omega <- matrix(0, nrow = numNode, ncol = numNode)
omega[2,1] <- 0.868
omega[3,1] <- 0.493
omega[3,4] <- 0.325
omega[5,2] <- 0.335
omega[5,3] <- 0.574
W <- solve(diag(numNode) - omega) # causal effect matrix

rm(bn.struct, g)

####################### create 32 shifting scenarios ###########################
numScene <- 32
scene.true <- t(combn(5,0,tabulate,nbins=5))
scene.true <- rbind(scene.true, t(combn(5,1,tabulate,nbins=5)))
scene.true <- rbind(scene.true, t(combn(5,2,tabulate,nbins=5)))
scene.true <- rbind(scene.true, t(combn(5,3,tabulate,nbins=5)))
scene.true <- rbind(scene.true, t(combn(5,4,tabulate,nbins=5)))
scene.true <- rbind(scene.true, t(combn(5,5,tabulate,nbins=5)))

sig <- 3 # mean shift scale
M <- 5000 # number of sample

###################### out of control ###########################
library(Matrix)
library(genlasso)

#type <- "VS-MSPC"
type <- "LASSO-BN"
scene.actual <- matrix(0, nrow = numScene, ncol = numNode)
i <- 2
while (i <= numScene) {
  print(i)
  emat <- mvrnorm(n = M, mu = scene.true[i,]*sig, Sigma = diag(numNode))
  dat <- t(W %*% t(emat))

  #guessNum <- nnzero(scene.true[i,]) # guess number of mean shift valuables
  #if (guessNum == 5) guessNum <- 4
  guessNum <- 2
  
  if (type == "LASSO-BN") {
    tmp.coef <- matrix(0, nrow = numNode, ncol = M)
    tmp.x <- diag(numNode)
    tmp.least <- solve(t(tmp.x) %*% tmp.x) %*% tmp.x
    for (j in 1:M) {
      tmp.y <- solve(W) %*% dat[j,]
      tmp.least2 <- tmp.least %*% tmp.y
      tmp.coef[sort(abs(tmp.least2), decreasing = T, index.return = T)$ix[1:guessNum], j] <- 1
    }
    for (j in 1:numNode) {
      scene.actual[i,j] <- nnzero(tmp.coef[j,])/M
    }  
  } else {
    tmp.coef <- matrix(0, nrow = numNode, ncol = M)
    tmp.x <- solve(W)
    for (j in 1:M) {
      tmp.y <- tmp.x %*% dat[j,]
      fit <- genlasso(tmp.y, tmp.x, diag(numNode))
      tmp.coef[sort(fit$beta[,guessNum+1], decreasing = T, index.return = T)$ix[1:guessNum], j] <- 1
    }
    for (j in 1:numNode) {
      scene.actual[i,j] <- nnzero(tmp.coef[j,])/M
    }
  }
  
  i <- i+1
}
scene.e <- abs(scene.true - scene.actual)
scene.actual2 <- round(scene.actual)

###################### trajectory of error rate ###########################
err.causationT2 <- read.csv(file = "./dat/Dr_Li_exmaple/CausationT2.csv", header = F)
err.MTY <- read.csv(file = "./dat/Dr_Li_exmaple/MTY.csv", header = F)
err.LassoBN <- read.csv(file = "./dat/Dr_Li_exmaple/LassoBN.csv", header = F)
err.VSMSPC <- read.csv(file = "./dat/Dr_Li_exmaple/VSMSPC.csv", header = F)

scene.err1 <- 1 - scene.true[-1,]
scene.err2 <- scene.true[-1,]

err1.causationT2 <- rowSums(scene.err1 * err.causationT2) / rowSums(scene.err1)
err1.MTY <- rowSums(scene.err1 * err.MTY) / rowSums(scene.err1)
err1.LassoBN <- rowSums(scene.err1 * err.LassoBN) / rowSums(scene.err1)
err1.VSMSPC <- rowSums(scene.err1 * err.VSMSPC) / rowSums(scene.err1)
err1 <- cbind(err1.LassoBN, err1.causationT2)
#err1[31,] <- 0

err2.causationT2 <- rowSums(scene.err2 * err.causationT2) / rowSums(scene.err2)
err2.MTY <- rowSums(scene.err2 * err.MTY) / rowSums(scene.err2)
err2.LassoBN <- rowSums(scene.err2 * err.LassoBN) / rowSums(scene.err2)
err2.VSMSPC <- rowSums(scene.err2 * err.VSMSPC) / rowSums(scene.err2)
err2 <- cbind(err2.LassoBN, err2.causationT2)

par(mfrow=c(1,2))
matplot(err1, type = c("l"), pch = 1, col = 1, ylab = "false positive rate")
legend("topright", c("LASSO-BN", "Causation-based T2"), lty = c(1,2,3))
matplot(err2, type = c("l"), pch = 1, col = 1, ylab = "false negative rate")
legend("topright", c("LASSO-BN", "Causation-based T2"), lty = c(1,2,3))
