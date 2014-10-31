rm(list = ls())

#setwd("C:\\Users\\Administrator\\Dropbox\\Research\\Monitoring\\code")
setwd("C:\\Users\\JIN\\Dropbox\\Research\\Monitoring\\code")

## Import packages
require("pcalg")      # Causal network

set.seed(2014)

bn.dag <- randomDAG(n = 8, prob = 0.25)  
plot(bn.dag)

cov.true <- trueCov(bn.dag)
cor.true <- cov2cor(cov.true)

ida(1, 4, cov.true, bn.dag) # estimated total causal effects
causalEffect(bn.dag, 4, 1) # true total causal effects
par.total <- matrix(data=NA, nrow=8, ncol=8)
for (i in 1:8) {
  par.total[,i] = round(idaFast(i, c(1:8), cov.true, bn.dag), digits = 3)
}

par.mat <- wgtMatrix(bn.dag) # path coefficient of BN
par.mat[3,1] * par.mat[4,3] # should be equal to true total causal effects

solve(diag(8) - par.mat) == par.total # two approaches are equal

%solve(cor.true[c(1,5,6), c(1,5,6)]) %*% cor.true[c(1,5,6),7]

# covariance of data should be equal W * W^T 
round(par.total %*% t(par.total), digits = 2) == round(cov.true, digits = 2)
