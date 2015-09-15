# It is simple: for each x_i, you just build a regression model using PA(x_i) as the predictors, 
# then calculate the residual vector, then use hypothesis testing (t-test?) to test 
# whether or not the mean of the residual is zero. If yes, then x_i is not a rout cause variable. 
# Otherwise, it is a root cause variable

# You can add the results of the causation based T2 method in Table 1 
rm(list=ls())

library(pcalg)
require(Rgraphviz)
library(graph)

load(file = "./dat/simu/bn100")
#load(file = "./dat/real/")
load(file = "./dat/simu/shift_100")
#plot(bn.dag)

edgeMat <- edgeMatrix(bn.dag)

shift.sig <- c(0.1, 0.3, 0.5, 0.7, 1, 1.5)
mean_pop <- 0

rate <- matrix(0, ncol = 3)
for (sig in shift.sig) {
  load(file = paste("./dat/simu/dat_100_", sig, sep = ""))
  for (pos in shift.pos) {
    pos_predictors <- edgeMat["from",which(edgeMat["to",] == pos)]
    dat4reg <- data.frame(resp = dat[,pos], dat[,pos_predictors])
    fit <- lm(resp ~ ., data = dat4reg)
    residual_vec <- residuals(fit)
    rate <- rbind(rate, c(sig, pos, sum(abs(residual_vec)>0.707)/nrow(dat)))
    #test <- t.test(residual_vec)
    #pvalues <- c(pvalues, test$p.value)
  } 
}
