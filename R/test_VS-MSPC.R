#############################################
# a test code for VS-MSPC
#############################################
rm(list=ls())

#############################################
# simulate data
# dim: 1000 observations X 100 features
#      where 800 observations are normal;
#      200 observations are abnormal 
#############################################
library(MASS)
num_var <- 100
num_obs_normal <- 800
num_obs_abnormal <- 200
dat_normal <- mvrnorm(n = num_obs_normal, mu = rep(0, num_var), Sigma = diag(num_var))
# 1st feature is out of control
dat_abnormal <- mvrnorm(n = num_obs_abnormal, mu = c(1, rep(0, num_var-1)), Sigma = diag(num_var)) 
dat <- rbind(dat_normal, dat_abnormal)

#############################################
# solve Equ.(9) to estimate mu 
# with forward variable selection method
# entering rule: F-to-enter rule
#############################################
X <- chol(solve(cov(dat)))
s <- 1 # number of nonzero variables
var_set <- c() # forward variable selection
t <- 1
while(s > 0) {
  Y <- dat[t,]; Y <- X%*%Y 
  fit <- lm(Y ~ X-1)
  R_sq <- summary(fit)$r.squared
  # compute R squared for k-1 features
  no_enter <- 0
  F_stat <- 0
  for(i in 1:num_var) {
    X2 <- X[,-i]
    fit2 <- lm(Y ~ X2-1)
    R_sq2 <- summary(fit2)$r.squared 
    # compute F
    F <- (R_sq^2-R_sq2^2)/(1-R_sq^2)
    # add the feature with largest F into the var_set
  }
  s <- s-1
}

#############################################
# alarm chart
# Equ. 10
#############################################
