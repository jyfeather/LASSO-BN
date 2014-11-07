rm(list = ls())

#setwd("C:\\Users\\Administrator\\Dropbox\\Research\\Monitoring\\code")
setwd("C:\\Users\\JIN\\Dropbox\\Research\\Monitoring\\code")

y <- c(3,1,2) # min{1/2||y - beta||_2^2 + lambda*||beta||_1}

# genlasso
require("genlasso") 
out.gen = genlasso(y, diag(3), diag(3))
coef(out.gen, df=1)$beta

# lars
require("lars")
out.lars = lars(diag(3), y)
coef(out.lars)[1+1,]

# lassoshooting
require("lassoshooting")
out.shoot = lassoshooting(diag(3), y, 2) # Lambda need be specified
beta = out.shoot$coefficients

# glmnet
require("glmnet")
out.glm = glmnet(diag(3), y, alpha = 1, intercept = F)
coef(out.glm)

# Compute objective value
obj = 0.5 * sum((y - beta)^2) + lambda * sum(abs(beta))