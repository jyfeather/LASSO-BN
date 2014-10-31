solver <- function(wei, yt, df, type = c("LASSO-BN", "VS-MSPC")) {
  num <- dim(wei)[1]
  y <- solve(wei) %*% yt
  X <- diag(num)
  if (type == "LASSO-BN") {
    D <- diag(num)
  } else if (type == "VS-MSPC") {
    D <- wei
  } else {
    stop("Wrong Type!")
  }
  
  fit <- genlasso(y, X, D, maxsteps = 50) # our df is always less than 50
  
  beta <- coef.genlasso(fit, df = df)$beta
  beta <- abs(beta)
  posi <- sort(beta, decreasing = T, index.return = T)$ix
  coef <- rep(0, num)
  coef[posi[1:df]] = 1

  obj <- summary(fit)[df,"rss"]
  return(list(obj.val = obj, coef.set = coef))
}