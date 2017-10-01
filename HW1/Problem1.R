rm(list=ls())
library(lasso2)
library(dplyr)

data(Prostate)
lpsa <- Prostate$lpsa
svi <- Prostate$svi
lcavol <- Prostate$lcavol

y <- lpsa
X <- cbind(1, lcavol, svi, lcavol*svi)

computeBetaCoefficients <- function(y, X) {
    betas <- solve(t(X) %*% X) %*% t(X) %*% y
    return(betas)
}

fittedValues <- function(X, beta) {
    yhat <- X %*% beta
    return(yhat)
}

meanSquaredError <- function(y, X) {
    n <- length(y)
    degrees_of_freedom <- ncol(X) - 1
    beta_hat <- computeBetaCoefficients(y, X)
    yhat <- fittedValues(X, beta_hat)
    mse <- ( (y - yhat)^2 %>% sum )/(n - degrees_of_freedom - 1)
    return(mse)
}

coefficientStandardErrors <- function(y, X) {
    beta_hat <- computeBetaCoefficients(y, X)
    yhat <- fittedValues(X, beta_hat)
    s2 <- meanSquaredError(y, X)
    se <- sqrt(s2*diag(solve(t(X) %*% X)))
    return(se)
}

tValue <- function(betas, se) {
    return(betas/se)
}

beta_hat <- computeBetaCoefficients(y, X)
se <- coefficientStandardErrors(y, X)
yhat <- fittedValues(X, beta_hat)
t_value <- tValue(beta_hat, se)

print(beta_hat %>% signif(5))
print(se %>% signif(5))
print(t_value %>% signif(5))
