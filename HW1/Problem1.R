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

pValue <- function(t_value, n, p) {
    p_val <- 2*(1 - pt(t_value, df=(n - p - 1)))
    return(p_val)
}


beta_hat <- computeBetaCoefficients(y, X)
se <- coefficientStandardErrors(y, X)
yhat <- fittedValues(X, beta_hat)
t_value <- tValue(beta_hat, se)
p_value <- pValue(t_value, n=length(y), p=(ncol(X) - 1))

options(digits=4)
print(beta_hat %>% signif(digits=5))
print(se %>% signif(digits=5))
options(digits=3)
print(t_value %>% signif(digits=5))
print(p_value %>% signif(digits=3))
