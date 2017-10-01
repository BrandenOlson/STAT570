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

reg_summary <- cbind(beta_hat, se, t_value, p_value)
colnames(reg_summary) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
row.names(reg_summary) <- c("(Intercept)", "lcavol", "svi", "lcavol:svi")

print(reg_summary)
