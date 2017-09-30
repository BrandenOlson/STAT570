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

fittedValues <- function(y, X) {
    beta <- computeBetaCoefficients(y, X)
    return( X %*% beta )
}

meanSquaredError <- function(y, yhat) {
    return( (y - yhat)^2 %>% mean )
}

# Not in agreement yet...
coefficientStandardErrors <- function(y, X) {
    yhat <- fittedValues(y, X)
    s2 <- meanSquaredError(y, yhat)
    se <- (s2*solve( t(X) %*% X )) %>% diag %>% sqrt
    return(se)
}

beta <- computeBetaCoefficients(y, X)
se <- coefficientStandardErrors(y, X)
yhat <- fittedValues(y, X)

print(beta)
print(se)
