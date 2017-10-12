library(emdbook)
library(rootSolve)

xs <- c(6.1, 4.2, 0.5, 8.8, 1.5, 9.2, 8.5, 8.7, 6.7, 6.5, 6.3, 6.7, 0.2, 8.7, 7.5)
ys <- c(0.8, 3.5, 12.4, 1.1, 8.9, 2.4, 0.1, 0.4, 3.5, 8.3, 2.6, 1.5, 16.6, 0.1, 1.3)

log_lik <- function(betas) {
    beta0 <- betas[1]
    beta1 <- betas[2]

    n <- length(xs)
    ell <- n*beta0 + sum(beta1*xs - exp(beta0 + beta1*xs)*ys)
    return(ell)
}

score <- function(betas) {
    beta0 <- betas[1]
    beta1 <- betas[2]
    n <- length(xs)
    s1 <- n - sum(exp(beta0 + beta1*xs)*ys)
    s2 <- sum(xs*(1 - beta1*exp(beta0 + beta1*xs)*ys))
    return(c(s1, s2))
}

roots <- multiroot(score, start=c(0, 0), rtol=1e-20)$root
beta0s <- seq(roots[1] - 1, roots[1] + 1, length.out=100)
beta1s <- seq(roots[2] - 1, roots[2] + 1, length.out=100)

reg <- function(beta0, beta1) {

}
