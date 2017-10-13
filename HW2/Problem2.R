library(emdbook)
library(rootSolve)

xs <- c(6.1, 4.2, 0.5, 8.8, 1.5, 9.2, 8.5, 8.7, 6.7, 6.5, 6.3, 6.7, 0.2, 8.7, 7.5)
ys <- c(0.8, 3.5, 12.4, 1.1, 8.9, 2.4, 0.1, 0.4, 3.5, 8.3, 2.6, 1.5, 16.6, 0.1, 1.3)

info_mat <- function(x) {
    mat <- matrix(NA, 2, 2)
    mat[1, 1] <- length(x)
    mat[1, 2] <- mat[2, 1] <- sum(x)
    mat[2, 2] <- sum(x^2)
    return(mat)
}

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
    s2 <- sum(xs*(1 - exp(beta0 + beta1*xs)*ys))
    return(c(s1, s2))
}

roots <- multiroot(score, start=c(0, 0), rtol=1e-20)$root
beta0_mle <- roots[1]
beta1_mle <- roots[2]
beta0s <- seq(roots[1] - 1, beta0_mle + 1, length.out=100)
beta1s <- seq(roots[2] - 0.5, beta1_mle + 0.5, length.out=100)
inv_info_mat <- solve(info_mat(xs))
asymp_1 <- dnorm(beta0s, mean=beta0_mle, sd=sqrt(inv_info_mat[1, 1]), log=T)
asymp_2 <- dnorm(beta1s, mean=beta1_mle, sd=sqrt(inv_info_mat[2, 2]), log=T)

l1 <- l2 <- {}
for(beta0 in beta0s) {
    l1 <- c(l1, log_lik(c(beta0, beta1_mle)))
}
for(beta1 in beta1s) {
    l2 <- c(l2, log_lik(c(beta0_mle, beta1)))
}

pdf("plot1.pdf", width=10, height=6)
plot((l1 - mean(l1))/sd(l1) ~ beta0s, type="l", 
     xlab="beta0", ylab="Log-likelihood", 
     main="Log-likelihood and asymptotic approximation for beta0" )
lines((asymp_1 - mean(asymp_1))/sd(asymp_1) ~ beta0s, col="blue", lty=2)
legend("bottomleft", c("Likelihood", "Normal approximation"), 
       col=c("black", "blue"), lty=c(1, 2))
dev.off()

pdf("plot2.pdf", width=10, height=6)
plot((l2 - mean(l2))/sd(l2) ~ beta1s, type="l", ylim=c(-0.5, 1.1),
     xlab="beta0", ylab="Log-likelihood", 
     main="Log-likelihood and asymptotic approximation for beta1")
lines((asymp_2 - mean(asymp_2))/sd(asymp_2) ~ beta1s, col="blue", lty=2)
legend("topright", c("Likelihood", "Normal approximation"),
       col=c("black", "blue"), lty=c(1, 2))
dev.off()
