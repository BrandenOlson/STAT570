library(data.table)
library(dplyr)
library(lasso2)
library(MASS)
library(matrixcalc)
library(xtable)

data(Prostate)

y <- Prostate$lpsa
x <- Prostate$lcavol
X <- cbind(1, x)
n <- length(y)

# Hyperparameter specification
m0 <- m1 <- 0
v00 <- v11 <- 2
v01 <- 0
a <- b <- 0

m <- c(m0, m1)
V <- matrix(c(v00, v01, v01, v11), 2, 2)

mle_model <- lm(y ~ x)
beta_mle <- mle_model$coef %>% unname
var_mle <- mle_model %>% vcov

RSS <- function(y, yhat) {
    res <- sum((y - yhat)^2)/(n - 1)
    return(res)
}

samplePiBeta <- function() {
    samp <- MASS::mvrnorm(n=1, mu=m, Sigma=V)    
    return(samp)
}

sampleBetas <- function(sigma2) {
    k <- ncol(X)
    W <- getW(sigma2)
    m_star <- W %*% beta_mle + (diag(rep(1, k)) - W) %*% m
    V_star <- W %*% var_mle
    samp <- MASS::mvrnorm(n=1, mu=m_star, Sigma=V_star)
    return(samp)
}

sampleSigma2 <- function(beta) {
    samp <- rgamma(n=1, a + n/2, b + crossprod(y - X %*% beta)/2)
    return(1/samp)
}

# Parameter specification
getW <- function(sigma2) {
    res <- solve((t(X) %*% X) + solve(V) * sigma2) %*% ( t(X) %*% X )
    return(res)
}

doMCMC <- function(beta_init, sigma2_init) {
    beta <- beta_init
    sigma2 <- sigma2_init
    trial_count <- 10000
    sample_mat <- matrix(NA, nrow=trial_count, ncol=3)
    for(trial in 1:trial_count) {
        beta_prev <- beta
        sigma2_prev <- sigma2
        beta <- sampleBetas(sigma2_prev)
        sigma2 <- sampleSigma2(beta_prev)
        sample_mat[trial, ] <- c(beta, sigma2)
    }
    burn_in <- trial.count/10
    return(sample_mat[(burn_in + 1):trial_count, ])
}

mcmc_1 <- doMCMC(beta_mle, RSS(y, mle_model$fit))
mcmc_2 <- doMCMC(samplePiBeta(), 100)

pdf("hists_1.pdf", width=10, height=14)
par(mfrow=c(3, 1))
hist(mcmc_1[, 1], xlab="beta0", main="Histogram of beta0")
hist(mcmc_1[, 2], xlab="beta1", main="Histogram of beta1")
hist(mcmc_1[, 3], xlab="sigma2", main="Histogram of sigma2")
dev.off()

pdf("scat_1.pdf", width=10, height=14)
par(mfrow=c(3, 1))
plot(mcmc_1[, 1] ~ mcmc_1[, 2], xlab="beta1", ylab="beta0")
plot(mcmc_1[, 1] ~ mcmc_1[, 3], xlab="sigma2", ylab="beta0")
plot(mcmc_1[, 2] ~ mcmc_1[, 3], xlab="sigma2", ylab="beta1")
dev.off()


pdf("trace_1.pdf", width=10, height=10)
par(mfrow=c(3, 1))
plot(mcmc_1[, 1], xlab="t", ylab="beta0(t)")
plot(mcmc_1[, 2], xlab="t", ylab="beta1(t)")
plot(log(mcmc_1[, 3]), xlab="t", ylab="log(sigma^2(t))")
dev.off()


pdf("hists_2.pdf", width=10, height=14)
par(mfrow=c(3, 1))
hist(mcmc_2[, 1], xlab="beta0", main="Histogram of beta0")
hist(mcmc_2[, 2], xlab="beta1", main="Histogram of beta1")
hist(mcmc_2[, 3], xlab="sigma2", main="Histogram of sigma2")
dev.off()

pdf("scat_2.pdf", width=10, height=14)
par(mfrow=c(3, 1))
plot(mcmc_2[, 1] ~ mcmc_2[, 2], xlab="beta1", ylab="beta0")
plot(mcmc_2[, 1] ~ mcmc_2[, 3], xlab="sigma2", ylab="beta0")
plot(mcmc_2[, 2] ~ mcmc_2[, 3], xlab="sigma2", ylab="beta1")
dev.off()

pdf("trace_2.pdf", width=10, height=10)
par(mfrow=c(3, 1))
plot(mcmc_2[, 1], xlab="t", ylab="beta0(t)")
plot(mcmc_2[, 2], xlab="t", ylab="beta1(t)")
plot(log(mcmc_2[, 3]), xlab="t", ylab="log(sigma^2(t))")
dev.off()

mcmc_dat_1 <- mcmc_1 %>% 
    apply(2, quantile, probs=c(0.1, 0.5, 0.9)) %>% 
    t %>% 
    data.table

sink("mcmc_1.tex")
print(xtable(mcmc_dat_1, digits=c(0, rep(3, 3))),
      include.rownames=FALSE)
sink()


mcmc_dat_2 <- mcmc_2 %>% 
    apply(2, quantile, probs=c(0.1, 0.5, 0.9)) %>% 
    t %>% 
    data.table

sink("mcmc_2.tex")
print(xtable(mcmc_dat_2, digits=c(0, rep(3, 3))),
      include.rownames=FALSE)
sink()

probGreaterThan05 <- function(beta) {
    beta_true <- beta[beta > 0.5]
    return( length(beta_true)/length(beta) )
}

p1 <- probGreaterThan05(mcmc_1[, 2])
p2 <- probGreaterThan05(mcmc_2[, 2])
