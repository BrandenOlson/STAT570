library(data.table)
library(dplyr)
library(rootSolve)
library(sandwich)
library(xtable)

ts <- c(94.3, 15.7, 62.9, 126, 5.2, 31.4, 1.1, 1.1, 2.1, 10.5)
ys <- c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22)
n <- length(ys)

# Problem 1
lambda_a <- sum(ys)/sum(ts)
se_a <- sqrt( lambda_a/(sum(ts)) )

dat_a <- data.table("$\\estim{\\lambda}$"=lambda_a, 
                    "SE($\\estim{\\lambda}$)"=se_a,
                    Estimator="MLE")
sink("a.tex")
print(xtable(dat_a, digits=c(0, 5, 5, 0)),
      include.rownames=FALSE,
      sanitize.text.function=function(x){x}
      )
sink()

# Problem 2
V_inv_tilde <- diag(1/ts)
lambda_b <- (t(ts) %*% V_inv_tilde %*% ys)/(t(ts) %*% V_inv_tilde %*% ts)
lambda_b <- lambda_b[1, 1]
mus <- lambda_b*ts
alpha <- sum( (ys - mus)^2/mus  )/(n - 1)
V <- lambda_b*diag(ts)
V_inv <- V_inv_tilde/lambda_b

se_b <- 1/sqrt( ( t(ts) %*% V_inv %*% ts )/alpha )

dat_b <- data.table("$\\estim{\\lambda}$"=lambda_b, 
                    "SE($\\estim{\\lambda}$)"=se_b,
                    Estimator="Quasi-likelihood")
sink("b.tex")
print(xtable(dat_b, digits=c(0, 5, 5, 0)),
      include.rownames=FALSE,
      sanitize.text.function=function(x){x}
      )
sink()

# Problem 3

A <- -sum(ys/lambda_a^2)/n
B <- sum((ys/lambda_a - ts)^2)/n
se_c <- (B/(n*A^2)) %>% sqrt

dat_c <- data.table("$\\estim{\\lambda}$"=lambda_a,
                    "SE($\\estim{\\lambda}$)"=se_c,
                    Estimator="Sandwich")

sink("c.tex")
print(xtable(dat_c, digits=c(0, 5, 5, 0)),
      include.rownames=FALSE,
      sanitize.text.function=function(x){x})
sink()

# Problem 4

log_lik <- function(params) {
    alpha <- params[1]
    beta <- params[2]
    
    ell_1 <- sum(log(beta) + digamma(alpha + ys) - digamma(alpha) - log(beta + ts))
    ell_2 <- sum(alpha/beta - (alpha + ys)/(beta + ts))
    
    return(c(ell_1, ell_2))
}

roots <- multiroot(log_lik, c(1, 1))$root
alpha_hat <- roots[1]
beta_hat <- roots[2]
ab_dat <- data.table("$\\estim{\\alpha}$"=alpha_hat,
                     "$\\estim{\\beta}$"=beta_hat)
sink("alphabeta.tex")
print(xtable(ab_dat, digits=c(0, 3, 3)),
      include.rownames=FALSE,
      sanitize.text.function=function(x){x})
sink()

sampleThetaI <- function(i, alpha, beta) {
    theta_i <- rgamma(n=1, ys[i] + alpha, ts[i] + beta)
    return(theta_i)
}

trial_count <- 1000
theta_dat <- matrix(NA, nrow=trial_count, ncol=n) %>% data.frame
for(trial in 1:trial_count) {
    sample_theta <- rep(NA, n)
    for(i in 1:n) {
        sample_theta[i] <- sampleThetaI(i, alpha_hat, beta_hat)
    }
    theta_dat[trial, ] <- sample_theta
}

q_1 <- quantile(theta_dat[, 1], probs=c(0.05, 0.5, 0.95))
q_10 <- quantile(theta_dat[, 10], probs=c(0.05, 0.5, 0.95))
q_dat <- data.table(rbind(q_1, q_10))
names(q_dat) <- c("5\\%", "50\\%", "95\\%")
q_dat$Estimator <- c("$\\estim{\\theta_1}$", "$\\estim{\\theta_{10}}$")

sink("q.tex")
print(xtable(q_dat, digits=c(0, rep(3, 3), 0)),
      include.rownames=FALSE,
      sanitize.text.function=function(x){x})
sink()


sampleAlpha <- function() {
    alpha <- rexp(n=1, 1)
    return(alpha)
}

sampleBeta <- function() {
    beta <- rgamma(n=1, 0.1, 0.1)
    return(beta)
}

sampleTheta <- function(alpha, beta) {
    theta <- rgamma(n=1, alpha, beta)
    return(theta)
}
