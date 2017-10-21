library(data.table)
library(dplyr)
library(rootSolve)
library(sandwich)
library(xtable)

qq_exp <- function(dat) {
    n <- length(dat)
    p <- ppoints(n)
    dat_quantile <- quantile(dat, p=p)
    plot(qexp(p), q, main="Exponential Q-Q-Plot", xlab="Theoretical quantiles", 
         ylab="Empirical quantiles")
    qqline(q, distribution=qexp, col="blue")
}

estimateAlpha <- function(ys, lambda_hat) {
    n <- length(ys)
    alpha <- sum((lambda_hat*ys - 1)^2)/(n - 2)
    return(alpha)
}

stress_dat <- read.csv("stress.csv", sep=' ', header=FALSE)
lengths <- stress_dat[, 1]
stress <- stress_dat[, -1]

names(stress) <- 1:13

estimateLambda <- function(ys) {
    return(1/mean(ys))
}

n <- ncol(stress)
lambda_hats <- apply(stress, 1, estimateLambda)
std_errs <- lambda_hats/sqrt(n)

dat_b <- data.table("$\\estim{\\lambda}$"=lambda_hats, 
                    "SE($\\estim{\\lambda}$)"=std_errs)

sink("mles.tex")
print(xtable(dat_b, digits=c(0, 3, 4)), sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()

pdf("qqplots.pdf", width=6, height=6)
par(mfrow=c(2, 2))
for(i in 1:4) {
    ys <- stress[i, ]
    pps <- ppoints(length(ys))
    q <- quantile(ys, pps)
    true_qs <- qexp(ppoints(length(ys)))
    plot(true_qs, q, 
           main=paste0("QQ plot for length = ", lengths[i]),
           xlab="Theoretical quantiles",
           ylab="Empirical quantiles"
           )
    qqline(quantile(ys, p=ppoints(length(ys))), distribution=qexp, col="blue")
}
dev.off()

pdf("histograms.pdf", width=6, height=6)
par(mfrow=c(2, 2))
for(i in 1:4) {
    hist(stress[i, ] %>% unlist, xlab="Stress", 
         main=paste0("Histogram for length = ", lengths[i]))
}
dev.off()

quasiScore <- function(lambda, ys) {
    return( ys - (1/lambda) )
}

lambdaSE <- function(lambda, ys) {
    n <- length(ys)
    std_err <- sqrt( lambda*sum((lambda*ys - 1)^2)/(n - 1) )
    return(std_err)
}


alphas <- rep(NA, 4)
for(i in 1:4) {
    ys <- stress[i, ]
    lambda_hat <- lambda_hats[i]
    alphas[i] <- estimateAlpha(ys, lambda_hat)
    std_errs[i] <- lambdaSE(lambda_hat, ys)
}

lambda_alpha_df <- data.table("Length"=lengths,
                              "$\\estim{\\lambda}$"=lambda_hats,
                              "$\\estim{\\alpha}$"=alphas,
                              "SE($\\estim{\\lambda}$)"=std_errs
                              )
sink("lambda_alpha.tex")
print(xtable(lambda_alpha_df, digits=c(0, 0, 3, 4, 4)), sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()

sand_ses <- rep(NA, 4)
for(i in 1:4) {
    ys <- stress[i, ] %>% unlist
    exp_model <- lm(log(ys) ~ 1)
    sand_ses[i] <- (sqrt(sandwich(exp_model))/length(ys))
}

getAlphaFromEta <- function(eta) {
    return( mean(ys^eta)^(1/eta) )
}

score <- function(eta, ys) {
    n <- length(ys)
    alpha <- getAlphaFromEta(eta)
    s1 <- n/eta - n*log(alpha) + sum(log(ys)) - 
        sum((ys/alpha)^eta * log(ys/alpha))
    return(s1)
}

etaAlphaSE <- function(eta, alpha, ys) {
    n <- length(ys)
    I11 <- n/eta^2 + sum( (ys/alpha)^eta * log(ys/alpha)^2 ) 
    I12 <- n/alpha - sum( (ys/alpha)^eta * (1/alpha) * ( 1 + eta*log(ys/alpha) ) )
    I22 <- eta*( sum( (1 + eta)*(ys^eta)/(alpha^(eta + 2) ) ) - (n/alpha^2) )
    info_mat <- matrix(c(I11, I12, I12, I22), 2, 2)
    inv_info_mat <- solve(info_mat)
    std_errs <- inv_info_mat %>% diag %>% sqrt
    return(std_errs)
}

eta_hats <- rep(NA, 4)
alpha_hats <- rep(NA, 4)
eta_ses <- rep(NA, 4)
alpha_ses <- rep(NA, 4)
for(i in 1:4) {
    ys <- stress[i, ]
    eta_hats[i] <- multiroot(score, 10, ys=ys)$root
    alpha_hats[i] <- getAlphaFromEta(eta_hats[i])
    std_errs <- etaAlphaSE(eta_hats[i], alpha_hats[i], ys)
    eta_ses[i] <- std_errs[1]
    alpha_ses[i] <- std_errs[2]
}

eta_alpha_df <- data.table("Length"=lengths, "$\\estim{\\eta}$"=eta_hats,
                           "SE($\\estim{\\eta}$)"=eta_ses,
                           "$\\estim{\\alpha}$"=alpha_hats,
                           "SE($\\estim{\\alpha}$)"=alpha_ses
                           )

sink("eta_alpha.tex")
print(xtable(eta_alpha_df, digits=c(0, 0, 2, 2, 2, 4)), 
      sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()
