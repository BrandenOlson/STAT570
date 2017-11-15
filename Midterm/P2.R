library(data.table)
library(dplyr)
library(rootSolve)
library(xtable)

srs <- read.table("MidtermSRS.txt")
strata_srs <- read.table("MidtermStrataSRS.txt")

ys <- srs$ysrs
zs <- srs$birthsrs

# Part (a)
mod <- lm(ys ~ zs)
beta1 <- mod$coef[2]
beta1_se <- summary(mod)[[4]][2, 2]

s_ys <- strata_srs$ystrata
s_zs <- strata_srs$birthstrata
s_xs <- strata_srs$racestrata

# Part (b)
mod_2 <- lm(s_ys ~ s_zs)
beta1_ssrs <- mod_2$coef[2]
beta1_ssrs_se <- summary(mod_2)[[4]][2, 2]

pdf("fit_b.pdf", width=10, height=6)
plot(s_ys ~ s_zs)
lines(mod_2$fitted.values ~ s_zs)
dev.off()

# Part (c)
mod_3 <- lm(s_ys ~ s_zs + s_xs)
beta1_c <- mod_3$coef[2]
beta1_c_se <- summary(mod_3)[[4]][2, 2]

# Part (d)
N <- 10000
N0 <- 2000
N1 <- 8000

s_ys_0 <- s_ys[s_xs == 0]
s_ys_1 <- s_ys[s_xs == 1]
n0 <- length(s_ys_0)
n1 <- length(s_ys_1)
ws <- ifelse(s_xs == 0, N/N0, N/N1)

mod_4 <- lm(s_ys ~ s_zs + ws)
beta1_d <- mod_4$coef[2]
beta1_d_se <- summary(mod_4)[[4]][2, 2]

# Part (e)
G <- function(beta) {
    beta <- beta %>% as.matrix
    res <- t(X) %*% W %*% (s_ys - X %*% beta)
    return(res)
}

estimG <- function(X, W, Y) {
    res <- solve(t(X) %*% W %*% X) %*% t(X) %*% W %*% Y
    return(res)
}

varBeta <- function(X, W, Y, beta) {
    sigma2 <- mean((Y - X %*% beta)^2) 
    xwx_inv <- solve(t(X) %*% W %*% X)
    res <- sigma2 * xwx_inv %*% t(X) %*% W %*% t(W) %*% X %*% t(xwx_inv)
    return(res) 
}

X <- cbind(1, s_zs)
W <- diag(ws) 
beta1_e <- estimG(X, W, s_ys)[2]
var_beta1_e <- varBeta(X, W, s_ys, betas)
beta1_e_se <- var_beta1_e %>% sqrt
