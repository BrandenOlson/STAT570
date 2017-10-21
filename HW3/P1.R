rm(list=ls())
library(dplyr)
library(data.table)
library(xtable)

data(Nile)

V <- function(phi, sigma, n) {
    mat <- matrix(NA, n, n)
    for(i in 1:n) {
        for(j in 1:n) {
            mat[i, j] <- sigma^2*phi^(abs(i - j))/(1 - phi^2)
        }
    }
    return(mat)
}

betaGLS <- function(y, X, V) {
    V_inv <- solve(V)
    beta_hat <- solve(t(X) %*% V_inv %*% X) %*% t(X) %*% V_inv %*% y
    return(beta_hat %>% c)
}

computeBetaCoefficients <- function(y, X) {
    betas <- solve(t(X) %*% X) %*% t(X) %*% y
    return(betas)
}

fittedValues <- function(X, beta) {
    yhat <- X %*% beta
    return(yhat)
}

meanSquaredError <- function(y, X, beta_hat) {
    n <- length(y)
    degrees_of_freedom <- ncol(X) - 1
    yhat <- fittedValues(X, beta_hat)
    mse <- ( (y - yhat)^2 %>% sum )/(n - degrees_of_freedom - 1)
    return(mse)
}

seOLS <- function(y, X) {
    beta_hat <- computeBetaCoefficients(y, X)
    s2 <- meanSquaredError(y, X, beta_hat)
    se <- sqrt(s2*diag(solve(t(X) %*% X)))
    return(se)
}

seGLS <- function(y, X, V, phi) {
	beta_hat <- betaGLS(y, X, V)
	s2 <- meanSquaredError(y, X, beta_hat)
    var_beta <- solve( t(X) %*% solve(V) %*% X )
    std_err <- sqrt(s2*diag(var_beta))
    return(std_err)
}

y <- Nile %>% c
x_bin <- c(rep(0, 28), rep(1, 72))
X <- cbind(1, x_bin)
n <- length(y)
sigma <- 1
phi <- 0.2
V_hat <- V(phi, sigma, n)



ls_mod <- lm(y ~ 1 + x_bin)
ls_beta <- ls_mod$coef
ls_fit <- ls_mod %>% fitted.values
ls_resid <- ls_mod$resid
ls_stderr <- summary(ls_mod)$sigma

pdf("lm_acf.pdf", width=10, height=6)
acf(ls_resid, main="Autocorrelation of OLS residuals")
dev.off()

gls_beta <- betaGLS(y, X, V_hat)
gls_fit <- X %*% gls_beta
gls_resid <- y - gls_fit

plot(y ~ x_bin)
lines(gls_fit ~ x_bin, col="blue")
lines(ls_fit ~ x_bin, col="red")

ls_df <- data.table("$\\estim{\\bbeta}_{\\text{OLS}}$"=ls_beta,
	"SE($\\estim{\\bbeta}_{\\text{OLS}}$)"=seOLS(y, X),
	"$\\estim{\\bbeta}_{\\text{GLS}}$"=gls_beta,
	"SE($\\estim{\\bbeta}_{\\text{GLS}}$)"=seGLS(y, X, V_hat, phi)
	)

sink("ls_df.tex")
print(xtable(ls_df), sanitize.text.function=function(x){x},
	include.rownames=FALSE)
sink()

computePhi <- function(y, X, beta) {
	n <- length(y)
	mu <- X %*% beta
	phi <- 0
	for(t in 2:n) {
		phi <- phi + (y[t] - mu[t])*(y[t-1] - mu[t-1])
	}
	phi <- phi/sum((y[1:(n-1)] - mu[1:(n-1)])^2)
	return(phi)
}
epsilon <- 1e-10
error <- Inf
beta_mle <- ls_beta
phi_mle <- computePhi(y, X, beta_mle)
V_mle <- V(phi_mle, 1, length(y))
beta_mle <- betaGLS(y, X, V_mle)
while(error > epsilon) {
	beta_old <- beta_mle
	phi_mle <- computePhi(y, X, beta_mle)
	V_mle <- V(phi_mle, 1, length(y))
	beta_mle <- betaGLS(y, X, V_mle)	
	error <- norm( as.matrix(beta_old) - as.matrix(beta_mle) )
}

beta_mle_df <- data.table("$\\estim{\\beta}_\\text{MLE}$"=beta_mle,
	"SE($\\estim{\\beta}_\\text{MLE}$)"=seGLS(y, X, V_mle, phi_mle))
sink("beta_mle.tex")
print(xtable(beta_mle_df), sanitize.text.function=function(x){x},
	include.rownames=FALSE)
sink()
	
