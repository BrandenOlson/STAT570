library(reshape2)
library(dplyr)
library(ggplot2)
library(sn)
library(xtable)

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

getIntervalWidths <- function(y, X) {
	alpha <- 0.05
	std_errs <- coefficientStandardErrors(y, X)
	widths <- qnorm(1 - alpha/2)*std_errs
	return(widths)
}
 
simulateData <- function(n, type) {
    if(type == "normal") {
        eps <- rnorm(n, 0, 3)
    } else if(type == "uniform") {
        eps <- runif(n, -3, 3)
    } else if(type == "skewnormal") {
        omega <- 1
        alpha <- 5
        delta <- alpha/sqrt(1 + alpha^2)
        xi <- -omega*delta*sqrt(2/pi)
        eps <- rsn(n=n, xi=xi, omega=omega, alpha=alpha)
    }

    ys <- beta_0 + beta_1*xs + eps
    return(ys)
}

OLSVariance <- function(type, X) {
    var_eps <- {}
    if(type == "normal") {
        var_eps <- 9
    } else if(type == "uniform") {
        var_eps <- 3
    } else if(type == "skewnormal") {
        var_eps <- 1 - 25/(13*pi)
    }
    var_matrix <- var_eps * solve(t(X) %*% X)
    return(var_matrix %>% diag %>% unname)
}

beta_0 <- 2
beta_1 <- -2
beta <- c(beta_0, beta_1)

x_mean <- 20
x_sd <- 3

ns <- c(10, 25)
trial.count <- 1000
i <- 1
beta.norms <- list()
beta.unifs <- list()
beta.skews <- list()
distribution_names <- c("Normal", "Uniform", "Skew-normal")
n_var_df <- list()
coverage_df <- data.frame(matrix(ncol=3, nrow=0)) %>% 
    setNames(c("Coefficient", "Distribution", "Coverage"))

for(n in ns) {
    xs <- rnorm(n, x_mean, x_sd)
    X <- cbind(1, xs)
    sim.df <- matrix(NA, n, 3) %>% data.frame
    true_normal_var <- OLSVariance("normal", X)
    true_uniform_var <- OLSVariance("uniform", X)
    true_skew_var <- OLSVariance("skewnormal", X)

    norm.covers <- {}
    unif.covers <- {}
    skew.covers <- {}

	beta.norms[[i]] <- beta.unifs[[i]] <- beta.skews[[i]] <- matrix(NA, trial.count, 2)
	for(trial in 1:trial.count) {
    	y.norm <- simulateData(n, type="normal")
    	y.unif <- simulateData(n, type="uniform")
    	y.skew <- simulateData(n, type="skewnormal")

		beta.norms[[i]][trial, ] <- beta.norm <- computeBetaCoefficients(y.norm, X)
		beta.unifs[[i]][trial, ] <- beta.unif <- computeBetaCoefficients(y.unif, X)
		beta.skews[[i]][trial, ] <- beta.skew <- computeBetaCoefficients(y.skew, X)

        norm.widths <- getIntervalWidths(y.norm, X)
        unif.widths <- getIntervalWidths(y.unif, X)
        skew.widths <- getIntervalWidths(y.skew, X)

        norm.covers.trial <- beta.norm - norm.widths < beta &
            beta < beta.norm + norm.widths
        unif.covers.trial <- beta.unif - unif.widths < beta_0 &
            beta < beta.unif + unif.widths
        skew.covers.trial <- beta.skew - skew.widths < beta_0 &
            beta < beta.skew + skew.widths

        # norm.covers.0 <- c(norm.covers.0, norm.covers.trial[1])
        # norm.covers.1 <- c(norm.covers.1, norm.covers.trial[2])

	}
	
	sample_normal_var <- beta.norms[[i]] %>% apply(2, var)
	sample_uniform_var <- beta.unifs[[i]] %>% apply(2, var)
	sample_skew_var <- beta.skews[[i]] %>% apply(2, var)

    
    print(norm.covers)
    norm.coverage <- mean(norm.covers)
    unif.coverage <- mean(unif.covers)
    skew.coverage <- mean(skew.covers)

	variance_mat <- rbind(true_normal_var, sample_normal_var,
		true_uniform_var, sample_uniform_var,
		true_skew_var, sample_skew_var)
	n_var_df[[i]] <- data.frame(variance_mat, 
		Type=c("Theoretical", "Empirical"),
		Distribution=rep(distribution_names, each=2),
	    n)

    i <- i + 1
}

var_df <- n_var_df %>% ldply(data.frame)
names(var_df)[c(1:2, 5)] <- c("$\\Var{\\bbeta_0}$", 
						"$\\Var{\\bbeta_1}$",
						"$n$")

sink("var_df.tex")
print(xtable(var_df, digits=c(0, 2, 3, 0, 0, 0)), sanitize.text.function=function(x){x})
sink()

df.melt <- sim.data.frames %>% first %>% melt
p <- df.melt %>% ggplot + geom_density(aes(x=value, fill=variable), alpha=0.5) + 
        xlim(-60, -10)

