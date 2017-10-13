library(dplyr)
library(MASS)
library(sandwich)
library(xtable)

beta_0 <- -2
beta_1 <- log(2)
bbeta <- c(beta_0, beta_1)
bs <- c(0.2, 1, 10, 1000)
ns <- c(10, 20, 50, 100, 250)
alpha <- 0.95
trial.count <- 1000

cov_df_names <- c("Parameter", "b", "n", "Model", "Coverage")
coverage_df <- matrix(ncol=length(cov_df_names), nrow=0) %>% data.frame %>% setNames(cov_df_names)
for(b in bs) {
    for(n in ns) {
        cat("b =", b, "n =", n, '\n')
        lik_covers_0 <- rep(NA, trial.count)
        lik_covers_1 <- rep(NA, trial.count)
        quasi_covers_0 <- rep(NA, trial.count)
        quasi_covers_1 <- rep(NA, trial.count)
        sandwich_covers_0 <- rep(NA, trial.count)
        sandwich_covers_1 <- rep(NA, trial.count)
        for(trial in 1:trial.count) {
            xs <- rnorm(n, mean=0, sd=1)
            mus <- exp(beta_0 + beta_1*xs)

            thetas <- rgamma(n, shape=b, rate=b)
            ys <- rpois(n, mus*thetas)

            lik_covers <- tryCatch({
                lik_mod <- glm(ys ~ 1 + xs, family="poisson")
                lik_ci <- suppressMessages(confint(lik_mod, level=alpha))
                lik_covers <- lik_ci[, 1] < bbeta & bbeta < lik_ci[, 2]
            }, error = function(e) { lik_covers <- c(0, 0) }
            )

            lik_covers_0[trial] <- ifelse(!is.na(lik_covers[1]), lik_covers[1], 0)
            lik_covers_1[trial] <- ifelse(!is.na(lik_covers[2]), lik_covers[2], 0)

            quasi_covers <- tryCatch({
                quasi_mod <- glm(ys ~ 1 + xs, family="quasipoisson")
                quasi_ci <- suppressMessages(confint(quasi_mod, level=alpha))
                quasi_covers <- quasi_ci[, 1] < bbeta & bbeta < quasi_ci[, 2]
            }, error = function(e) { quasi_covers <- c(0, 0) }
            )
            quasi_covers_0[trial] <- quasi_covers[1]
            quasi_covers_1[trial] <- quasi_covers[2]

            sandwich_mod <- sandwich(lik_mod)
            sandwich_se <- sandwich_mod %>% diag %>% sqrt
            sandwich_left <- lik_mod$coef - 1.96*sandwich_se
            sandwich_right <- lik_mod$coef + 1.96*sandwich_se
            sandwich_covers <- sandwich_left < bbeta & bbeta < sandwich_right
            sandwich_covers_0[trial] <- sandwich_covers[1]
            sandwich_covers_1[trial] <- sandwich_covers[2]
        }
        coverage_df <- rbind(coverage_df,
                             data.frame(Parameter="$\\beta_0$", n=n, b=b, 
                                        Model="Poisson", Coverage=mean(lik_covers_0)),
                             data.frame(Parameter="$\\beta_1$", n=n, b=b, 
                                        Model="Poisson", Coverage=mean(lik_covers_1)),
                             data.frame(Parameter="$\\beta_0$", n=n, b=b, 
                                        Model="Quasi-likelihood", 
                                        Coverage=mean(quasi_covers_0)),
                             data.frame(Parameter="$\\beta_1$", n=n, b=b, 
                                        Model="Quasi-likelihood", 
                                        Coverage=mean(quasi_covers_1)),
                             data.frame(Parameter="$\\beta_0$", n=n, b=b, 
                                        Model="Sandwich", 
                                        Coverage=mean(sandwich_covers_0)),
                             data.frame(Parameter="$\\beta_1$", n=n, b=b, 
                                        Model="Sandwich", 
                                        Coverage=mean(sandwich_covers_1))
                             )
    }
}

beta0_df <- coverage_df[coverage_df$Parameter == "$\\beta_0$", names(coverage_df) != "Parameter"]
beta1_df <- coverage_df[coverage_df$Parameter == "$\\beta_1$", names(coverage_df) != "Parameter"]
beta0_df_lik <- beta0_df[beta0_df$Model == "Poisson", names(beta0_df) != "Model"]
beta0_df_quasi <- beta0_df[beta0_df$Model == "Quasi-likelihood", 
                           names(beta0_df) != "Model"]
beta0_df_sand <- beta0_df[beta0_df$Model == "Sandwich", names(beta0_df) != "Model"]
beta1_df_lik <- beta1_df[beta1_df$Model == "Poisson", names(beta1_df) != "Model"]
beta1_df_quasi <- beta1_df[beta1_df$Model == "Quasi-likelihood", 
                           names(beta1_df) != "Model"]
beta1_df_sand <- beta1_df[beta1_df$Model == "Sandwich", names(beta1_df) != "Model"]

sink("beta0_lik.tex")
print(xtable(beta0_df_lik), digits=c(0, 1, 0, 3), sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()

sink("beta1_lik.tex")
print(xtable(beta1_df_lik), digits=c(0, 1, 0, 3), sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()

sink("beta0_quasi.tex")
print(xtable(beta0_df_quasi), digits=c(0, 1, 0, 3), sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()

sink("beta1_quasi.tex")
print(xtable(beta1_df_quasi), digits=c(0, 1, 0, 3), sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()

sink("beta0_sand.tex")
print(xtable(beta0_df_sand), digits=c(0, 1, 0, 3),
      include.rownames=FALSE ,
      sanitize.text.function=function(x){x})
sink()

sink("beta1_sand.tex")
print(xtable(beta1_df_sand), digits=c(0, 1, 0, 3), 
      include.rownames=FALSE,
      sanitize.text.function=function(x){x})
sink()
