library(data.table)
library(dplyr)
library(xtable)

mmEstimates <- function(ys, N, n) {
    estimates <- N*ys/n
    std_errs <- ys %>% sapply(getSE, N=N, n=n)
    return(data.table("$\\estim{X_k}_\\text{, MME}$"=estimates, 
                      "SE($\\estim{X_k}_\\text{, MME}$)"=std_errs))
}

getSE <- function(y, N, n) {
    std_err <- N^2*(N - n)/(n^2*(N - 1))*y*(1 - y/n)
    return(std_err %>% sqrt)
}

posterior <- function(ys, N, n, alphas) {
    alpha_plus <- sum(alphas)
    estimates <- (N - n)*(alphas + ys)/(alpha_plus + n) + ys
    vars <- (N - n)*(alphas + ys)*(1 - (alphas + ys)/(alpha_plus + n))*
        ((N + alpha_plus)/(alpha_plus + n + 1))/(alpha_plus + n)
    std_errs <- vars %>% sqrt
    return(data.table("$\\estim{X_k}_\\text{, Bayes}$"=estimates, 
                      "SE($\\estim{X_k}_\\text{, Bayes}$)"=std_errs))
}

N <- 800
n <- 50
ys <- c(36, 14, 0)
alphas <- rep(1, 3)

mm_df <- mmEstimates(ys, N, n)
sink("mm_df.tex")
print(xtable(mm_df,
             digits=c(0, 0, 2)),
      sanitize.text.function=function(x){x})
sink()

bayes_df <- posterior(ys, N, n, alphas)
sink("bayes_df.tex")
print(xtable(bayes_df,
             digits=c(0, 0, 2)),
      sanitize.text.function=function(x){x})
sink()
