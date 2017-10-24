library(data.table)
library(dplyr)
library(rootSolve) 
library(xtable)

chiRoots <- function(params) {
    a <- params[1]
    b <- params[2]
    res <- 0.9 - pchisq(2/a, df=2*b) + pchisq(0.2/a, df=2*b)
    return(abs(res))
}

roots <- optim(c(20, 1), chiRoots)$par

a <- roots[1]
b <- roots[2]
q1 <- pchisq(2/a, df=2*b) 
q2 <- pchisq(0.2/a, df=2*b)

cat( "Pr(0.1 <= Gamma(", a, ", ", b, ")) = ", pgamma(1, shape=a, scale=b), '\n', sep='')

stress_df <- read.table("stress.csv")
stress <- stress_df[, -1]
lengths <- stress_df[, 1]

getPosteriorMeanAndSD <- function(ys, a, b) {
    alpha <- length(ys) + a
    beta <- b + sum(ys)
    post_mean <- alpha/beta
    post_var <- alpha/beta^2
    post_sd <- post_var %>% sqrt
    return(list(mean=post_mean, sd=post_sd))
}

mean_sd <- stress %>% apply(1, getPosteriorMeanAndSD, a=a, b=b) %>%
    sapply(rbind) %>%
    as.data.table %>%
    setNames(sapply(lengths, as.character)) %>%
    t %>% 
    as.data.table %>%
    setNames(c("Mean$(\\lambda | \\by)$", 
               "SD$(\\lambda | \\by)$"))

sink("mean_sd.tex")
print(xtable(mean_sd, digits=c(0, 3, 3)),
      sanitize.text.function=function(x) { x },
      include.rownames=FALSE)
sink()

normalRoots <- function(params) {
    mu <- params[1]
    sigma <- params[2]
    res <- 0.9 - pnorm( (log(30) - mu)/sigma ) -
            pnorm( (log(0.5) - mu)/sigma )
    return(res)
}

normal_roots <- optim(c(exp(0), exp(1)), normalRoots)$par

