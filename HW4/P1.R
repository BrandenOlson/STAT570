library(data.table)
library(dplyr)
library(invgamma)
library(nleqslv)
library(rstan)
library(xtable)

gammaRoots <- function(params) {
    a <- params[1]
    b <- params[2]
    f1 <- 0.05 - pgamma(0.1, a, b)
    f2 <- 0.95 - pgamma(1, a, b)
    return(c(f1, f2))
}

roots <- nleqslv(c(3, 2), gammaRoots)$x

a <- roots[1]
b <- roots[2]

cat("a =", a, ", b =", b, '\n')


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

plotHistogram <- function(ys, a, b, length, type) {
    n <- length(ys)
    sample_size <- 10000
    if(type == "gamma") {
        lambdas <- rgamma(n=sample_size, n + a, b + sum(ys))
    } else if(type == "invgamma") {
        lambdas <- rinvgamma(n=sample_size, n + a, 1/(b + sum(ys)))
    }
    hist(lambdas, breaks=30, main=paste0("Posterior for length = ", length))
}

pdf("histograms.pdf", width=10, height=10)
par(mfrow=c(2, 2))
for(i in 1:4) {
    ys <- stress[i, ] 
    plotHistogram(ys, a, b, lengths[i], "gamma")
}
dev.off()

pdf("inv_histograms.pdf", width=10, height=10)
par(mfrow=c(2, 2))
for(i in 1:4) {
    ys <- stress[i, ] 
    plotHistogram(ys, a, b, lengths[i], "invgamma")
}
dev.off()

normalRoots <- function(params, left_prob, right_prob) {
    mu <- params[1]
    sigma <- params[2]
    f1 <- 0.05 - pnorm( (log(left_prob) - mu)/sigma )
    f2 <- 0.95 - pnorm( (log(right_prob) - mu)/sigma )
    return(c(f1, f2))
}

eta_roots <- nleqslv(c(1, 1), normalRoots, left_prob=0.5, right_prob=30)$x
alpha_roots <- nleqslv(c(1, 1), normalRoots, left_prob=1, right_prob=4)$x

