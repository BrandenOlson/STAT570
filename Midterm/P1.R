library(data.table)
library(dplyr)
library(xtable)

srs <- read.table("MidtermSRS.txt")
strata_srs <- read.table("MidtermStrataSRS.txt")

ys <- srs$ysrs
zs <- srs$birthsrs
n <- length(ys)

mu_y <- mean(ys)
mu_z <- mean(zs)
mu_var <- var(ys)/length(ys)
mu_se <- mu_var %>% sqrt

z_true <- 2003.898
mu_rat <- z_true*mu_y / mu_z

mu_rat_var <- (z_true^2 / (n * mu_z^2))*(var(ys) - 2*mu_y*cov(ys, zs)/mu_z +
                                         mu_y^2*var(zs)/mu_z^2) 
mu_rat_se <- mu_rat_var %>% sqrt

mod <- lm(ys ~ zs)
beta1_hat <- mod$coef[2]
mu_reg <- mean(ys) + beta1_hat*(mean(zs) - z_true)
Z <- cbind(1, zs)
mu_reg_var <- (var(ys) - 2*beta1_hat*cov(ys, zs) + beta1_hat^2 * var(zs))/n
mu_reg_se <- mu_reg_var %>% sqrt

mu_dat <- data.table("$\\estim{\\E{\\mu}}$"=c(mu_y, mu_rat, mu_reg),
                     "$\\SE{\\mu}$"=c(mu_se, mu_rat_se, mu_reg_se),
                     Estimator=c("$\\estim{\\mu}$",
                                 "$\\estim{\\mu}_\\text{rat}$",
                                 "$\\estim{\\mu}_\\text{reg}$"))

sink("mus.tex")
print(xtable(mu_dat, digits=c(0, 2, 4, 0)),
      sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()

# Part (g)
N <- 10000
N0 <- 2000
N1 <- 8000
s_ys <- strata_srs$ystrata
s_zs <- strata_srs$birthstrata
s_xs <- strata_srs$racestrata
s_ys_0 <- s_ys[s_xs == 0]
s_ys_1 <- s_ys[s_xs == 1]
n0 <- length(s_ys_0)
n1 <- length(s_ys_1)
ws <- ifelse(s_xs == 0, N/N0, N/N1)

mu_ssrs <- mean(s_ys)
mu_ssrs_var <- var(s_ys)/length(s_ys)
mu_ssrs_se <- mu_ssrs_var %>% sqrt


mu_w <- sum(ws*s_ys)/sum(ws)
mu_w_var <- (n0*N1^2*var(s_ys_0) + n1*N0^2*var(s_ys_1))/(n0*N1 + n1*N0)^2
mu_w_se <- mu_w_var %>% sqrt

mu_dat_2 <- data.table("$\\estim{\\E{\\mu}}$"=c(mu_ssrs, mu_w),
                       "$\\SE{\\mu}$"=c(mu_ssrs_se, mu_w_se),
                       Estimator=c("$\\estim{\\mu}$",
                                   "$\\estim{\\mu}_w$"))

sink("mu2.tex")
print(xtable(mu_dat_2, digits=c(0, 2, 4, 0)),
      sanitize.text.function=function(x){x},
      include.rownames=FALSE)
sink()

