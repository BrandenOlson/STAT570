library(dplyr)

srs <- read.table("MidtermSRS.txt")
strata_srs <- read.table("MidtermStrataSRS.txt")

ys <- srs$ysrs
zs <- srs$birthsrs

mu <- mean(ys)
mu_var <- var(ys)/length(ys)

z_true <- 2003.898
mu_rat <- z_true*mean(ys)/mean(zs)
mu_rat_var <- z_true^2 * var(ys) / var(zs)

mod <- lm(ys ~ zs)
beta1_hat <- mod$coef[2]
mu_reg <- mean(ys) + beta1_hat*(mean(zs) - z_true)

# Problem 2

beta1_se <- summary(mod)[[4]][2, 2]

strata_ys <- strata_srs$ystrata
strata_zs <- strata_srs$birthstrata
strata_xs <- strata_srs$racestrata

mod_2 <- lm(strata_ys ~ strata_zs)
mod_3 <- lm(strata_ys ~ strata_zs+strata_xs)
