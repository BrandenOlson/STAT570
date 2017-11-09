library(dplyr)

srs <- read.table("MidtermSRS.txt")
strata_srs <- read.table("MidtermStrataSRS.txt")

ys <- srs$ysrs
zs <- srs$birthsrs
n <- length(ys)

mu_y <- mean(ys)
mu_z <- mean(zs)
mu_var <- var(ys)/length(ys)

z_true <- 2003.898
mu_rat <- z_true*mu_y / mu_z
mu_rat_var <- z_true^2 * (mu_y/mu_z)^2 * ( var(ys)/(n*mu_y^2) - 
                                          2*cov(ys, zs)/(n*mu_y*mu_z) +
                               var(zs)/(n*mu_z^2) )
 # mu_rat_var <- z_true^2 * sum((ys/zs - mean(ys/zs))^2)/(n - 1)

mod <- lm(ys ~ zs)
beta1_hat <- mod$coef[2]
mu_reg <- mean(ys) + beta1_hat*(mean(zs) - z_true)
Z <- cbind(1, zs)
mu_reg_var <- (var(ys) + vcov(mod)[2,2] * var(zs))/n

# Problem 2

beta1_se <- summary(mod)[[4]][2, 2]

strata_ys <- strata_srs$ystrata
strata_zs <- strata_srs$birthstrata
strata_xs <- strata_srs$racestrata

mod_2 <- lm(strata_ys ~ strata_zs)
mod_3 <- lm(strata_ys ~ strata_zs+strata_xs)
