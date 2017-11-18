library(data.table)
library(dplyr)
library(xtable)

hd <- c(396, 568, 1212, 171, 554, 1104, 257, 435, 295, 397, 288, 1004, 431, 
        795, 1621, 1378, 902, 958, 1283, 2415)
non_hd <- c(375, 375, 752, 208, 151, 116, 736, 192, 315, 1252, 675, 700, 440,
            771, 688, 426, 410, 979, 377, 503)

pdf("cdfs.pdf", width=10, height=6)
plot(ecdf(hd), 
     main="CDFs for T4 cell counts for Hodgkin's disease and non-Hodgkin's disease")
lines(ecdf(non_hd), col="blue")
dev.off()

dat <- rbind(summary(hd),
             summary(non_hd)) %>% data.table
dat$Group <- c("Hodgkin's", "Non-Hodgkin's")

sink("summary.tex")
print(xtable(dat, digits=c(0, 0, rep(1, 3), 0, 0, 0)),
      include.rownames=FALSE)
sink()

n <- length(hd)
x1 <- rep(1, 2*n)
x2 <- c(rep(0, n), rep(1, n))
X <- cbind(x1, x2)
y <- c(hd, non_hd)

mod_b1 <- lm(y ~ 1 + x2)
ci_b1 <- confint(mod_b1, level=0.9) %>% as.data.table %>% "["(2)

mod_b2 <- lm(log(y) ~ 1 + x2)
ci_b2 <- confint(mod_b2, level=0.9) %>% as.data.table %>% "["(2)

mod_b3 <- lm(sqrt(y) ~ 1 + x2)
ci_b3 <- confint(mod_b3, level=0.9) %>% as.data.table %>% "["(2)

ci_dat <- rbind(ci_b1, ci_b2, ci_b3)
ci_dat$Model <- c("Original", "Log(Y)", "Sqrt(Y)")

sink("ci.tex")
print(xtable(ci_dat, digits=c(0, 3, 3, 0)),
      include.rownames=FALSE)
sink()

pdf("resids.pdf", width=10, height=4)
par(mfrow=c(1, 3))
plot(mod_b1$resid, ylab="Residual", main="Original")
plot(mod_b2$resid, ylab="Residual", main="Log")
plot(mod_b3$resid, ylab="Residual", main="Sqrt")
dev.off()

pdf("qqplots.pdf", width=10, height=4)
par(mfrow=c(1, 3))
qqnorm(mod_b1$resid)
qqline(mod_b1$resid)
qqnorm(mod_b2$resid)
qqline(mod_b2$resid)
qqnorm(mod_b3$resid)
qqline(mod_b3$resid)
dev.off()

mod_c1 <- glm(y ~ 1 + x2, family=poisson(link="log"))
ci_c1 <- confint(mod_c1, level=0.9) %>% as.data.table %>% "["(2)

mod_c2 <- glm(y ~ 1 + x2, family=Gamma(link="inverse"))
ci_c2 <- confint(mod_c2, level=0.9) %>% as.data.table %>% "["(2)

mod_c3 <- glm(y ~ 1 + x2, family=inverse.gaussian(link="1/mu^2"))
ci_c3 <- confint(mod_c3, level=0.9) %>% as.data.table %>% "["(2)

ci_dat_2 <- rbind(ci_c1, ci_c2, ci_c3)
ci_dat_2$Model <- c("Poisson", "Gamma", "Inverse Gaussian")

sink("ci2.tex")
print(xtable(ci_dat_2, digits=c(0, 7, 7, 0)),
      include.rownames=FALSE)
sink()
