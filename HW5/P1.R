library(data.table)
library(dplyr)
library(xtable)

X0Y1 <- 104
X0Y0 <- 666
X1Y1 <- 96
X1Y0 <- 109

n1 <- X0Y1 + X1Y1
n2 <- X0Y0 + X1Y0

# Part b
p1 <- X1Y1/n1
p2 <- X1Y0/n2

theta <- p1/(1 - p1)*(1 - p2)/p2
se_theta <- sqrt(p1*(1 - p2)*(p2*(1 - p2) + p1*(1 - p1))/(p2^3*(1 - p1)^3*(n1 + n2)))
ci_left <- theta - qnorm(0.95)*se_theta
ci_right <- theta + qnorm(0.95)*se_theta

# Case
a <- 1
b <- 1
post_a <- a + X1Y1
post_b <- b + n1 - X1Y1
cred_left_case <- qbeta(0.05, post_a, post_b)
cred_right_case <- qbeta(0.95, post_a, post_b)

mean_case <- (X1Y1 + 1)/(n1 + 2)
mode_case <- X1Y1/n1
sd_case <- sqrt((X1Y1 + 1)*(n1 - X1Y1 + 1))/((n1 + 2)*sqrt(n1 + 3))

post_a_control <- a + X1Y0
post_b_control <- b + n2 - X1Y0
cred_left_control <- qbeta(0.05, post_a_control, post_b_control)
cred_right_control <- qbeta(0.95, post_a_control, post_b_control)

mean_control <- (X1Y0 + 1)/(n2 + 2)
mode_control <- X1Y0/n2
sd_control <- sqrt((X1Y0 + 1)*(n2 - X1Y0 + 1))/((n2 + 2)*sqrt(n2 + 3))
# Control

dat <- data.table(Mean=c(mean_case, mean_control),
                  Mode=c(mode_case, mode_control),
                  SD=c(sd_case, sd_control),
                  "0.05 Quantile"=c(cred_left_case, cred_left_control),
                  "0.95 Quantile"=c(cred_right_case, cred_right_control))
dat$Group <- c("Case", "Control")

sink("bayes.tex")
print(xtable(dat, digits=c(0, rep(3, 5), 0)),
      include.rownames=FALSE)
sink()

asymp_left_case <- qnorm(0.05, p1, p1*(1 - p1)/sqrt(n1))
asymp_right_case <- qnorm(0.95, p1, p1*(1 - p1)/sqrt(n1))

asymp_left_control <- qnorm(0.05, p2, p2*(1 - p2)/sqrt(n2))
asymp_right_control <- qnorm(0.95, p2, p2*(1 - p2)/sqrt(n2))

dat_asymp <- data.table("0.05 Quantile"=c(asymp_left_case,
                                          asymp_left_control),
                        "0.95 Quantile"=c(asymp_right_case,
                                          asymp_right_control))
dat_asymp$Group <- c("Case", "Control")

sink("asymp.tex")
print(xtable(dat_asymp, digits=c(0, rep(3, 2), 0)),
      include.rownames=FALSE)
sink()

sample_size <- 1000
case_samp <- rbeta(n=sample_size, post_a, post_b)
case_samp_qs <- quantile(case_samp, probs=c(0.05, 0.95))

pdf("hist_case.pdf", width=10, height=6)
hist(case_samp, xlab="p1", main="Histogram of case posterior", breaks=20)
dev.off()

control_samp <- rbeta(n=sample_size, post_a_control, post_b_control)
control_samp_qs <- quantile(control_samp, probs=c(0.05, 0.95))

pdf("hist_control.pdf", width=10, height=6)
hist(control_samp, xlab="p2", main="Histogram of control posterior", breaks=20)
dev.off()

sample_dat <- data.table(rbind(case_samp_qs, control_samp_qs))
sample_dat$Group <- c("Case", "Control")

sink("samp.tex")
print(xtable(sample_dat, digits=c(0, rep(3, 2), 0)),
      include.rownames=FALSE)
sink()

theta_samp <- (case_samp/(1 - case_samp))/(control_samp/(1 - control_samp))
theta_qs <- quantile(theta_samp, probs=c(0.05, 0.5, 0.95))

pdf("theta.pdf", width=10, height=6)
hist(theta_samp, xlab="theta", main="Histogram of odds ratio", breaks=20)
dev.off()

sink("theta.tex")
print(xtable(data.table(rbind(theta_qs)), digits=c(0, rep(3, 3))),
      include.rownames=FALSE)
sink()
