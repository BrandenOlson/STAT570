library(data.table)
library(dplyr)
library(ggplot2)
library(lmtest)
library(MASS)
library(rootSolve)
library(xtable)

dat <- read.csv("mental.csv", header=TRUE)
y <- dat$MI - 1 # THe data come in values 1-4, but the problem assumes values 0 - 3
x1 <- dat$LE
x2 <- dat$SES

# Problem 1
z <- ifelse(y %in% c(0, 1), 0, 1)

matrix_1 <- matrix(NA, 0, 3)
    
x1_vals <- 0:9
for(x in x1_vals) {
    x_prob_0 <- mean(z[x1 == x & x2 == 0])
    x_prob_1 <- mean(z[x1 == x & x2 == 1])
    matrix_1 <- rbind(matrix_1, c(x_prob_0, x, 0))
    matrix_1 <- rbind(matrix_1, c(x_prob_1, x, 1))
}

dat_1 <- matrix_1 %>%
    data.table %>%
    setNames(c("Frequency", "LE", "SES"))
dat_1$SES <- as.factor(dat_1$SES)
p <- ggplot(dat_1, aes(x=LE, y=Frequency, colour=SES)) + 
    geom_line() +
    ggtitle("Empirical frequency of Z by life event severity and socioeconomic status") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("plot1.pdf", width=10, height=6)

# Problem 3
mod_3 <- glm(z ~ 1 + x1 + x2, family=binomial(link="logit"))
ci_3 <- confint(mod_3)

dat_3 <- cbind(mod_3$coef, ci_3) %>%
    data.table(keep.rownames=TRUE) %>%
    setNames(c("Covariate", "Estimate", "2.5\\%", "97.5\\%"))

sink("mod3.tex")
print(xtable(dat_3), digits=c(rep(0, 2), rep(3, 3)),
      include.rownames=FALSE,
      sanitize.text.function=function(x) { x })
sink()

odds_ratio_ci <- ci_3 %>% exp
odds_dat <- odds_ratio_ci %>%
    data.table(keep.rownames=TRUE) %>%
    setNames(c("Covariate", "2.5\\%", "97.5\\%"))

sink("odds.tex")
print(xtable(odds_dat), digits=c(rep(0, 2), rep(3, 2)),
      include.rownames=FALSE,
      sanitize.text.function=function(x) { x })
sink()

# Problem 6
matrix_6 <- matrix(NA, 0, 4)
x1_vals <- 0:9
x2_vals <- c(0, 1)
y_vals <- 0:3
for(x1_val in x1_vals) {
    for(x2_val in x2_vals) {
        y_sub <- y[x1 == x1_val & x2 == x2_val]
        y_sub_len <- length(y_sub)
        for(y_val in y_vals) {
            y_val_prob <- length(y_sub[y_sub == y_val])/y_sub_len
            matrix_6 <- rbind(matrix_6, c(y_val_prob, y_val, x1_val, x2_val))
        }
    }
}

dat_6 <- matrix_6 %>%
    data.table %>%
    setNames(c("Frequency", "MI", "LE", "SES"))
dat_6$MI <- as.factor(dat_6$MI)
dat_6$SES <- as.factor(dat_6$SES)
p2 <- ggplot(dat_6, aes(x=LE, y=Frequency, colour=MI, linetype=SES)) + 
    geom_line() +
    ggtitle("Empirical frequency of MI by life event severity and socioeconomic status") +
    theme(plot.title = element_text(hjust = 0.5))
ggsave("plot2.pdf", width=10, height=6)

# Problem 7
y <- as.factor(y)

mod_7a <- polr(y ~ 1)
mod_7b <- polr(y ~ 1 + x1)
mod_7c <- polr(y ~ 1 + x2)
mod_7d <- polr(y ~ 1 + x1 + x2)

getPValue <- function(lr_test) {
    return(lr_test$Pr[2])
}

test_ab <- lrtest(mod_7a, mod_7b)
test_ac <- lrtest(mod_7a, mod_7c)
test_ad <- lrtest(mod_7a, mod_7d)
test_bd <- lrtest(mod_7b, mod_7d)
test_cd <- lrtest(mod_7c, mod_7d)

lr_mat <- matrix(NA, 0, 3)
lr_mat <- rbind(lr_mat, 
                c(1, 2, test_ab %>% getPValue),
                c(1, 3, test_ac %>% getPValue),
                c(1, 4, test_ad %>% getPValue),
                c(2, 4, test_bd %>% getPValue),
                c(3, 4, test_cd %>% getPValue)
                )

lr_dat <- lr_mat %>%
    data.table %>%
    setNames(c("Reduced", "Full", "P-value"))

sink("lr.tex")
print(xtable(lr_dat, digits=c(0, 0, 0, 4)),
      include.rownames=FALSE)
sink()

betas <- mod_7d$coef
alphas <- mod_7d$zeta
ses <- summary(mod_7d)[[1]][, 2] %>% unname

coef_dat <- data.table(Parameter=c(
                                   "$\\estim\\beta_1^{(12)}$",
                                   "$\\estim\\beta_2^{(12)}$",
                                   "$\\estim\\alpha_0^{(12)}$",
                                   "$\\estim\\alpha_1^{(12)}$",
                                   "$\\estim\\alpha_2^{(12)}$"
                                   ),
                       Estimate=c(betas, alphas),
                       SE=ses)

sink("coef.tex")
print(xtable(coef_dat, digits=c(0, 0, 3, 3)),
      include.rownames=FALSE,
      sanitize.text.function=function(x) { x })
sink()


# Problem 8
plot_positions <- c("topright", "bottomleft", "bottomright", "topleft")
pdf("probs.pdf", width=10, height=8)
par(mfrow=c(2, 2))
for(i in 1:4) {
    plot(mod_7d$fitted.values[, i] ~ x1, col=as.factor(x2), 
         pch=ifelse(x2, 19, 15), xlab="LE", ylab="Estimated probability",
         main=paste("Probability that Y =", i, "by LE and SES"))
    legend(plot_positions[i], c("Low SES", "High SES"),
           col=as.factor(c(0, 1)),
           pch=c(15, 19),
           )
}
dev.off()
