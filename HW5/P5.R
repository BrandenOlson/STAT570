library(dplyr)
library(INLA)
library(xtable)

# Data pre-processing - grabbed from book website
lung <- read.table("MNlung.txt", header=TRUE, sep="\t")
radon <- read.table("MNradon.txt", header=TRUE)
Obs <- apply(cbind(lung[,3], lung[,5]), 1, sum)
Exp <- apply(cbind(lung[,4], lung[,6]), 1, sum)
rad.avg <- rep(0, length(lung$X))
for(i in 1:length(lung$X)) {
            rad.avg[i] <- mean(radon[radon$county==i,2])
}
x <- rad.avg
rad.avg[26]<-0
rad.avg[63]<-0
x[26] <- NA
x[63] <- NA
newy <- Obs[is.na(x)==F]
newx <- x[is.na(x)==F]
newE <- Exp[is.na(x)==F]

inla_df <- data.frame(cancer=newy, radon=newx, Expected=newE)

lin.mod <- inla(cancer ~ 1 + radon, data=inla_df, family="poisson", E=Expected)

plot(lin.mod)

sum_dat <- lin.mod$summary.fixed[, -c(6, 7)]

sink("inla.tex")
print(xtable(sum_dat, digits=c(0, rep(3, 5))),
      include.rownames=FALSE)
sink()
