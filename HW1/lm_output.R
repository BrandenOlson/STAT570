library(lasso2)
data(Prostate)
attach(Prostate)
y <- Prostate$lpsa
lcavol <- Prostate$lcavol
svi <- Prostate$svi
lmod <- lm(y~lcavol+svi+svi:lcavol)
summary(lmod)
