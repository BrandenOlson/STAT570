data{
  int<lower=0> n;
  real<lower=0> y[n];
}
parameters{
  real<lower=0> eta;
  real<lower=0> alpha;
}
model{
  eta ~ lognormal(1.354025, 1.244592);
  alpha ~ lognormal(0.6931472, 0.4214036);
  y ~ weibull(eta, alpha);
}
