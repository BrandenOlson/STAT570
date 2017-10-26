data{
  int<lower=0> n;
  real<lower=0> y[n];
}
parameters{
  real<lower=0> eta;
  real<lower=0> alpha;
}
model{
  target += lognormal_lpdf(eta |  1.354025, 1.244592);
  target += lognormal_lpdf(alpha | 0.6931472, 0.4214036);
  target += weibull_lpdf(y | eta, alpha);
}
