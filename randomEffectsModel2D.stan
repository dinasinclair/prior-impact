// A two layer hierarchical model with I studies of fixed size
// J. Used https://discourse.mc-stan.org/t/simple-question-on-hierarchical-non-centered-parameterization/2216 to help think through the correct coding strategy.

data {
  int<lower=0> I; // number of studies
  real Y[I]; // data points from study i individual j
  real<lower=0> sigma[I]; // var of effect estimates 
  int<lower=0> g[I]; // Group assignment of each study
  int<lower=0> G; // Number of groups
}
parameters {
  real mu; // population mean
  real<lower=0> tau; // population s.d.
  real eta[G]; // group level errors
}
transformed parameters {
real theta[G]; // Group effects
for (x in 1:G)
  theta[x] = mu + tau * eta[x];
}
model {
  target += normal_lpdf(eta | 0,1);
  for (i in 1:I){
    y[i] ~ normal(theta[g[i]], sigma[i])
  }
}
