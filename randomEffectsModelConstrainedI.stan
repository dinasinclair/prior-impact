/*
saved as randomEffectsModelConstrainedI.stan
Dina Sinclair June 9 2018
This implements a randomized effects bayesian model under the constraint that
I (number of individuals per study) is a constant.
*/


data {
  int<lower=0> J; // number of individuals per study
  int<lower=0> I; // number of studies
  real Y[I, J]; // data points from study i individual j
  real<lower=0> sigmaSq[I, J]; // var of effect estimates 
}
parameters {
  real mu; 
  real<lower=0> tau;
  real theta[I];
}
transformed parameters {
}
model {
  for (i in 1:I){
    theta[i] ~ normal(mu, tau^2);
    for (j in 1:J){
      Y[i, j] ~ normal(theta[i], sigmaSq[i, j]);
    }
  }
}