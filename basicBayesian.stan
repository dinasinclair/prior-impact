// saved as basicBayesian.stan
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
