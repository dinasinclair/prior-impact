// A one layer hierarchical model with I studies of constant
// but alterable size J.

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
    target += normal_lpdf(theta[i] | mu, tau^2);
    for (j in 1:J){
      target += normal_lpdf( Y[i, j] | theta[i], sigmaSq[i, j]);
    }
  }
}
