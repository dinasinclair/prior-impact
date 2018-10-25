// A two layer hierarchical model with I studies of fixed size
// J. TODO this is currently starter code from the 1D model.

data {
  int<lower=0> I; // number of studies
  real Y[I]; // data points from study i individual j
  real<lower=0> sigmaSq[I]; // var of effect estimates 
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
    target += normal_lpdf( Y[i] | theta[i], sigmaSq[i]);
  }
}
