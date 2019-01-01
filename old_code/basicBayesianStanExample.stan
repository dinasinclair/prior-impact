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

/*
Unclear if should include or if generating this in R is fine
functions {
    /*
   		Return thetas drawn from a normal distribution determined by mu, tauSq
   */

   vector hierarchical_rng(real mu, real tauSq, int n) {
      matrix[n, n] Y; // the data to be generated
      vector [n] theta; // 
      
      // find theta
      for (i in 1:n)
        theta[n] <- normal(mu, tauSq);
      return theta;

      // Note that we don't have to generate the sigma, since those are assumed to be known.
   }
 }*/