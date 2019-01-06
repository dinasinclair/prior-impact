// A two layer hierarchical model with I studies of fixed size
// J. Used https://discourse.mc-stan.org/t/simple-question-on-hierarchical-non-centered-parameterization/2216 to help think through the correct coding strategy.

data {
  int<lower=0> I; // number of studies
  int<lower=0> N; // Number of groups
  real Y_mean[I]; // mean from study i 
  real<lower=0> Y_sd[I]; // sd from study i
  int<lower=0> groups[I]; // Group assignment of each study i

}
parameters {
  real mu; // population mean
  real<lower=0> tau; // population s.d.
  real G_mean[N]; // mean from group g
  real<lower=0> G_sd[N]; // sd from group g
  real theta[I];
}
transformed parameters {
}
model {
  // First calculate group mean/sd
  for (n in 1:N){
    target += normal_lpdf(G_mean[n]| mu, tau);
    target += lognormal_lpdf(G_sd[n] | 0, 1);
  }

  // Then calculate study mean based on group results (we assume Y_sd is a known constant)
  for (i in 1:I){
      target += normal_lpdf(theta[i] | G_mean[groups[i]], G_sd[groups[i]]);
      target += normal_lpdf( Y_mean[i] | theta[i],  Y_sd[I]);
  }
}
