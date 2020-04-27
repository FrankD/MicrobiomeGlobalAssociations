// Only linear regression
// Mixture model to capture effect a, -a or 0
// Marginalization over components as required by stan
// Uses the REBAR approach for approximating
// the categorical latent variable

functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
      real alpha_plus = sum(alpha);
      
      return(lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y))) - 
               lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha)));
  }
}

data {
  int<lower=0> N; // number of samples
  int<lower=0> S; // number of species
  int y[N,S]; // observed species abundances
}



parameters {
  vector[S] gamma;
}

model {
  gamma ~ gamma(1.0, 1.0);
  
  for(n in 1:N) {
    y[n,] ~ dirichlet_multinomial(gamma);
  }
}
