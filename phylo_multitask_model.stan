// First version of the model
// only incorporates linear regression

data {
  int<lower=0> N; // number of samples
  int<lower=0> S; // number of species
  int<lower=0> p; // number of covariates
  matrix[N,p] X; // covariate data
  vector[S] y[N]; // observed species abundances
  matrix[S,S] Sigma; // covariance matrix of the species 
                            //(based on Brownian phylogenetics model)
}
parameters {
  vector[p] mu; 
  matrix[S,p] a;
  real a_sgn[S,p];
  real noise;
}

transformed parameters {
  matrix[S,p] a_sgn_trans;
  row_vector[S] mu_resp[N];
  
  for(n in 1:N) {
    mu_resp[n] = X[n,] * a[,]';
  }
  
  for(s in 1:S) {
    for(j in 1:p) {
      a_sgn_trans[s,j] = step(a_sgn[s,j])*2-1;
    }
  }
}

model {
  
  mu ~ gamma(1,1);
  
  for(j in 1:p) {
    a_sgn[,j] ~ normal(0,1);
    a[,j] ~ normal(a_sgn_trans[,j] .* rep_vector(mu[j],S), noise);
  }
  
  
  y ~ multi_normal(mu_resp, Sigma);
  //for(s in 1:S) {
  //  y[,s] ~ multi_normal(X * a[s,]', noise);
  //}
}
