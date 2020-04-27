// Only linear regression
// Mixture model to capture effect a, -a or 0
// Marginalization over components as required by stan

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
  simplex[3] theta[p]; // mixing proportions
  matrix[S,p] a;
  real a_global[p];
  real<lower=0> noise[2];
}

transformed parameters {
  row_vector[S] mu_resp[N];
  
  // Mean for each 
  for(n in 1:N) {
    mu_resp[n] = X[n,] * a[,]';
  }
}

model {
  
  a_global ~ gamma(1,1);
  noise ~ exponential(10);
  
  for(j in 1:p) {
    vector[3] log_theta = log(theta[j]);
    
    for(s in 1:S) {
      vector[3] lps = log_theta; 
      
      // Three mixture components (negative, zero, positive)
      lps[1] = lps[1] + normal_lpdf(a[s,j] | -a_global[j], noise[1]);
      lps[2] = lps[2] + normal_lpdf(a[s,j] | 0, noise[2]);
      lps[3] = lps[3] + normal_lpdf(a[s,j] | a_global[j], noise[1]);
      
      target += log_sum_exp(lps);
    }
  }
  
  
  y ~ multi_normal(mu_resp, Sigma);
  //for(s in 1:S) {
  //  y[,s] ~ multi_normal(X * a[s,]', noise);
  //}
}
