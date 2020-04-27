// Log-ratio Gaussian model for microbiome associations
// Full model: inference of effect sizes and directions

data {
  int<lower=0> N; // number of samples
  int<lower=0> S; // number of species
  int<lower=0> p; // number of covariates
  matrix[N,p] X; // covariate data
  vector[S] z[N]; // log-ratio of species proportions
  matrix[S,S] Sigma; // covariance matrix of the species 
                            //(based on Brownian phylogenetics model)
  real<lower=0> noise; // noise parameter for effects (currently fixed)
  matrix[S,p] a_sgn; // Direction of species effects
}

transformed data {
  matrix[S,S] L;
  L = cholesky_decompose(Sigma);
}

parameters {
  vector<lower=0>[p] a_global; // Global effects
  matrix[S,p] a_raw;
}

transformed parameters {
  vector[S] mu[N];
  matrix[S,p] a; // Species effects

  a = a_sgn .* rep_matrix(a_global, S)' + noise*a_raw;
  
  for(i in 1:N) {
    mu[i] = a * X[i,]';
  }

}

model {
  a_global ~ normal(0,1);
  
  for(j in 1:p) {
    a_raw[,j] ~ normal(0, 1);
  }
  
  z ~ multi_normal_cholesky(mu, L);
  
}
