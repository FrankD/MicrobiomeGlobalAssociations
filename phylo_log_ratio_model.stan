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
  real tau; // concentration parameter
  vector[2] prior_mix_belief; // Prior belief about the mixture frequency 
                              // of -1, 1 effect direction. 
}

transformed data {
  matrix[S,S] L;
  L = cholesky_decompose(Sigma);
}

parameters {
  vector<lower=0>[p] a_global; // Global effects
  matrix[S,p] a_raw;
  simplex[2] theta[p]; // mixing proportions
  real G[2,S,p]; // Gumbel parameters
}

transformed parameters {
  vector[S] mu[N];
  matrix[S,p] a; // Species effects
  vector[2] one_hot[S,p]; // Direction of species effects
  matrix[S,p] mean_eff; // Mean of species effects

  for(j in 1:p) {
    for (s in 1:S){
      one_hot[s,j] = softmax((log(theta[j]) + to_vector(G[,s,j]))/tau);
      mean_eff[s,j] = one_hot[s,j][1]*a_global[j] - one_hot[s,j][2]*a_global[j];
    }
     
    a[,j] = mean_eff[,j] + noise*a_raw[,j];
  }
  
  for(i in 1:N) {
    mu[i] = a * X[i,]';
  }

}

model {
  a_global ~ normal(0,1);
  
  for(j in 1:p) {
    theta[j] ~ dirichlet(prior_mix_belief);
    a_raw[,j] ~ normal(0, 1);
    
    for(cat_i in 1:2) {
      G[cat_i,,j] ~ gumbel(0,1);
    }
  }
  
  z ~ multi_normal_cholesky(mu, L);
  
}
