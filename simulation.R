# Functions for generating simulated phylogenetic trees and microbiome data, 
# using the ape and MultiRNG packages.

library(ape)
library(mvtnorm)
library(MultiRNG)
library(dirmult)

# Generate a covariance matrix based on a random phylogenetic tree, for a 
# given number of species S, under the assumption of a Brownian noise model.
generate_cov_matrix <- function(S) {
  # Random tree
  tree = rtree(n = S)
  
  # Calculate Brownian covariance matrix
  nb.tip = length(tree$tip.label)
  dis = dist.nodes(tree)
  MRCA = mrca(tree, full = FALSE)
  M = dis[as.character(nb.tip + 1), MRCA]
  dim(M) = rep(sqrt(length(M)), 2)
  
  # branch lengths all 0.1
  return(M)
}

# Generate parameters (effect sizes), given a 
# covariance matrix and global effects for the covariates.
generate_hier_effect_sizes <- function(n, p, S, cov.matrix, 
                                       a.global, noise) {
  
  x = matrix(rnorm(p*n), n, p)
  
  a = matrix(NA, p, S)
  
  a.sgn = matrix(NA, p, S)
  
  for(j in 1:p) {
    a.sgn[j,] = rbinom(S,2,0.5)-1
    a[j,] = rmvnorm(1, a.sgn[j,]*a.global[j], noise*cov.matrix)
  }
  
  return(list(X=x, a.sgn=a.sgn, a=a))
  
}

# Generate random data from the log ratio model, given a 
# covariance matrix and global effects for the covariates.
# Note that phylogenetics is determining the covariance of the 
# log ratios, not the effect parameters.
generate_data_phylo_log_ratio <- function(n, p, S, cov.matrix, 
                                          a.global, noise, 
                                          total.num=100*S) {
  x = cbind(1, matrix(rnorm(p*n), n, p))
  
  y = matrix(NA, n, S-1)
  
  a = array(NA, c(n, p+1, S-1))
  
  a.sgn = matrix(NA, S-1, p+1)
  
  I = diag(S-1)
  
  for(j in 1:(p+1)) {
    a.sgn[,j] = rbinom(S-1,2,0.5)-1
    a[,j,] = rmvnorm(n, a.sgn[,j]*a.global[j], noise*I)
  }
  
  for(i in 1:n) {
    mu = x[i,] %*% a[i,,]
    y[i,] = rmvnorm(1, mu, cov.matrix)
  }
  
  y_exp = exp(y)
  y_exp_sum = rowSums(y_exp) + 1
  probs = y_exp/y_exp_sum
  probs = cbind(probs, 1/y_exp_sum)

  counts = matrix(NA, n, S)
  for(i in 1:n) {
    counts[i,] = rmultinom(1, total.num, probs[i,])
  }
  
  return(list(X=x, y=y, probs=probs, counts=counts, a.sgn=a.sgn, a=a))
  
}



# Generate random data from the multi-nomial dirichlet model, given a 
# covariance matrix and global effects for the covariates.
generate_data_phylo_dirichlet <- function(n, p, S, cov.matrix, 
                                          a.global, noise, alpha, 
                                          intercept.effect=1,
                                          total.num=100*S) {
  #x = cbind(1, matrix(rnorm(p*n), n, p))
  x = matrix(rnorm(p*n), n, p)
  
  y = matrix(NA, n,S)
  
  a = array(NA, c(n, p, S))
  
  a.sgn = matrix(NA, S, p)
  
  # intercept
  #a[,1,] = rmvnorm(n, rep(intercept.effect, S), noise*cov.matrix)
  
  for(j in 1:p) {
    a.sgn[,j] = rbinom(S,2,0.5)-1
    a[,j,] = rmvnorm(n, a.sgn[,j]*a.global[j], noise*cov.matrix)
  }
  
  for(i in 1:n) {
    gamma = alpha*exp(x[i,] %*% a[i,,])/sum(exp(x[i,] %*% a[i,,]))
    #alpha = exp(x[i,] %*% a[i,,])
    #y[i,] = draw.dirichlet.multinomial(1, S, alpha, 1, 1000)
    probs.temp = rdirichlet(1, gamma)
    y[i,] = rmultinom(1, total.num, probs.temp)
  }  
  
  return(list(X=x, y=y, a.sgn=a.sgn, a=a))
  
}

# Generate random data from the negative binomial model, given a 
# covariance matrix and global effects for the covariates.
generate_data_phylo_neg_binom <- function(n, p, S, cov.matrix, 
                                          a.global, noise, alpha, 
                                          intercept.effect=1,
                                          total.num=100*S) {
  #x = cbind(1, matrix(rnorm(p*n), n, p))
  x = matrix(rnorm(p*n), n, p)
  
  y = matrix(NA, n,S)
  
  a = array(NA, c(n, p, S))
  
  a.sgn = matrix(NA, S, p)
  
  # intercept
  #a[,1,] = rmvnorm(n, rep(intercept.effect, S), noise*cov.matrix)
  
  for(j in 1:p) {
    a.sgn[,j] = rbinom(S,2,0.5)-1
    a[,j,] = rmvnorm(n, a.sgn[,j]*a.global[j], noise*cov.matrix)
  }
  
  for(i in 1:n) {
    mu = exp(x[i,] %*% a[i,,])
    y[i,] = rnbinom(S, 1, mu=mu)
  }  
  
  return(list(X=x, y=y, a.sgn=a.sgn, a=a))
  
}

# Generate random data from the multinomial dirichlet model,
# with fixed gamma parameters. Uses either MultiRNG package
# or combination of rdirichlet and rmultinom.
generate_data_dirichlet_mn <- function(n, S, gamma, use.package=FALSE) {
  y = matrix(NA, n,S)
  
  for(i in 1:n) {
    if(use.package) {
      y[i,] = draw.dirichlet.multinomial(1, S, gamma, 1, 500)
    } else {
      probs.temp = rdirichlet(1, gamma)
      y[i,] = rmultinom(1, 100*S, probs.temp)
    }
  }  
  
  return(y)
  
}

# Generate random data from the multinomial dirichlet regression model,
# with gamma determined by regression on set of covaraites. Uses either 
# MultiRNG package or combination of rdirichlet and rmultinom.
generate_data_dirichlet_mn_regression <- function(n, S, x, a, b, use.package=FALSE) {
  y = matrix(NA, n,S)
  
  xa = exp(x %*% a)
  
  for(i in 1:n) {
    gamma = b*xa[i,]/sum(xa[i,])
    #gamma = xa[i,]
    
    if(use.package) {
      y[i,] = draw.dirichlet.multinomial(1, S, gamma, 1, 500)
    } else {
      probs.temp = rdirichlet(1, gamma)
      y[i,] = rmultinom(1, 100*S, probs.temp)
    }
  }  
  
  return(y)
  
}


# Generate random data from the multinomial dirichlet model,
# with fixed alpha parameters.
generate_data_dirichlet <- function(n, S, alpha) {
  y = matrix(NA, n,S)
  
  for(i in 1:n) {
    y[i,] = draw.dirichlet.multinomial(1, S, alpha, 1, 500)
  }  
  
  return(y)
  
}


# Generate random data from the multi-nomial dirichlet model,
# with alpha defined by a linear combination of covariates
generate_data_dirichlet_regression <- function(n, S, X, a) {
  y = matrix(NA, n,S)
  
  for(i in 1:n) {
    alpha = exp(X[i,] %*% a)
    y[i,] = draw.dirichlet.multinomial(1, S, alpha, 1, 500)
  }  
  
  return(y)
  
}

# Generate random data from the multi-nomial dirichlet model,
# with alpha defined by a linear combination of covariates
# and the regression parameters drawn from independent Gaussian hyperpriors
generate_data_dirichlet_hierarchical_mv_regression <- 
  function(n, S, p, X, a.global, cov.matrix, noise=0.1) {
    
  y = matrix(NA, n,S)
  a = array(NA, c(n, p, S))
  
  for(j in 1:p) {
    a[,j,] = rmvnorm(n, rep(a.global[j],S), noise*cov.matrix)
  }
  
  for(i in 1:n) {
    alpha = exp(X[i,] %*% a[i,,])
    y[i,] = draw.dirichlet.multinomial(1, S, alpha, 1, 500)
  }  
  
  return(y)
  
}

# Generate random data from the multi-nomial dirichlet model,
# with alpha defined by a linear combination of covariates
# and the regression parameters drawn from a multivariate Gaussian hyperprior
generate_data_dirichlet_hierarchical_regression <- function(n, S, p, X, a.global,
                                                            noise=0.1) {
  y = matrix(NA, n,S)
  a = array(NA, c(n, p, S))
  
  for(j in 1:p) {
    a[,j,] = rnorm(n*S, a.global[j], noise)
  }
  
  for(i in 1:n) {
    alpha = exp(X[i,] %*% a[i,,])
    y[i,] = draw.dirichlet.multinomial(1, S, alpha, 1, 500)
  }  
  
  return(y)
  
}
