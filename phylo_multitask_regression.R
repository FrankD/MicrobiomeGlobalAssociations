# Simple test for the multitask model with phylogenetic information
# Idea is to have regression over outcomes, where outcomes are abundances of species in the microbiome.
# The species are related, so we want to use the degree of relatedness in the regression.
# This is done by putting a prior over the rows in the matrix of regression coefficients
# over all species, and using the phylogenetic tree branch lengths to specify a 
# covariance matrix for the prior.

library('mvtnorm')

# Try three-species model
n = 100
noise = 0.01

x = rnorm(n)

# branch lengths all 0.1
cov.matrix = matrix(0, 3, 3)
cov.matrix[1,1] = 0.2; cov.matrix[2,2] = 0.2
cov.matrix[3,3] = 0.1
cov.matrix[1,2] = 0.1; cov.matrix[2,1] = 0.1

a.global = abs(rnorm(1))
a.sgn = rbinom(3,1,0.5)*2-1

a = rmvnorm(1, a.sgn*a.global, cov.matrix)

y = x %*% a + matrix(rnorm(3*n,0,noise), n, 3)


model_func <- function(params,x,y,sigma,noise) {
  prob.a = dmvnorm(abs(params[1:3]), rep(params[4], 3), sigma, log=TRUE)
  prob.y = sum(sapply(1:dim(y)[1], function(i) 
    dmvnorm(t(y[i,]), x[i]*t(params[1:3]), noise*diag(3), log=TRUE)))
  
  
  return(-prob.a - prob.y)
}

test = optim(rnorm(4), model_func, gr=NULL, x, y, cov.matrix, noise)
