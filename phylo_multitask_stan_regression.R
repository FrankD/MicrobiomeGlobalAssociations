# Multitask model with phylogenetic information using stan
# Idea is to have regression over outcomes, where outcomes are abundances of species in the microbiome.
# The species are related, so we want to use the degree of relatedness in the regression.
# This is done by putting a prior over the rows in the matrix of regression coefficients
# over all species, and using the phylogenetic tree branch lengths to specify a 
# covariance matrix for the prior.

library(mvtnorm)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

# Try three-species model
n = 30
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

mb_dat = list(S=3, N=n,
              p=1,
              X=matrix(x), y=y, Sigma=cov.matrix)

fit = stan(file='phylo_multitask_model.stan', data=mb_dat,
           iter=2000, chains=1)

model_func <- function(params,x,y,sigma,noise) {
  prob.a = dmvnorm(abs(params[1:3]), rep(params[4], 3), sigma, log=TRUE)
  prob.y = sum(sapply(1:dim(y)[1], function(i) 
    dmvnorm(t(y[i,]), x[i]*t(params[1:3]), noise*diag(3), log=TRUE)))
  
  
  return(-prob.a - prob.y)
}

test = optim(rnorm(4), model_func, gr=NULL, x, y, cov.matrix, noise)


# Now for a bigger example
library(ape)


S = 40

# Random tree
tree = rtree(n = S)

# Calculate Brownian covariance matrix
nb.tip = length(tree$tip.label)
dis = dist.nodes(tree)
MRCA = mrca(tree, full = FALSE)
M = dis[as.character(nb.tip + 1), MRCA]
dim(M) = rep(sqrt(length(M)), 2)

# branch lengths all 0.1
cov.matrix = M

# Note, model is now using covariance matrix
# for the responses, assuming univariate Gaussian 
# N(mu,sigma) for parameters - seems more sensible.
#
# Trying it out in stan, but suspect it needs a different
# sampler for indicator variables and/or sparsity prior on
# a; want to find parameter mu such that it quantifies the 
# overall effect in all species.


a.global = abs(rnorm(1))
a.sgn = rbinom(S,1,0.5)*2-1
a.sgn = rbinom(S,2,0.5)-1

x = rnorm(n)

y = matrix(NA, n,S)

a = rnorm(S, a.sgn*a.global, noise)

for(i in 1:n) {
  y[i,] = rmvnorm(1, x[i] * a, cov.matrix)
}

mb_dat = list(S=S, N=n,
              p=1,
              X=matrix(x), y=y, Sigma=cov.matrix)

fit = stan(file='phylo_multitask_model_mixture.stan', data=mb_dat,
           iter=4000, chains=4, include=FALSE, pars=c('mu_resp', 'a_sgn'))


