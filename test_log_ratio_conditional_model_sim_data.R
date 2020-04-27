# Apply the conditional log_ratio model to simulated datasets with small S and small p
# Conditional here means that we specify the direction of the effects (this can be
# estimated in a pre-processing step) and only infer the size and the global effect.

# TODO: 
# - Add in preprocessing step to determine effect sign (probably just fit linear models independently).
# - Generate a number of datasets with varying S, and estimate a.global 
# - Make plot showing how estimation error of a.global changes with S, n and noise
#
# NB: This still takes a long time for S>100; it may be necessary to set some parameters to zero
# during pre-processing.

library(rstan)
library(here)
source(here('Code/simulation.R'))
source(here('Code/linear_model.R'))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

set.seed(10)

S = 6
p = 3
n = 500
a.global = c(1,2,0,0.5)
noise = 0.8
exact = FALSE # Use true sign of a parameters, otherwise fit linear models.

# Generate sim data under log ratio model
cov.matrix = generate_cov_matrix(S-1)
cov.matrix = diag(S-1)
sim.data = generate_data_phylo_log_ratio(n, p, S, cov.matrix, a.global, noise)

# Compile model (only need to do this once if sampling for multiple datasets)
st.model = stan_model(file=here('Code/phylo_conditional_log_ratio_model.stan'))
  
# Get data 
x = sim.data$X
y = sim.data$counts

# Estimation of proportions and log ratio transformation
props = (sim.data$counts / rowSums(sim.data$counts)) + 0.0001 # added small constant to avoid zero counts
z = log(props[,1:(S-1)]/matrix(props[,S], dim(props)[1], dim(props)[2]-1))

if(exact) {
  # Use true sign (except where zero)
  a_sgn = sim.data$a.sgn
  # Set sign of zero parameters at random
  a_sgn[a_sgn==0] = (runif(sum(a_sgn==0))<0.5)*2-1 
} else {
  # Use linear model to determine signs marginally for each species
  a_sgn = fit_linear_models(z, x)
}

mb_dat = list(S=S-1, N=n, p=p+1, noise=0.1, prior_mix_belief=c(0.5,0.5),
              X=x, z=z, Sigma=cov.matrix, tau=0.1, a_sgn=sim.data$a.sgn,
              a_global=a.global)

# Sampling
fit = sampling(st.model, data=mb_dat,
             iter=2000, chains=4, include=TRUE, 
             pars=c('a_global', 'a'))
  


