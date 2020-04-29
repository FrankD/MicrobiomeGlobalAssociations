# Applying the full log-ratio model on some generated data

library(rstan)
library(here)
source(here('Simulation.R'))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()) # Parallel run on multiple cores

set.seed(2)

# Inputs
n = 100 # Sample size
S = 10 # Number of species
p = 2 # Number of covariates
a.global = c(1,0,0.5) # True values for covariates
noise = 0.1 # Noise 

# Generate some simulated data under log-ratio model
cov.matrix = diag(S-1)
sim.data = generate_data_phylo_log_ratio(n, p, S, cov.matrix, a.global, noise)

# Compile model (only need to do this once if sampling for multiple datasets)
st.model = stan_model(file=here('phylo_log_ratio_model.stan'))

# Get data 
x = sim.data$X
y = sim.data$counts

# Estimation of proportions and log ratio transformation
props = (sim.data$counts / rowSums(sim.data$counts)) + 0.0001 # added small constant to avoid zero counts
z = log(props[,1:(S-1)]/matrix(props[,S], dim(props)[1], dim(props)[2]-1))

mb_dat = list(S=S-1, N=n, p=p+1, noise=0.1, prior_mix_belief=c(0.5,0.5),
              X=x, z=z, Sigma=cov.matrix, tau=0.1)

# Sampling
fit <- sampling(st.model, data=mb_dat,
                iter=12000, chains=4, include=TRUE, 
                pars=c('a_global', 'a'))

# Save the output
saveRDS(paste0('result_full_model_n_',n,'_S_',S,'_p_',p,'_noise_',noise,'.RData'))
