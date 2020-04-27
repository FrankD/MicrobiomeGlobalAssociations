# Apply the full log_ratio model to simulated datasets with small S and small p

# TODO: 
# - Generate a number of datasets with small S and small p, and estimate a.global 
# - Make plot showing how estimation error of a.global changes with n and noise

library(rstan)
library(here)
source(here('Code/simulation.R'))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

set.seed(10)

S = 6
p = 3
n = 400
a.global = c(1,2,0,0.5)
noise = 0.5

# Generate sim data under log ratio model
#cov.matrix = generate_cov_matrix(S-1)
cov.matrix = diag(S-1)
sim.data = generate_data_phylo_log_ratio(n, p, S, cov.matrix, a.global, noise)

# Compile model (only need to do this once if sampling for multiple datasets)
st.model = stan_model(file=here('Code/phylo_log_ratio_model.stan'))
  
# Get data 
x = sim.data$X
y = sim.data$counts

# Estimation of proportions and log ratio transformation
props = (sim.data$counts / rowSums(sim.data$counts)) + 0.0001 # added small constant to avoid zero counts
z = log(props[,1:5]/matrix(props[,6], dim(props)[1], dim(props)[2]-1))

mb_dat = list(S=S-1, N=n, p=p+1, noise=0.1, prior_mix_belief=c(0.5,0.5),
              X=x, z=z, Sigma=cov.matrix, tau=0.1, 
              a_global=a.global)

# Sampling
fit = sampling(st.model, data=mb_dat,
             iter=2000, chains=4, include=TRUE, 
             pars=c('a_global', 'a'))
  
fit = sampling(st.model, data=mb_dat, init=get_inits(fit2),
               seed=get_seed(fit2),
               iter=2000, chains=4, include=TRUE, 
               pars=c('a_global', 'a'))

