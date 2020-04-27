# Apply the full model to simulated datasets with small S and small p

library(rstan)
library(here)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

set.seed(10)

on.hec = FALSE

alpha.vals = c(1, 5, 10)
n.vals = c(50, 100, 150, 200, 250)
dataset.vals = 1:10

settings = expand.grid(alpha.vals, n.vals, dataset.vals)

alpha = settings[i,1]
n = settings[i,2]
dataset.i = settings[i,3]

dataset = readRDS(here(paste0('Data/Simulation/small_s_alpha_', alpha, '.rds')))

data.indices = sample(500, n)

st.model = stan_model(file=here('Code/phylo_dirichlet_multitask_model_full.stan'))

cat('Dataset:', dataset.i, '\n',
    'n=', n, '\n',
    'alpha=', alpha, '\n')

sim.data = dataset[[dataset.i]]
  
x = cbind(1, sim.data$x[data.indices,])
y = sim.data$y[data.indices,]
cov.matrix = sim.data$cov.matrix
  
S = dim(y)[2]
p = dim(x)[2]

k = 20
a = 4

mb_dat = list(S=S, N=n, p=p, noise=0.01, prior_mix_belief=c(a,k-2*a,a),
              X=x, y=y, Sigma=cov.matrix, tau=0.5)

fit = sampling(st.model, data=mb_dat,
             iter=8000, chains=4, include=TRUE, 
             pars=c('a_global', 'alpha', 'theta', 'a', 'one_hot'))
  
if(!on.hec) {
  saveRDS(fit, file=here(paste0('Results/Simulation/Full/small_s_alpha_', alpha, '_n_', n, '_dataset_', dataset.i, '.rds')))
} else {
  saveRDS(fit, file=paste0('/storage/users/dondelin/Work/Projects/IBDProject/Results/Simulation/Full/small_s_alpha_', alpha, '_n_', n, '_dataset_', dataset.i, '.rds'))
}
  



