# Apply the full model to simulated datasets with small S and small p

library(rstan)
library(here)

source(here('Code/apply_DESeq2.R'))

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

set.seed(10)

on.hec = FALSE

alpha.vals = c(1, 5, 10)
n.vals = c(50, 100, 150, 200, 250)
S.vals = c(20, 50, 100, 250, 500, 1000, 2500)
dataset.vals = 1:10

settings = expand.grid(n.vals, S.vals, dataset.vals, alpha.vals)

alpha = settings[i,4]
n = settings[i,1]
S = settings[i,2]
dataset.i = settings[i,3]

dataset = readRDS(here(paste0('Data/Simulation/large_s_', S, '_alpha_', alpha, '.rds')))

data.indices = sample(500, n)

st.model = stan_model(file=here('Code/phylo_dirichlet_multitask_model_fixed.stan'))

cat('Dataset:', dataset.i, '\n',
    'n=', n, '\n',
    'alpha=', alpha, '\n',
    'S=', S, '\n')

sim.data = dataset[[dataset.i]]$sim.data
  
x = sim.data$X[data.indices,]
y = sim.data$y[data.indices,]
  
deseq.results = estimate_signs(x,y,marginal=TRUE)
a.sgn = deseq.results$a.sgn
  
cov.matrix = dataset[[dataset.i]]$cov.matrix
  
S = dim(y)[2]
p = dim(x)[2]

mb_dat = list(S=S, N=n, p=p, noise=0.01, one_hot=a.sgn,
                X=x, y=y, Sigma=cov.matrix)

fit = sampling(st.model, data=mb_dat,
             iter=8000, chains=4, include=TRUE, 
             pars=c('a_global', 'alpha', 'a'), thin=10)

sampled.params = extract(fit, c('a_global', 'alpha', 'a'))
sampling.summary = summary(fit)$summary
  
if(!on.hec) {
  saveRDS(list(sampled.params=sampled.params, summary=sampling.summary, deseq.results=deseq.results), 
          file=here(paste0('Results/Simulation/Approximate/large_s_', S, '_alpha_', alpha, '_n_', n, '_dataset_', dataset.i, '.rds')))
} else {
  saveRDS(list(sampled.params=sampled.params, summary=sampling.summary, deseq.results=deseq.results), 
          file=paste0('/storage/users/dondelin/Work/Projects/IBDProject/Results/Simulation/Approximate/large_s_', S, '_alpha_', alpha, '_n_', n, '_dataset_', dataset.i, '.rds'))
}




