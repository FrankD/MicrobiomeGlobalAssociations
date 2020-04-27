# Apply the full model to simulated datasets with small S and small p

library(rstan)
library(here)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores()-2)

set.seed(10)

on.hec = FALSE

tax_level = 'Rank5'

dataset = load_in_phyloseq_reanalysis()
subset_data = tax_glom(dataset$data, tax_level)
y = t(otu_table(subset_data))
x = t(ifelse(dataset$status=='Control', 0, 1))
x = cbind(1, x[rownames(y),,drop=FALSE])

st.model = stan_model(file=here('Code/phylo_dirichlet_multitask_model_full.stan'))

cov.matrix = ape::vcv(phy_tree(subset_data))
  
S = dim(y)[2]
p = dim(x)[2]

k = 20
a = 4

mb_dat = list(S=S, N=dim(y)[1], p=p, noise=0.01, prior_mix_belief=c(a,k-2*a,a),
              X=x, y=y, Sigma=cov.matrix, tau=0.5)

fit = sampling(st.model, data=mb_dat,
             iter=8000, chains=4, include=TRUE, 
             pars=c('a_global', 'alpha', 'theta', 'a'))
  
if(!on.hec) {
  saveRDS(fit, file=here(paste0('Results/Simulation/Full/small_s_alpha_', alpha, '_n_', n, '_dataset_', dataset.i, '.rds')))
} else {
  saveRDS(fit, file=paste0('/storage/users/dondelin/Work/Projects/IBDProject/Results/Simulation/Full/small_s_alpha_', alpha, '_n_', n, '_dataset_', dataset.i, '.rds'))
}
  



