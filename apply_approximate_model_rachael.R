# Apply the full model to Rachael and Emma's dataset

library(rstan)
library(here)

source(here('Code/apply_DESeq2.R'))
source(here('Code/util.R'))

rstan_options(auto_write = TRUE)
#options(mc.cores = parallel::detectCores()-2)

set.seed(10)

on.hec = TRUE

tax_level = 'Rank5'

dataset = load_in_phyloseq_reanalysis()
subset_data = tax_glom(dataset$data, tax_level)
y = otu_table(subset_data)
x = t(ifelse(dataset$status=='Control', 0, 1))

st.model = stan_model(file=here('Code/phylo_dirichlet_multitask_model_fixed.stan'))

deseq.results = estimate_signs(x,y,marginal=TRUE)
a.sgn = deseq.results$a.sgn
  
cov.matrix = dataset[[dataset.i]]$cov.matrix
  
S = dim(y)[2]
p = dim(x)[2]

mb_dat = list(S=S, N=n, p=p, noise=0.01, one_hot=a.sgn,
                X=x, y=y, Sigma=cov.matrix)

fit = sampling(st.model, data=mb_dat,
             iter=8000, chains=4, include=TRUE, 
             pars=c('a_global', 'alpha', 'a'), thin=10, cores=1)

sampled.params = extract(fit, c('a_global', 'alpha', 'a'))
sampling.summary = summary(fit)$summary
  
if(!on.hec) {
  saveRDS(list(sampled.params=sampled.params, summary=sampling.summary, deseq.results=deseq.results), 
          file=here(paste0('Results/Simulation/Approximate/large_s_prop_05_', S, '_alpha_', alpha, '_n_', n, '_dataset_', dataset.i, '.rds')))
} else {
  saveRDS(list(sampled.params=sampled.params, summary=sampling.summary, deseq.results=deseq.results), 
          file=paste0('/storage/users/dondelin/Work/Projects/IBDProject/Results/Simulation/Approximate/large_s_prop_05_', S, '_alpha_', alpha, '_n_', n, '_dataset_', dataset.i, '.rds'))
}




