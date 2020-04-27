# Generate simulated microbiome data.

library(readr)
source('simulation.R')

set.seed(10)

# Try three-species model
n = 99
noise = 0.001
S = 5
p = 2
a.global = c(1.75, 0.5)

cov.matrix = generate_cov_matrix(S)

sim.data = generate_data_phylo_dirichlet(n, p, S, cov.matrix, a.global)

write_csv(as.data.frame(sim.data$X), paste('../Data/Simulation/mb', n, noise, S, p, 'covariates.csv', sep='_'))
write_csv(as.data.frame(sim.data$y), paste('../Data/Simulation/mb', n, noise, S, p, 'response.csv', sep='_'))
write_csv(as.data.frame(sim.data$a.sgn), paste('../Data/Simulation/mb', n, noise, S, p, 'signs.csv', sep='_'))
write_csv(as.data.frame(cov.matrix), paste('../Data/Simulation/mb', n, noise, S, p, 'covmat.csv', sep='_'))


