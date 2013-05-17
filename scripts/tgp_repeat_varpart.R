library(vegan)

setwd('~/tgp_management/')

source('./scripts/tgp_functions.R')

dir.create('./results/')

env = read.csv('./data/tgp_utm_env_complete.csv')
env = env[env$repeat_plot == 1, ]

comm = read.csv('./data/tgp_comm_mat_all.csv')
comm = comm[match(env$plot_yr, comm$plot.yr), ]

## drop species that don't occur
comm = comm[ , colSums(comm) > 0 ] 

row.names(comm) = comm$plot.yr
comm = comm[ , -1]
comm_sqr = sqrt(comm)

env$sr = rowSums(comm > 0)

## create dummy matrices
plot_id = sort(unique(env$plot))
year_id = sort(unique(env$yr))
plot_mat = matrix(0, ncol=length(plot_id), nrow=nrow(env))
year_mat = matrix(0, ncol=length(year_id), nrow=nrow(env))

for(i in 1:nrow(env)) {
  plot_mat[i, match(env$plot[i], plot_id)] = 1
  year_mat[i, match(env$yr[i], year_id)] = 1 
}

mang_vars = c('YrsOB', 'BP5Yrs', 'YrsSLB')
mang_mat = env[ , mang_vars]

## drop first columns
plot_mat = plot_mat[ , -1]
year_mat = year_mat[ , -1]

## carry out partitioning analysis
ols_part = ordi_part(env$sr, plot_mat, year_mat, mang_mat, method='rda')

rda_part = ordi_part(comm_sqr, plot_mat, year_mat, mang_mat, method='rda')

cca_part = ordi_part(comm_sqr, plot_mat, year_mat, mang_mat, method='cca',
                     nperm=999)

## print results
ols_part$part

rda_part$part

cca_part$part

save(ols_part, rda_part, cca_part, file='./results/repeat_varpart.Rdata')
