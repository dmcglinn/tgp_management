library(vegan)

source('./scripts/tgp_functions.R')
source('./scripts/tgp_repeat_data_import.R')

dir.create('./results/')

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


