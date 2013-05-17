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

## create explanatory modeling variables
## soil variables
soil_vars = c("P","CA","MG","K","NA","B","FE","MN","CU","ZN","AL")
soil_vars = paste('log', soil_vars, sep='')
soil_pca = princomp(scale(env[ , soil_vars]))
soil_mat = as.data.frame(soil_pca$scores[ , 1:3])

## rain variables
env$rain2 = ifelse(is.na(env$rain2), 0, env$rain2)
sum_rain = apply(env[ , paste('rain', 6:9, sep='')], 1, sum)
win_rain = apply(env[ , paste('rain', c(10:12, 1), sep='')], 1, sum)
spr_rain = apply(env[ , paste('rain', 2:5, sep='')], 1, sum)
rain_mat = cbind(sum_rain, win_rain, spr_rain)

## management variables
mang_vars = c('YrsOB', 'BP5Yrs', 'YrsSLB')
mang_mat = env[ , mang_vars]

## carry out partitioning analysis
ols_part = ordi_part(env$sr, soil_mat, rain_mat, mang_mat, method='rda')

rda_part = ordi_part(comm_sqr, soil_mat, rain_mat, mang_mat, method='rda')

cca_part = ordi_part(comm_sqr, soil_mat, rain_mat, mang_mat, method='cca',
                     nperm=999)

## print results
ols_part$part

rda_part$part

cca_part$part

save(ols_part, rda_part, cca_part, file='./results/repeat_varpart.Rdata')
