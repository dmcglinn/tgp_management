
env = read.csv('./data/tgp_utm_env_complete.csv')
comm = read.csv('./data/tgp_comm_mat_all.csv')

## define grassland plots
env[is.na(env$waterpct), ]$waterpct = 0 ## checked coordinate
env[is.na(env$rockpct), ]$rockpct = 0

grassland = env$waterpct == 0 &
            env$rockpct <= 20 &
            env$woodypct <= 20
env = env[grassland, ]

pl_yr = env$plot_yr[env$repeat_plot == 0]
pl_yr = sort(c(pl_yr, env$plot_yr[env$repeat_plot == 1 & env$yr == 1998]))

env = env[match(pl_yr, env$plot_yr), ]

comm = comm[match(env$plot_yr, comm$plot.yr), ]
## drop species that don't occur
comm = comm[ , apply(comm, 2, sum) > 0 ] 
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

## management variables
mang_vars = c('YrsOB', 'BP5Yrs', 'YrsSLB')
mang_mat = env[ , mang_vars]
