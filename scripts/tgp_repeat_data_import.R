
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

## site and year dummy variables
plot_id = sort(unique(env$plot))
year_id = sort(unique(env$yr))
plot_mat = matrix(0, ncol=length(plot_id), nrow=nrow(env))
year_mat = matrix(0, ncol=length(year_id), nrow=nrow(env))

for(i in 1:nrow(env)) {
  plot_mat[i, match(env$plot[i], plot_id)] = 1
  year_mat[i, match(env$yr[i], year_id)] = 1 
}

## drop first columns so no singular variables in models
plot_mat = plot_mat[ , -1]
year_mat = year_mat[ , -1]
