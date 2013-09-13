setwd('~/tgp_management/')

env = read.csv('./data/tgpall_env.csv')
env$date_samp = as.Date(env$date_samp, format='%m/%d/%Y')
env$date_of_burn = as.Date(env$date_of_burn, format='%m/%d/%Y')

## drop the plots not on the UTM grid
env = subset(env, subset=grid_plot == 1)

soilvars = c("P","CA","MG","K","NA.","B","FE","MN","CU","ZN","AL")
## add 1.0 to Boron so that there are no zero values
env$B = env$B + 1
logSoil = log10(env[ , soilvars])
names(logSoil) = paste('log', soilvars, sep='')

eastness = sin(env$aspect) 
northness = cos(env$aspect)

woody_vars = c('woodyht', 'woodypct', 'basal_area', 'density', 'avgdens')
env[ , woody_vars] = ifelse(is.na(as.matrix(env[ , woody_vars])), 
                            0, 
                            as.matrix(env[ , woody_vars]))

env = env[ , !(names(env) %in% soilvars)]
env = cbind(env, logSoil, eastness, northness)

## write out filtered env file
write.csv(env, file='./data/tgp_utm_env.csv', row.names=F)

plot(env$easting, env$northing)
plot(env$easting, env$northing, cex = log10(env$woodypct))

plot_id = sort(unique(env$plot))

plot_yr = sort(unique(env$plot_yr))

length(plot_id)

## read in species community data
comm = read.csv('./data/tgp_comm_mat_all.csv')

sp_plot_id = as.integer(
              sapply(strsplit(as.character(comm[,1]),'.', fixed=TRUE), 
                     function(x) x[[1]]))

sum(unique(sp_plot_id) %in% plot_id)
sum(comm[ ,1] %in% plot_yr)

## so every plot we have enviornmental data for we also have 
## species data

## generate plot year UTM output that will be used to query the 
## GIS layers

uni_plots = unique(env$plot)
out = data.frame(plot = uni_plots, 
                 easting = env$easting[match(uni_plots, env$plot)],
                 northing = env$northing[match(uni_plots, env$plot)])
write.csv(out, './data/tgp_plot_utm.csv', row.names=F)

### examine distance decay rel
comm = comm[match(plot_yr, comm[ , 1]), ]
library(vegan)
gdist = dist(env[ , c('easting', 'northing')])
vsim = 1 - vegdist(comm[ , -1])

plot(gdist, vsim, xlim=c(0, max(gdist)/2))
lines(lowess(gdist, vsim), col='red')