setwd('~/Lab data/tgp_management/')

env = read.csv('./data/tgpall_env.csv')
## drop the plots not on the UTM grid
env = subset(env, subset=grid_plot == 1)

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