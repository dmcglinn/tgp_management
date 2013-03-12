library(sp)
library(rgdal)

setwd('~/Lab data/tgp_management/')

source('./scripts/tgp_functions.R')

load('./data/tgp_shpfiles.Rdata')

env = read.csv('./data/tgp_utm_env_complete.csv')

prj = proj4string(burns[[1]])

tgp_utm = SpatialPointsDataFrame(
         coords = coordinates(cbind(env$easting, env$northing)), 
         data = env, proj4string = CRS(prj))

plot(pasture)
points(env$easting, env$northing, pch=19, cex=.25)

plot(tgp_utm, pch=1)
par(bg='transparent')
par(new=TRUE)
plot(pasture)

path = getwd()
writeOGR(tgp_utm, dsn=file.path(path, 'tgpburn/coords'), 
         layer= 'tgp_utm', driver='ESRI Shapefile',
         dataset_options='RESIZE=no', overwrite_layer=TRUE)
