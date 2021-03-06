library(sp)
library(rgdal)
library(ggplot2)

setwd('~/Lab data/tgp_management/')

source('./scripts/tgp_functions.R')

## pasture folder
fencep = read.shape('fencep', path='./gis/pasture')
plot(fencep, axes=T)
summary(fencep)
head(fencep@data)
## this layer is UTM projected
## I assumed this is NAD 83

pasture_nad27 = read.shape('pasture_nad_27', path='./gis/pasture')
plot(pasture_nad27, axes=T)
summary(pasture_nad27)
head(pasture_nad27@data)
## in UTM coordinates
## name suggests its nad 27 projected layer
## this will have to be compared to other layers in nad27
## to verify

fencepasture = read.shape('fencepasture', path='./gis/pasture')
plot(fencePasture, axes=T)
summary(fencePasture)
head(fencePasture@data)
## this layer is Lat/long projected

pasture = read.shape( 'pasture', path='./gis/pasture')
plot(pasture, axes=T)
summary(pasture)
head(pasture@data)
## this layer is Lat/long projected
## has information on bison and cattle grazing
## transform this pasture layer to nad27 
prj = proj4string(read.shape('1990', path='./gis/burn'))
pasture = spTransform(pasture, CRS(prj))

xlims = c(7.2e5, 7.4e5)
ylims = c(4.068e6, 4.088e6)
plot(pasture, xlim=xlims, ylim=ylims, axes=T)
plot(fencep, xlim=xlims, ylim=ylims, add=T)
##zoom in
xlims = c(728798, 731413)
ylims = c(4081505, 4079789)
plot(pasture, xlim=xlims, ylim=ylims, axes=T)
plot(fencep, xlim=xlims, ylim=ylims, add=T)
##
burn00 = read.shape('2000', './gis/burn/')
xlims = c(7.2e5, 7.4e5)
ylims = c(4.068e6, 4.088e6)
plot(pasture_nad27, xlim=xlims, ylim=ylims, axes=T)
plot(burn00, xlim=xlims, ylim=ylims, lty=2, add=T)
##zoom in
xlims = c(728798, 731413)
ylims = c(4081505, 4079789)
plot(pasture_nad27, xlim=xlims, ylim=ylims, axes=T)
plot(burn00, xlim=xlims, ylim=ylims, lty=2, add=T)

## from these graphs we can see that the transformed layer
## pasture is truely in nad27 while the layer referred to as
## pasture_nad27 is mostly likely actually in nad28 or something
## very close

## tgp boundary----------------------------------------------------------------------
## TGP boundary
tgpBnd = read.shape('tgpBnd', path='./gis/tract')

ylims = c(4.070e6, 4.08e6)
xlims = c(7.240e5, 7.248e5)
plot(tgpBnd, axes=T, ylim=ylims, xlim=xlims)
plot(pasture, ylim=ylims, xlim=xlims, border='red', add=T)

## so the boundary is in NAD27
proj4string(tgpBnd) = CRS(proj4string(pasture))

## grazing folder--------------------------------------------------------------------

bison_unit_2008 = read.shape('bison_unit_2008', path='./gis/grazing')
plot(bison_unit_2008, col='green', xlim=c(7.2e5, 7.45e5), ylim=c(4.07e6, 4.082e6), 
     axes=T)
summary(bison_unit_2008)
head(bison_unit_2008@data)
## NAD83 zone 14 projection per email from Matt Poole

bison_update = read.shape('bisonupdate', path='./gis/grazing')
par(bg="transparent")
plot(bison_unit_2008, col='green', xlim=c(7.2e5, 7.45e5), ylim=c(4.07e6, 4.082e6),
     axes=T)
par(new=TRUE)
plot(bison_update, col="transparent", xlim=c(7.2e5, 7.45e5), ylim=c(4.07e6, 4.082e6))
## guessing NAD27 zone 14 projection b/c it appears slightly offset from
## the bison_unit_2008 shapefile

## let's check if it is actually nad27
xlims = c(7.2e5, 7.4e5)
ylims = c(4.068e6, 4.088e6)
plot(bison_update, xlim=xlims, ylim=ylims, axes=T)
plot(pasture, xlim=xlims, ylim=ylims, lty=2, add=T)
##zoom in
xlims = c(728798, 731413)
ylims = c(4081505, 4079789)
plot(bison_update, xlim=xlims, ylim=ylims, axes=T)
plot(pasture, xlim=xlims, ylim=ylims, lty=2, add=T)

## ok so its safe to associate the nad 27 projection info
## with this shapefile, let's rename it and do this
bison = bison_update
proj4string(bison) = CRS(proj4string(pasture))

## burn folder-----------------------------------------------------------------------
## burn layers are nad27 default

burn_files = dir('./gis/burn')
burn_shps = burn_files[grep('.shp', burn_files)]
burn_shps_nad27 = burn_shps[c(1:4, 6:9, 11:12, 14, 16, 18:19,
                               20, 22, 34, 26, 35, 32)] 
                              
burns = vector("list", length(burn_shps_nad27))
names(burns) = 1990:2009
prj = proj4string(read.shape('1990', path='./gis/burn'))

for (i in seq_along(burns)) {
  shpName = sub('.shp','', burn_shps_nad27[i])
  burns[[i]] = read.shape(shpName, path='./gis/burn')
  if (is.na(proj4string(burns[[i]]))) {
    proj4string(burns[[i]]) = CRS(prj)
  }  
}
names(burns)

## examine column names of burn shapefiles
sapply(burns, names)

sapply(burns, function(x) grep('DATE', names(x)))
## so in every case there is at least one column with the 
## word date in it, which should have the year/month/day of
## the burn

## consider fixing the 2005 burn date formate here so that it is 
## consistent with the other burn layers

## shapefile export------------------------------------------------------------------
## output burn shape files
save(pasture, bison, burns, tgpBnd,
     file='~/Lab data/tgp_management/data/tgp_shpfiles.Rdata')


