library(vegan)

setwd('~/Lab data/tgp_management/')

env = read.csv('./data/tgp_utm_env_complete.csv')
comm = read.csv('./data/tgp_comm_mat_all.csv')

pl_yr = env$plot_yr[env$repeat_plot == 0]
pl_yr = sort(c(pl_yr, env$plot_yr[env$repeat_plot == 1 & env$yr == 1998]))

env = env[match(pl_yr, env$plot_yr), ]

comm = comm[match(env$plot_yr, comm$plot.yr), ]

nrow(env)
nrow(comm)

all.equal(env$plot_yr, comm$plot.yr)

row.names(comm) = comm$plot.yr
comm = comm[ , -1]
comm_sqr = sqrt(comm)

## define site vars 

## examine woodyness variables-------------------------------------------------------
woody_vars = c('woodyht', 'woodypct', 'basal_area', 'density', 'avgdens')
woody_pca = princomp(scale(env[ , woody_vars]))
summary(woody_pca)
sum(woody_pca$sdev)
par(mfrow=c(1,2))
plot(woody_pca)
biplot(woody_pca)
## so if woodyness is to be used as an explanatory variable
## it appears that woodypct is the best variable for the job

## examine soil variables------------------------------------------------------------
soil_vars = c("P","CA","MG","K","NA.","B","FE","MN","CU","ZN","AL")
soil_vars = paste('log', soil_vars, sep='')
soil_pca = princomp(scale(env[ , soil_vars]))
summary(soil_pca)
sum(soil_pca$sdev)
par(mfrow=c(1,2))
plot(soil_pca)
biplot(soil_pca)
par(mfrow=c(1,1))
biplot(soil_pca, col=c('white', 'black'))
biplot(soil_pca, choices = c(3,4), col=c('white', 'black'))

soil_sc = data.frame(soil_pca$scores, ca = env$logCA, ph=env$PH)
coordinates(soil_sc) = env[ , c('easting', 'northing')]
spplot(soil_sc, 'Comp.1', col.regions=rev(terrain.colors(6))[-1],
       cex=3.25, pch=rep(15, 5))

## 1st axis is predominately pH/CA/Fe gradient
## 2nd axis is MG / P gradient - not as interpretable
## 4 axes are required to get 77 % of the variance

## examine site variables------------------------------------------------------------
soil_vars = c("P","CA","MG","K","NA.","B","FE","MN","CU","ZN","AL")
soil_vars = paste('log', soil_vars, sep='')
site_vars = c(soil_vars, 'eastness', 'northness', 'aspect')
site_pca = princomp(scale(env[ , site_vars]))
summary(site_pca)
sum(site_pca$sdev)
par(mfrow=c(1,2))
plot(site_pca)
biplot(site_pca)
par(mfrow=c(1,1))
biplot(site_pca, col=c('white', 'black'))
biplot(site_pca, choices = c(3,4), col=c('white', 'black'))

## the topographic variables are too similar between plots to 
## dominate much of the site differences

## examine management variables -----------------------------------------------------
## examine soil variables
mang_vars = c('YrsOB', 'BP5Yrs', 'YrsSLB')
mang_pca = princomp(scale(env[ , mang_vars]))
summary(mang_pca)
sum(mang_pca$sdev)
par(mfrow=c(1,2))
plot(mang_pca)
biplot(mang_pca)
par(mfrow=c(1,1))
biplot(mang_pca, col=c('white', 'black'))

mang_sc = data.frame(mang_pca$scores, env[ , mang_vars])
coordinates(mang_sc) = env[ , c('easting', 'northing')]
spplot(mang_sc, 'YrsSLB', col.regions=rev(terrain.colors(6))[-1],
       cex=3.25, pch=rep(15, 5))

## these maps indicate how the management variables will act
## to capture coarser scale regional spatial variation between
## plots thus also acting to large degree as spatial predictiors

soil_mat = soil_pca$scores[ , 1:3]
mang_mat = as.matrix(env[ , mang_vars])

## basic question:
## how much variance in composition is related to 
## soil differences vs management differences

## pca
pca = rda(comm_sqr)
pca_mso = mso(pca, env[ , c('easting', 'northing')])
## rda
rda_full = rda(comm_sqr, cbind(soil_mat, mang_mat))
plot(rda_full)
rda_mso = mso(rda_full, env[ , c('easting', 'northing')])

par(mfrow=c(1,2))
msoplot(pca_mso, ylim=c(20, 40))
msoplot(rda_mso, ylim=c(20, 40))

varpart(comm_sqr, soil_mat, mang_mat)
## 8% for soil vs 2% for management

## cca



## spatial sampling date examination ------------------------------------------------
library(sp)
library(dichromat)
plots = data.frame(plot=env$plot, yr = env$yr)
coordinates(plots) = cbind(env$easting, env$northing)
plots = plots[env$repeat_plot == 0, ]
plots@data$yr = as.factor(plots@data$yr)
cols = colorschemes$Categorical.12[1:5]
spplot(plots, 'yr', cex = 2, col.regions = cols)
