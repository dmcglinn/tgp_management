library(vegan)
library(sp)
library(nlme)


setwd('~/Lab data/tgp_management/')

source('./scripts/tgp_functions.R')

load('./data/tgp_shpfiles.Rdata')

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

nrow(env)
nrow(comm)

all.equal(env$plot_yr, comm$plot.yr)

row.names(comm) = comm$plot.yr
comm = comm[ , -1]
comm_sqr = sqrt(comm)

env$sr = rowSums(comm > 0)


prj = proj4string(pasture)
env_sp = SpatialPointsDataFrame(coords = cbind(env$easting, env$northing),
                                data=env, proj4string=CRS(prj))
                                  

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
soil_vars = c("P","CA","MG","K","NA","B","FE","MN","CU","ZN","AL")
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
soil_vars = c("P","CA","MG","K","NA","B","FE","MN","CU","ZN","AL")
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

soil_mat = as.data.frame(soil_pca$scores[ , 1:3])
mang_mat = env[ , mang_vars]

## basic question:
## how much variance in composition is related to 
## soil differences vs management differences

tgp_xy = env[ , c('easting', 'northing')]

## pca
pca = rda(comm_sqr)
pca_mso = mso(pca, tgp_xy)
## rda
rda_full = rda(comm_sqr, cbind(soil_mat, mang_mat))
plot(rda_full, display='bp')
rda_mso = mso(rda_full, tgp_xy)

par(mfrow=c(1,2))
msoplot(pca_mso, ylim=c(20, 40))
msoplot(rda_mso, ylim=c(20, 40))

rda_soil = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
               Condition(as.matrix(mang_mat)))
rda_mang = rda(comm_sqr ~ mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB +
               Condition(as.matrix(soil_mat)))

soil_tst = anova(rda_soil, step = 1000)
          Df     Var      F N.Perm Pr(>F)    
Model      3  2.3181 4.5891    999  0.001 ***
Residual 121 20.3736                         

soil_term_tst = anova(rda_soil, step=1000, by='margin')
                Df     Var      F N.Perm Pr(>F)    
soil_mat[, 1]   1  1.3312 7.9059    999  0.001 ***
soil_mat[, 2]   1  0.6598 3.9184    999  0.001 ***
soil_mat[, 3]   1  0.3375 2.0044    999  0.001 ***
Residual      121 20.3736                        

mang_tst = anova(rda_mang, step = 1000)
          Df     Var      F N.Perm Pr(>F)    
Model      3  1.2114 2.3214    999  0.001 ***
Residual 144 25.0491

mang_term_tst = anova(rda_mang, step=1000, by='margin')
                  Df     Var      F N.Perm Pr(>F)    
mang_mat$YrsOB    1  0.4307 2.5578    999 0.0010 ***
mang_mat$BP5Yrs   1  0.2095 1.2440    999 0.0820 .  
mang_mat$YrsSLB   1  0.2327 1.3817   1999 0.0335 *  
Residual        121 20.3736

varpart(comm_sqr, soil_mat, mang_mat)
[a] = X1|X2           3                 0.07683     TRUE
[b]                   0                 0.01494    FALSE
[c] = X2|X1           3                 0.02341     TRUE
[d] = Residuals                         0.88482    FALSE
## 8% for soil vs 2% for management
## same analysis but let's control for year effect
rda_yr = rda(comm_sqr ~ env$yr)
varpart(residuals(rda_yr), soil_mat, mang_mat)
## 7.6% for soil vs 1.3% for management

## ca
ca = cca(comm_sqr)
ca_mso = mso(ca, tgp_xy)
## cca
cca_full = cca(comm_sqr, cbind(soil_mat, mang_mat))
plot(cca_full, display=c('bp'))
cca_mso = mso(cca_full, tgp_xy)

par(mfrow=c(1,2))
msoplot(ca_mso, ylim=c(4, 10))
msoplot(cca_mso, ylim=c(4, 10))

cca_soil = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
               Condition(as.matrix(mang_mat)))
cca_mang = cca(comm_sqr ~ mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB +
               Condition(as.matrix(soil_mat)))

soil_tst = anova(cca_soil, step = 1000)
          Df  Chisq      F N.Perm Pr(>F)    
Model      3 0.3342 3.1180    999  0.001 ***
Residual 121 4.3228

soil_term_tst = anova(cca_soil, step=1000, by='margin')
                Df  Chisq      F N.Perm Pr(>F)    
soil_mat[, 1]   1 0.1649 4.6148    999  0.001 ***
soil_mat[, 2]   1 0.1075 3.0086    999  0.001 ***
soil_mat[, 3]   1 0.0630 1.7640    999  0.001 ***
Residual      121 4.3228

mang_tst = anova(cca_mang, step=1000)
          Df  Chisq      F N.Perm Pr(>F)    
Model      3 0.1781 1.6613    999  0.001 ***
Residual 121 4.3228

mang_term_tst = anova(cca_mang, step=1000, by='margin')
                  Df  Chisq      F N.Perm Pr(>F)    
mang_mat$YrsOB    1 0.0659 1.8451    999  0.001 ***
mang_mat$BP5Yrs   1 0.0405 1.1329    999  0.204    
mang_mat$YrsSLB   1 0.0511 1.4301    999  0.003 ** 
Residual        121 4.3228

full.r2<-r2_adj(comm_sqr,cbind(soil_mat, mang_mat),reps=1000,method='cca')
so.r2<-r2_adj(comm_sqr,soil_mat,reps=1000,method='cca')
ma.r2<-r2_adj(comm_sqr,mang_mat,reps=1000,method='cca')

par(mfrow=c(1,2))
hist(so.r2[-(1:2)])
hist(ma.r2[-(1:2)])

##fractions (legendre style)##much more complex looking results in almost identical results to the more straightforward palmer style
#soil | mang
full.r2[1] - ma.r2[1]
full.r2[2] - ma.r2[2]
#mang | soil
full.r2[1] - so.r2[1]
full.r2[2] - so.r2[2]
#soil + mang
full.r2[1] - (full.r2[1] - ma.r2[1]) - (full.r2[1] - so.r2[1])
full.r2[2] - (full.r2[2] - ma.r2[2]) - (full.r2[2] - so.r2[2])
#residuals
1 - full.r2[1]
1 - full.r2[2]

R2s<-cbind(c(full.r2[1],
             full.r2[1] - ma.r2[1],
             full.r2[1] - so.r2[1],
             full.r2[1] - (full.r2[1] - ma.r2[1]) - (full.r2[1] - so.r2[1])),
           c(full.r2[2],
             full.r2[2] - ma.r2[2],
             full.r2[2] - so.r2[2],
             full.r2[2] - (full.r2[2] - ma.r2[2]) - (full.r2[2] - so.r2[2]))
           )

colnames(R2s)<-c('R2','R2adj')
rownames(R2s)<-c('all','soil','mang','soil+mang')
round(R2s,3)
             R2 R2adj
all       0.111 0.065
soil      0.069 0.047
mang      0.037 0.014
soil+mang 0.006 0.004


## how does adding spatial predictors change our outcome

tgp_pcnm = pcnm(dist(tgp_xy))
par(mfrow=c(1,3))
ordisurf(tgp_xy, scores(tgp_pcnm, choi=1), bubble=4)
ordisurf(tgp_xy, scores(tgp_pcnm, choi=2), bubble=4)
ordisurf(tgp_xy, scores(tgp_pcnm, choi=3), bubble=4)

varpart(comm_sqr, soil_mat, mang_mat, scores(tgp_pcnm)[ , 1:3])
## not much change at all
varpart(comm_sqr, scores(tgp_pcnm)[ , 1:106], mang_mat)



rs = rowSums(comm) / sum(comm)
pcnmw = pcnm(dist(tgp_xy), w = rs)
ord = cca(comm ~ scores(pcnmw))
ord_mso = mso(ord, tgp_xy)
msoplot(ord_mso)
plot(ord)
par(mfrow=c(1,3))
ordisurf(tgp_xy, scores(pcnmw, choices=1), bubble = 3, main = "PCNM 1")
ordisurf(tgp_xy, scores(pcnmw, choices=2), bubble = 3, main = "PCNM 2")
ordisurf(tgp_xy, scores(pcnmw, choices=3), bubble = 3, main = "PCNM 3")

## Richness: model building ---------------------------------------------------------

mod1 = lm(sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
          YrsOB + BP5Yrs + YrsSLB, data=env)

mod2 = update(mod1, ~ . + slope + northness + eastness)
anova(mod1, mod2)
summary(mod2)

mod3 = update(mod1, ~ . + slope)
anova(mod3, mod2)
## only slope seems important for topography

mod4 = update(mod3, ~ . + woodypct)
anova(mod3, mod4)
summary(mod4)

mod_full = lm(sr ~ logCA + eastness + northness + slope + 
              YrsOB + YrsSLB + BP5Yrs + woodypct +
              avgdens + basal_area + rockpct,
              data = env)
              
mod_red = step(mod_full)
summary(mod_red)


mod5 = gls(sr ~ logCA + YrsOB + BP5Yrs +
            YrsSLB + slope + eastness,
            data=env)
mod5b = gls(sr ~ logCA + YrsOB + BP5Yrs +
             YrsSLB + slope + eastness,
             data=env, cor=corLin(form= ~ easting + northing))
anova(mod5, mod5b)

plot(Variogram(gls(sr ~ 1, data=env), form = ~ easting + northing)) 
plot(Variogram(mod5, resType='r'))
plot(Variogram(mod5b, resType='n')) 
par(mfrow=c(1,3))
plot(dist(tgp_xy), dist(env$sr))
lines(lowess(dist(tgp_xy), dist(env$sr)), lwd=2, col='red')
plot(dist(tgp_xy), dist(env$sr), ylim=c(12, 16))
lines(lowess(dist(tgp_xy), dist(env$sr)), lwd=2, col='red')
plot(lowess(dist(tgp_xy), dist(env$sr)), lwd=2, col='red')
abline(v=8e3)

plot(Variogram(gls(sr ~ 1, data=env), form = ~ easting + northing,
               maxDist=8e3),ylim = c(.8, 1)) 


## so there is very little spatial autocorrelation in richness,
## it is not necessary to use more complex models.

## examine partial effects for job talk
attach(env)
mod = lm(sr ~ logCA + slope + northness +
          YrsOB + BP5Yrs + YrsSLB)
termplot(mod, terms='YrsOB', partial=T, se=T,
         lwd.term=3,lwd.se=3,
         col.term='blue', col.se='blue',
         col.res = 'grey', col.smth = "red",
         frame.plot=F, axes=F, xlab='', ylab='',
         ylim=c(-30, 30))
axis(side=1, cex.axis=2, at=1:5)
axis(side=2, cex.axis=2)
res = residuals(mod, 'partial')
res = res[ , 'YrsOB', drop=FALSE]
lines(lowess(YrsOB, res), col='red', lty=2, lwd=3)

detach(env)
## Richness: variation partioning ---------------------------------------------------
varpart(env$sr, soil_mat, mang_mat)
[a] = X1|X2           3                 0.14920     TRUE
[b]                   0                -0.01084    FALSE
[c] = X2|X1           3                -0.00937     TRUE
[d] = Residuals                         0.87101    FALSE

varpart(env$sr, 
        cbind(env$logCA, env$slope, env$woodypct),
        mang_mat)
[a] = X1|X2           3                 0.25824     TRUE
[b]                   0                -0.02529    FALSE
[c] = X2|X1           3                 0.00508     TRUE
[d] = Residuals                         0.76197    FALSE

## final model to examine - b/c this was the model used
## for the repeat plots
varpart(env$sr, 
        cbind(env$logCA, env$slope, env$northness),
        mang_mat)
[a] = X1|X2           3                 0.21188     TRUE
[b]                   0                -0.02160    FALSE
[c] = X2|X1           3                 0.00139     TRUE
[d] = Residuals                         0.80832    FALSE

mod = lm(sr ~ logCA + slope + northness + YrsSLB + 
           YrsOB + BP5Yrs, data=env)
summary(mod)
Estimate Std. Error t value Pr(>|t|)    
(Intercept) 150.9279    18.3455   8.227 2.53e-13 ***
logCA       -24.3230     4.7970  -5.070 1.45e-06 ***
slope         0.5413     0.1564   3.461 0.000744 ***
northness     0.1798     1.4676   0.123 0.902683    
YrsSLB       -1.4955     1.0263  -1.457 0.147656    
YrsOB         0.7278     0.8023   0.907 0.366091    
BP5Yrs       -1.0704     1.2295  -0.871 0.385710
Multiple R-squared: 0.2299,  Adjusted R-squared: 0.1917 

## spatial sampling date examination ------------------------------------------------
library(sp)
library(dichromat)
plots = data.frame(plot=env$plot, yr = env$yr)
coordinates(plots) = cbind(env$easting, env$northing)
plots = plots[env$repeat_plot == 0, ]
plots@data$yr = as.factor(plots@data$yr)
cols = colorschemes$Categorical.12[1:5]
spplot(plots, 'yr', cex = 2, col.regions = cols)
