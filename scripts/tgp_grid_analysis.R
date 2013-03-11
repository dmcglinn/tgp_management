library(vegan)
library(sp)

setwd('~/Lab data/tgp_management/')

source('./scripts/tgp_functions.R')

env = read.csv('./data/tgp_utm_env_complete.csv')
comm = read.csv('./data/tgp_comm_mat_all.csv')

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
Model      3  2.8348 5.4321    999  0.001 ***
Residual 144 25.0491

soil_term_tst = anova(rda_soil, step=1000, by='margin')
               Df     Var      F N.Perm Pr(>F)    
soil_mat[, 1]   1  1.6704 9.6029    999  0.001 ***
soil_mat[, 2]   1  0.8394 4.8254    999  0.001 ***
soil_mat[, 3]   1  0.3538 2.0340    999  0.004 ** 
Residual      144 25.0491                         

mang_tst = anova(rda_mang, step = 1000)
          Df     Var      F N.Perm Pr(>F)    
Model      3  1.2114 2.3214    999  0.001 ***
Residual 144 25.0491

mang_term_tst = anova(rda_mang, step=1000, by='margin')
                 Df     Var      F N.Perm Pr(>F)    
mang_mat$YrsOB    1  0.4361 2.5073    999  0.002 ** 
mang_mat$BP5Yrs   1  0.4346 2.4981    999  0.001 ***
mang_mat$YrsSLB   1  0.3107 1.7861    999  0.006 ** 
Residual        144 25.0491

varpart(comm_sqr, soil_mat, mang_mat)
[a] = X1|X2           3                 0.08003     TRUE
[b]                   0                 0.01133    FALSE
[c] = X2|X1           3                 0.02386     TRUE
[d] = Residuals                         0.88478    FALSE
## 8% for soil vs 2% for management
## same analysis but let's control for year effect
rda_yr = rda(comm_sqr ~ env$yr)
varpart(residuals(rda_yr), soil_mat, mang_mat)
## 7.8% for soil vs 1.7% for management

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
               env$eastness + env$northness + 
               Condition(as.matrix(mang_mat)))
cca_mang = cca(comm_sqr ~ mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB +
               Condition(as.matrix(soil_mat)))

soil_tst = anova(cca_soil, step = 1000)
          Df  Chisq      F N.Perm Pr(>F)    
Model      3 0.3897 3.3705    999  0.001 ***
Residual 144 5.5497

soil_term_tst = anova(cca_soil, step=1000, by='margin')
               Df  Chisq      F N.Perm Pr(>F)    
soil_mat[, 1]   1 0.2001 5.1916    999  0.001 ***
soil_mat[, 2]   1 0.1276 3.3113    999  0.001 ***
soil_mat[, 3]   1 0.0642 1.6658    999  0.015 *  
Residual      144 5.5497

mang_tst = anova(cca_mang, step=1000)
          Df  Chisq      F N.Perm Pr(>F)    
Model      3 0.1960 1.6951    999  0.001 ***
Residual 144 5.5497

mang_term_tst = anova(cca_mang, step=1000, by='margin')
                 Df  Chisq      F N.Perm Pr(>F)   
mang_mat$YrsOB    1 0.0678 1.7595    999  0.004 **
mang_mat$BP5Yrs   1 0.0700 1.8164    999  0.002 **
mang_mat$YrsSLB   1 0.0590 1.5313    999  0.005 **
Residual        144 5.5497

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
            R2  R2adj
all       0.100 0.062
soil      0.063 0.045
mang      0.032 0.013
soil+mang 0.005 0.004



## how does adding spatial predictors change our outcome

tgp_pcnm = pcnm(dist(tgp_xy))
varpart(comm_sqr, soil_mat, mang_mat, scores(tgp_pcnm)[ , 1:3])
## not much change at all


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



## spatial sampling date examination ------------------------------------------------
library(sp)
library(dichromat)
plots = data.frame(plot=env$plot, yr = env$yr)
coordinates(plots) = cbind(env$easting, env$northing)
plots = plots[env$repeat_plot == 0, ]
plots@data$yr = as.factor(plots@data$yr)
cols = colorschemes$Categorical.12[1:5]
spplot(plots, 'yr', cex = 2, col.regions = cols)
