library(vegan)
library(sp)

r2_adj<-function(Y,X,Z,reps,method,dummy=0) {
  ##purpose: returns R2, R2adj, and all R2adj replicates that result from permuations
  ##Y,X,Z are spdata, expl mat, covar mat
  ##reps the number of permutations to perform, if no reps then only r2 and r2adj returned 
  ##if reps is not provided then it calculates an analytical R2adj (only ok for RDA)
  ##method specifies "cca" or "rda"
  ##dummy is a number 0, 1 or 2 depending on how many collinear variables are in the explanatory matrix
  ##dummy is only necessary for the analytical R2adj calculation
  Y<-as.matrix(Y)
  X<-as.matrix(X)
  if(missing(Z)){
    cca.emp<-eval(parse(text= paste(method,'(Y,X)')))
    r2<-summary(cca.emp)$constr.chi/cca.emp$tot.chi 
    if(missing(reps)){
      n<-nrow(Y)
      p<-ncol(X)-dummy
      out = c(r2, 1-(((n-1)/(n-p-1))*(1-r2)))
    }
    else{
      rand.r2<-rep(NA,reps)
      for(i in 1:reps){
        Xrand<-X[sample(nrow(X)),]
        rand.r2[i]<-summary( eval(parse(text= paste(method,'(Y,Xrand)'))))$constr.chi
        print(i)
      }
      out = c(r2, 
              1-(1/(1-mean(rand.r2/cca.emp$tot.chi)))*(1-r2),
              1-(1/(1-rand.r2/cca.emp$tot.chi))*(1-r2))
    }
  }  
  else{
    Z<-as.matrix(Z)
    cca.emp<-eval(parse(text= paste(method,'(Y,X,Z)')))
    r2<-summary(cca.emp)$constr.chi/cca.emp$tot.chi
    if(missing(reps)){
      n<-nrow(Y)
      p<-ncol(X)-dummy
      out = c(r2,1-(((n-1)/(n-p-1))*(1-r2)))
    }
    else{
      rand.r2<-rep(NA,reps)
      for(i in 1:reps){
        rhold<-sample(nrow(X))
        Xrand<-X[rhold,]
        Zrand<-Z[rhold,]
        rand.r2[i]<-summary( eval(parse(text= paste(method,'(Y,Xrand,Zrand)'))))$constr.chi
        print(i)
      }
      out = c(r2,
              1-(1/(1-mean(rand.r2/cca.emp$tot.chi)))*(1-r2),
              1-(1/(1-rand.r2/cca.emp$tot.chi))*(1-r2))
    }  
  }
  return(out)
}

setwd('~/Lab data/tgp_management/')

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

tgp_xy = env[ , c('easting', 'northing')]

## pca
pca = rda(comm_sqr)
pca_mso = mso(pca, tgp_xy)
## rda
rda_full = rda(comm_sqr, cbind(soil_mat, mang_mat))
plot(rda_full)
rda_mso = mso(rda_full, tgp_xy)

par(mfrow=c(1,2))
msoplot(pca_mso, ylim=c(20, 40))
msoplot(rda_mso, ylim=c(20, 40))

varpart(comm_sqr, soil_mat, mang_mat)
## 8% for soil vs 2% for management

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
