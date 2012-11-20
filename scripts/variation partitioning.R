
spdat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)

endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env.csv',sep=',',header=T)


names(endat)
 [1] "plot"      "yr"        "plotnum"   "yr.plot"   "X205"      "X206"     
 [7] "X208"      "X220"      "X222"      "X226"      "X238"      "X244"     
[13] "X254"      "X259"      "X303"      "X307"      "X308"      "X309"     
[19] "X317"      "X319"      "X331"      "X343"      "X346"      "X350"     
[25] "X1998"     "X1999"     "X2000"     "X2001"     "X2002"     "X2003"    
[31] "X2004"     "X2005"     "X2006"     "X2007"     "X2008"     "OrBis"    
[37] "ChBis"     "bison"     "cattle"    "YrsOB"     "BP5Yrs"    "YrsSLB"   
[43] "burn"      "spring"    "summer"    "fall"      "slope"     "northness"
[49] "eastness"  "OM"        "pH"        "SolS"      "logP"      "logCa"    
[55] "logFe"     "logMg"     "logK"      "logNa"     "logB"      "logMn"    
[61] "logCu"     "logZn"     "logAl"     "PDSIavg"   "SodTemp"   "SPI1"     
[67] "SPI12"     "SPI24"

library(vegan)

sp<-as.matrix(sqrt(spdat[,-1]))
site<-as.matrix(endat[,5:24])
yr<-as.matrix(endat[,25:35])

sp.dw<-downweight(sp,frac=5)

#based on site pca
#site.vars is soil and topo
site.vars<-as.matrix(endat[,48:63])
colnames(site.vars)
site.vars<-site.vars[,-c(6,9,14,16)]
rownames(dist)<-endat[,3]

#based on dist pca
dist<-as.matrix(cbind(endat[,38:42],endat[,44]))
colnames(dist)<-c(colnames(endat)[38:42],colnames(endat)[44])
rownames(dist)<-endat[,3]

#based on clim pca
clim<-as.matrix(endat[,64:68])
rownames(clim)<-endat[,3]


###DCA###################
dca<-decorana(sp.dw)
names<-colnames(sp.dw)
pl<-plot(dca,dis='sp')
identify(pl,'sp',labels=names)
###############################
###############################
(cca.mike<-cca(sp.dw~spring+cattle+YrsOB+BP5Yrs+YrsSLB+Condition(site)+Condition(yr),data=endat))
              Inertia Rank
Total         1.47585     
Conditional   0.85788   20
Constrained   0.02795    5
Unconstrained 0.59002  194
anova.cca(cca.mike,strata=rep(1:20,11),by='terms')
anova.cca(cca.mike,strata=endat[,2],step=1000)
anova.cca(cca.mike,by='terms')
plot(cca.mike)
ordiplot3d(cca.mike)
###############################
###############################
##all vars
(cca.all<-cca(sp.dw~site.vars+clim+dist))
names<-colnames(sp.dw)
pl<-plot(cca.all,dis=c('sp','bp'))
identify(pl,'sp',labels=names,cex=.75)

##clim only
(cca.clim<-cca(sp.dw~clim))
pl<-plot(cca.clim,dis=c('sp','bp'))
identify(pl,'sp',labels=names,cex=.75)


##try to expl between site variation##
(cca.yr.site.vars<-cca(sp.dw~site.vars+Condition(dist)+Condition(yr)))
              Inertia Rank
Total          1.4758     
Conditional    0.2187   15
Constrained    0.4292   12 ##29.1%
Unconstrained  0.8279  192
##try to expl within site variation##
(cca.site.site.vars<-cca(sp.dw~site.vars+Condition(clim)+Condition(dist)+Condition(site)))
              Inertia Rank
Total         1.47585     
Conditional   0.90940   29
Constrained   0.04089   10 ##2.8%
Unconstrained 0.52555  180
#
(cca.site.clim<-cca(sp.dw~clim+Condition(site.vars)+Condition(dist)+Condition(site)))
              Inertia Rank
Total         1.47585     
Conditional   0.92116   34
Constrained   0.02914    5 ##2.0%
Unconstrained 0.52555  180
#
(cca.site.dist<-cca(sp.dw~dist+Condition(site.vars)+Condition(clim)+Condition(site)))
              Inertia Rank
Total         1.47585     
Conditional   0.92334   34
Constrained   0.02696    5 ##1.8%
Unconstrained 0.52555  180
#
(cca.site.dist.site.vars<-cca(sp.dw~dist+site.vars+Condition(clim)+Condition(site)))
              Inertia Rank
Total          1.4759     
Conditional    0.8821   24
Constrained    0.0682   15 
Unconstrained  0.5255  180
# shared var between site.vars and dist approx. 0
(cca.site.dist.clim<-cca(sp.dw~dist+clim+Condition(site.vars)+Condition(site)))
             Inertia Rank
Total         1.47585     
Conditional   0.89313   29
Constrained   0.05717   10
Unconstrained 0.52555  180
# shared var between dist and clim approx. 0
(cca.site.clim.site.vars<-cca(sp.dw~clim+site.vars+Condition(dist)+Condition(site)))
              Inertia Rank
Total          1.4759     
Conditional    0.8757   24
Constrained    0.0746   15
Unconstrained  0.5255  180
# shared var between clim and site.vars approx. 0

###############################
###############################
(rda.full<-rda(sp.dw,as.matrix(cbind(site,yr))))
              Inertia Rank
Total           25.50     
Constrained     15.45   29
Unconstrained   10.05  190
part<-varpart(sp.dw,site,yr)
plot(part)
Partition table:
                     Df R.squared Adj.R.squared Testable
[a+b] = X1           19   0.53229       0.48786     TRUE
[b+c] = X2           10   0.07304       0.02868     TRUE
[a+b+c] = X1+X2      29   0.60533       0.54509     TRUE
Individual fractions                                    
[a] = X1|X2          19                 0.51641     TRUE
[b]                   0                -0.02855    FALSE
[c] = X2|X1          10                 0.05723     TRUE
[d] = Residuals                         0.45491    FALSE

anova(rda(sp,yr,site),step=200)
anova(rda(sp,yr,site),step=200,strata=as.factor(endat[,1])) ##within sites
anova(rda(sp,site,yr),step=200,strata=as.factor(endat[,2])) ##within yrs
##both are significicant

##CCA##############################
(cca.full<-cca(sp.dw,as.matrix(cbind(site,yr))))
              Inertia Rank
Total          1.4758     
Constrained    0.9133   29
Unconstrained  0.5625  190


(pcca.site<-cca(sp.dw,site,yr))
Conditional   0.07052   10 ##4.8% of tot interia expl by yr alone
Constrained   0.84280   19 ##57.1% of tot interia expl by site alone
##shared variation is 0

anova(pcca.site,strata=as.factor(endat[,2]))
anova(pcca.yr,strata=as.factor(endat[,1]))



dist<-as.matrix(endat[,38:41])
colnames(dist)
dist[,1]<-as.factor(ifelse(dist[,1]==1,'cattle','bison'))

soil<-as.matrix(endat[,46:62])

##replace soil with pca axis 1 to 4, so equal # of vars between dist and mang
soil<-site.pca$x[,1:4]

##Investigate between site variation
(pcca.site.full<-cca(sp.dw~dist+soil+Condition(yr)))
0.65196/tot #44.2%
(pcca.site.soil<-cca(sp.dw~soil+Condition(dist)+Condition(yr)))
0.5123/tot #34.7%
(pcca.site.dist<-cca(sp.dw~dist+Condition(soil)+Condition(yr)))
0.04653/tot #3.1%
.65196-.5123-.05653 #8.3% shared var

##Investigate within site variation
(pcca.yr.full<-cca(sp.dw~dist+soil+Condition(site)))
 0.0903/tot #6.1%
(pcca.yr.soil<-cca(sp.dw~soil+Condition(dist)+Condition(site)))
 0.06128/tot #4.2%
(pcca.yr.dist<-cca(sp.dw~dist+Condition(soil)+Condition(site)))
0.02313/tot #1.6%
0.0903-0.06128-0.02313 #0.6% shared var

names<-colnames(sp.dw)
pl<-plot(pcca.site.soil,dis=c('sp','bp'),main='between sites')
identify(pl,'sp',labels=names,cex=.75)

pl<-plot(pcca.yr.soil,dis=c('sp','bp'),main='within sites')
identify(pl,'sp',labels=names,cex=.75)

################################################################
################################################################
##adjusted CCA fractions####
##Peres-Neto et al. 2006####
(1) randomly permute entire rows of data matrix X (i.e., no substantial difference was
found if regressors were permuted separately), leading to Xperm; 
(2) Calculate R2 CCAjX for a CCA based on Xperm;
(3) repeat steps 1 and 2 m times (in this study we used m 1000);
(4) calculate the mean Xperm across all 1000 R2 

m<-20
r2.adj<-function(m,Y,X,method){
 rand.r2<-rep(NA,m)
 if(method=='cca'){
  for(i in 1:m){
   Xrand<-site[sample(dim(X)[1]),]
   rand.r2[i]<-summary( cca(Y,Xrand))$constr.chi
  }
  cca.emp<-cca(Y,X)
  r2<-summary(cca.site)$constr.chi/cca.site$tot.chi 
 }
 1-(1/(1-mean(rand.r2)))*(1-r2)
}
r2.adj(20,sp.dw,site,'cca')

cca.site<-cca(sp.dw,site)
R2<-summary(cca.site)$constr.chi/cca.site$tot.chi
1-(1/(1-mean(rand.r2)))*(1-R2)

##for pcca##






