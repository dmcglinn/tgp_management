spdat<-read.table('C:/Users/dmcglinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
endat<-read.table('C:/Users/dmcglinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
#spdat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
#endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
attach(endat)
##begin working only at the 100 m2 scale so no corners
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)

library(vegan)

##prepare variables
spdat.n<-spdat[,-1]
rownames(spdat.n)<-spdat[,1]
##sqrt transform and downweight rare species species of sp mat 
spdat.n<-downweight(sqrt(spdat.n))

plot(colSums(sqrt(spdat[,-1]))[order(colSums(sqrt(spdat[,-1])),decreasing=T)],type='l',lwd=2)
points(colSums(spdat.n)[order(colSums(sqrt(spdat[,-1])),decreasing=T)],col='red',type='l')


yr.f<-factor(yr)
site.m<-endatlv1[,8:27]
yr.m<-endatlv1[,28:38]
dist.m<-cbind(YrsOB,YrsSLB,BP5Yrs)

uni.yrs<-unique(yr)
trues<-matrix(0,ncol=length(uni.yrs),nrow=length(yr))
trues<-t(sapply(uni.yrs,function(x)apply(spdat.n[yr==x,],2,sum)>0))
par(mfrow=c(1,2))
plot(uni.yrs,sapply(1:length(uni.yrs),function(x)cca(spdat.n[yr==uni.yrs[x],trues[x,]])$tot.chi),type='o',xlab='year',ylab='community variability',main='Chi sqr. intertia')
plot(uni.yrs,sapply(1:length(uni.yrs),function(x)rda(spdat.n[yr==uni.yrs[x],trues[x,]])$tot.chi),type='o',xlab='year',ylab='community variability',main='Euclid intertia (i.e., Total Variance)')
plot(uni.yrs,sapply(1:length(uni.yrs),function(x)rda(spdat.n[yr==uni.yrs[x],])$tot.chi),type='o',xlab='year',ylab='community variability')
##turns out including species with no abundance doesn't change anything, duh


##DCA
tgp.dca<-decorana(spdat.n)
tgp.dca
Call:
decorana(veg = spdat.n) 
Detrended correspondence analysis with 26 segments.
Rescaling of axes with 4 iterations.
Downweighting of rare species from fraction 1/5.

                  DCA1   DCA2    DCA3    DCA4
Eigenvalues     0.1519 0.1083 0.08495 0.06401
Decorana values 0.1586 0.1084 0.07202 0.04639
Axis lengths    2.2343 1.8212 1.53508 1.16142
##the eigenvalues are slightly different than produced with CANOCO
##axis lengths match though - these are lengths of gradient in canoco

sit.sco<-scores(tgp.dca)

##import CANOCO sol file
dca.sites<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/DCA/dcasiteresult.csv',sep=',',header=T)
names(dca.sites)
[1] "N"      "Site"   "AX1"    "AX2"    "AX3"    "AX4"    "WEIGHT" "N2"    

plot(sit.sco,pch=1,xlab='',ylab='',xlim=c(-1,2),ylim=c(-1,2))
points(dca.sites[,3]-1.05,dca.sites[,4]-.7,pch=19)
##but the plots between CANOCO and vegan look identical


##ploting Figure 3 for paper
x<-1 #axis1 #
y<-2 #axis2 #
##ploting Figure 3 for paper
plot(sit.sco[,x],sit.sco[,y],type='n',xlab='',ylab='',cex.axis=1.5)
uni.plots<-unique(name)
for(i in 1:20){
 bsum<-sum(bison[name==uni.plots[i]])
 if(bsum>0){
  ##this algo works b/c once a plot becomes bison it stays bison
  if(bsum<11){
   points(sit.sco[name==uni.plots[i],x][1:(12-bsum)],sit.sco[name==uni.plots[i],y][1:(12-bsum)],type='l',col='grey',lwd=3)
  }    
  points(sit.sco[name==uni.plots[i],x][(12-bsum):11],sit.sco[name==uni.plots[i],y][(12-bsum):11],type='l',col='black',lwd=3)
 }
 else{
  points(sit.sco[name==uni.plots[i],x],sit.sco[name==uni.plots[i],y],type='l',col='grey',lwd=3)
}}
points(sit.sco[yr==1998,x],sit.sco[yr==1998,y],pch=19,cex=1.25)
points(sit.sco[yr==2008,x],sit.sco[yr==2008,y],pch=0,cex=1.25)
#abline(h=0,lty=2,col='grey')
#abline(v=0,lty=2,col='grey')

##ploting Figure 3 for Talk in color
plot(sit.sco,type='n',xlab='',ylab='',cex.axis=1.5)
uni.plots<-unique(name)
for(i in 1:20){
 bsum<-sum(bison[name==uni.plots[i]])
 if(bsum>0){
  ##this algo works b/c once a plot becomes bison it stays bison
  if(bsum<11){
   points(sit.sco[name==uni.plots[i],1:2][1:(12-bsum),],type='l',col='pink3',lwd=3)
  }    
  points(sit.sco[name==uni.plots[i],1:2][(12-bsum):11,],type='l',col='lightblue',lwd=3)
 }
 else{
  points(sit.sco[name==uni.plots[i],],type='l',col='pink3',lwd=3)
}}
points(sit.sco[yr==1998,1:2],pch=19,cex=1.25)
points(sit.sco[yr==2008,1:2],pch=0,cex=1.25)
plot(1:10,1:10,type='n')
legend('center',c('',''),col=c('pink3','lightblue'),lwd=3,bty='n')


##test out ANOSIM
tgp.dist<-vegdist(spdat.n)
##basing the dist measures off non-sqrt transformed covers 
##ends up showing a lot more within plot variance in composition

tgp.ano<-anosim(tgp.dist,name)
plot(tgp.ano)
tgp.ano

##test out nmds
library(MASS)
tgp.mds0 <- isoMDS(tgp.dist)
stressplot(tgp.mds0, tgp.dist)
sit.sco<-scores(tgp.mds0)
plot(sit.sco,type='n',xlab='',ylab='')
uni.plots<-unique(name)
for(i in 1:20){
 bsum<-sum(bison[name==uni.plots[i]])
 if(bsum>0){
  ##this algo works b/c once a plot becomes bison it stays bison
  if(bsum<11){
   points(sit.sco[name==uni.plots[i],1:2][1:(12-bsum),],type='l',col='grey',lwd=2)
  }    
  points(sit.sco[name==uni.plots[i],1:2][(12-bsum):11,],type='l',col='black',lwd=2)
 }
 else{
  points(sit.sco[name==uni.plots[i],],type='l',col='grey',lwd=2)
}}
points(sit.sco[yr==1998,1:2],pch=19)
points(sit.sco[yr==2008,1:2],pch=0)

tgp.nmds1<-metaMDS(spdat[,-1], trace = FALSE);alarm()
names(tgp.nmds1)
stressplot(tgp.nmds1)
plot(scores(tgp.nmds1))


##first use raw data approach, cca and rda
##start with rda
rda1<-rda(spdat.n~name+yr.f)
rda1
test<-anova(rda1,by="terms")
test
plot(rda1,type='n')
text(rda1, dis="cn",cex=.5)

pl<-plot(rda2,dis=c('sp','bp','cn'),label="F")
identify(pl,'sp',labels=colnames(spdat.n),cex=.75))
ordiplot3d(rda2)

##varpart requires that the factor variables are matrices
##will only work for RDA
varpart(spdat.n,site.m,yr.m,dist.m)
plot(varpart(spdat.n,yr.m,site.m,dist.m))

###########################################################
cca.full<-cca((spdat.n)~name+yr.f+YrsOB+YrsSLB+BP5Yrs)
cca.full
lcca.full <- as.mlm(cca.full)
## Influential observations
influence.measures(lmod)
plot(cca.full, type = "n")
points(cca.full, cex = 10*hatvalues(lcca.full), pch=16, xpd = TRUE)
text(cca.full, display = "cn", col = "blue",cex=.75) 
##
mso.full<-mso(cca.full,yr)
msoplot(mso.full)
##
ca<-cca(spdat.n)
pca<-rda(spdat.n)
mso.ca<-mso(ca,yr,permutations=300)
mso.pca<-mso(pca,yr,permutations=300)
msoplot(mso.ca)
msoplot(mso.pca)
##
################################################################
################################################################
##adjusted CCA fractions####
##Peres-Neto et al. 2006####
(1) randomly permute entire rows of data matrix X (i.e., no substantial difference was
found if regressors were permuted separately), leading to Xperm; 
(2) Calculate R2 CCAjX for a CCA based on Xperm;
(3) repeat steps 1 and 2 m times (in this study we used m 1000);
(4) calculate the mean Xperm across all 1000 R2 

r2.adj<-function(Y,X,Z,reps){
 rand.r2<-rep(NA,reps)
 Y<-as.matrix(Y)
 X<-as.matrix(X)
 if(missing(Z)){
  for(i in 1:reps){
   Xrand<-X[sample(nrow(X)),]
   rand.r2[i]<-summary( cca(Y,Xrand))$constr.chi
   print(i)
  }
  cca.emp<-cca(Y,X)
  r2<-summary(cca.emp)$constr.chi/cca.emp$tot.chi 
  c(r2,1-(1/(1-mean(rand.r2/cca.emp$tot.chi)))*(1-r2))
 }
 else{
  Z<-as.matrix(Z)
  for(i in 1:reps){
   rhold<-sample(nrow(X))
   Xrand<-X[rhold,]
   Zrand<-Z[rhold,]
   rand.r2[i]<-summary( cca(Y,Xrand,Zrand))$constr.chi
   print(i)
  }
  cca.emp<-cca(Y,X,Z)
  r2<-summary(cca.emp)$constr.chi/cca.emp$tot.chi 
  c(r2,1-(1/(1-mean(rand.r2/cca.emp$tot.chi)))*(1-r2))
}}

r2.adj(spdat.n,cbind(site.m,yr.m,dist.m),reps=10)
r2.adj(spdat.n,cbind(site.m,yr.m),dist.m,reps=10)
##seems to work

r2.adj.confi<-function(Y,X,Z,reps){
 rand.r2<-rep(NA,reps)
 Y<-as.matrix(Y)
 X<-as.matrix(X)
 Z<-as.matrix(Z)
 if(missing(Z)){
  for(i in 1:reps){
   Xrand<-X[sample(nrow(X)),]
   rand.r2[i]<-summary( cca(Y,Xrand))$constr.chi
   print(i)
  }
  cca.emp<-cca(Y,X)
  r2<-summary(cca.emp)$constr.chi/cca.emp$tot.chi 
  c(r2,1-(1/(1-mean(rand.r2/cca.emp$tot.chi)))*(1-r2),1-(1/(1-rand.r2/cca.emp$tot.chi))*(1-r2))
 }
 else{
  for(i in 1:reps){
   rhold<-sample(nrow(X))
   Xrand<-X[rhold,]
   Zrand<-Z[rhold,]
   rand.r2[i]<-summary( cca(Y,Xrand,Zrand))$constr.chi
   print(i)
  }
  cca.emp<-cca(Y,X,Z)
  r2<-summary(cca.emp)$constr.chi/cca.emp$tot.chi 
  c(r2,1-(1/(1-mean(rand.r2/cca.emp$tot.chi)))*(1-r2),1-(1/(1-rand.r2/cca.emp$tot.chi))*(1-r2))
}}

test<-r2.adj.confi(spdat.n,cbind(site.m,yr.m,dist.m),reps=10)
plot(test[-(1:2)])
hist(test[-(1:2)])
test<-r2.adj.confi(spdat.n,dist.m,cbind(site.m,yr.m),100)
plot(sort(test[-(1:2)]))
hist(test[-(1:2)])
hist(test[-(1:2)],xlim=c(-.1,0.1))
points(test[1:2],c(0,0))

##
cca.full<-cca(spdat.n~name+yr.f+YrsOB+YrsSLB+BP5Yrs)
full.r2<-r2.adj.confi(spdat.n,cbind(site.m,yr.m,dist.m),reps=499)

cca.pl<-cca(spdat.n,site.m,cbind(yr.m,dist.m))
pl.r2<-r2.adj.confi(spdat.n,site.m,cbind(yr.m,dist.m),499)

cca.yr<-cca(spdat.n,yr.m,cbind(site.m,dist.m))
yr.r2<-r2.adj.confi(spdat.n,yr.m,cbind(site.m,dist.m),499)

cca.ds<-cca(spdat.n,dist.m,cbind(site.m,yr.m))
ds.r2<-r2.adj.confi(spdat.n,dist.m,cbind(site.m,yr.m),499)

cca.plyr<-cca(spdat.n,cbind(site.m,yr.m),dist.m)
plyr.r2<-r2.adj.confi(spdat.n,cbind(site.m,yr.m),dist.m,499)

cca.plds<-cca(spdat.n,cbind(site.m,dist.m),yr.m)
plds.r2<-r2.adj.confi(spdat.n,cbind(site.m,dist.m),yr.m,499)

cca.yrds<-cca(spdat.n,cbind(yr.m,dist.m),site.m)
yrds.r2<-r2.adj.confi(spdat.n,cbind(yr.m,dist.m),site.m,499)

par(mfrow=c(2,3))
hist(pl.r2[-(1:2)])
hist(yr.r2[-(1:2)])
hist(ds.r2[-(1:2)])
hist(plyr.r2[-(1:2)])
hist(plds.r2[-(1:2)])
hist(yrds.r2[-(1:2)])

##fractions
showvarparts(3)
#plot & yr
plyr.r2[1]-pl.r2[1]-yr.r2[1]
plyr.r2[2]-pl.r2[2]-yr.r2[2]
#plot & ds
plds.r2[1]-pl.r2[1]-ds.r2[1]
plds.r2[2]-pl.r2[2]-ds.r2[2]
#yr & ds
yrds.r2[1]-yr.r2[1]-ds.r2[1]
yrds.r2[2]-yr.r2[2]-ds.r2[2]
#pl&yr&ds
full.r2[1]-pl.r2[1]-yr.r2[1]-ds.r2[1]-(plyr.r2[1]-pl.r2[1]-yr.r2[1])-(plds.r2[1]-pl.r2[1]-ds.r2[1])-(yrds.r2[1]-yr.r2[1]-ds.r2[1])
full.r2[2]-pl.r2[2]-yr.r2[2]-ds.r2[2]-(plyr.r2[2]-pl.r2[2]-yr.r2[2])-(plds.r2[2]-pl.r2[2]-ds.r2[2])-(yrds.r2[2]-yr.r2[2]-ds.r2[2])

R2s<-cbind(c(full.r2[1],pl.r2[1],yr.r2[1],ds.r2[1],plyr.r2[1]-pl.r2[1]-yr.r2[1],plds.r2[1]-pl.r2[1]-ds.r2[1],yrds.r2[1]-yr.r2[1]-ds.r2[1],full.r2[1]-pl.r2[1]-yr.r2[1]-ds.r2[1]-(plyr.r2[1]-pl.r2[1]-yr.r2[1])-(plds.r2[1]-pl.r2[1]-ds.r2[1])-(yrds.r2[1]-yr.r2[1]-ds.r2[1])),
c(full.r2[2],pl.r2[2],yr.r2[2],ds.r2[2],plyr.r2[2]-pl.r2[2]-yr.r2[2],plds.r2[2]-pl.r2[2]-ds.r2[2],yrds.r2[2]-yr.r2[2]-ds.r2[2],full.r2[2]-pl.r2[2]-yr.r2[2]-ds.r2[2]-(plyr.r2[2]-pl.r2[2]-yr.r2[2])-(plds.r2[2]-pl.r2[2]-ds.r2[2])-(yrds.r2[2]-yr.r2[2]-ds.r2[2])))

colnames(R2s)<-c('R2','R2adj')
rownames(R2s)<-c('all','pl','yr','ds','plyr','plds','yrds','plyrds')
round(R2s,2)
          R2 R2adj
all     0.63  0.57
pl      0.50  0.46
yr      0.04  0.00
ds      0.01  0.00
plyr    0.01  0.03
plds    0.07  0.08
yrds    0.00  0.00
plyrds -0.01  0.00

rda(rich,site.m,cbind(yr.m,dist.m))
varpart(rich,site.m,yr.m,dist.m)

##AIC in CCA
extractAIC(cca.full)
extractAIC(cca(spdat.n~name+yr.f+YrsOB+BP5Yrs))
extractAIC(cca(spdat.n~name+yr.f+YrsOB))
extractAIC(cca(spdat.n~name+yr.f))
##CCA full has lowest AIC

###################################################
###################################################
anova(cca.full)
anova(cca.full,by="margin",step=200);alarm()
Permutation test for cca under reduced model
Marginal effects of terms
Model: cca(formula = spdat.n ~ name + yr.f + YrsOB + YrsSLB + BP5Yrs)
           Df   Chisq       F N.Perm Pr(>F)   
name      19  0.7441 13.4723    199  0.005 **
yr.f      10  0.0653  2.2465    199  0.005 **
YrsOB      1  0.0067  2.3164    199  0.005 **
YrsSLB     1  0.0061  2.1079    199  0.005 **
BP5Yrs     1  0.0047  1.6154    199  0.005 **
Residual 187  0.5436 
cca.dist<-cca(spdat.n~YrsOB+YrsSLB+Condition(name)+Condition(yr.f))
anova(cca.dist,by="m");alarm()
plot(cca.dist)
##same result as above for mangement vars

anova(cca.full,by="margin",step=200,strata=name);alarm()
Permutation test for cca under reduced model
Permutations stratified within `name'
Marginal effects of terms
Model: cca(formula = (spdat.n) ~ name + yr.f + YrsOB + bison + YrsSLB + BP5Yrs)
          Df   Chisq       F N.Perm Pr(>F)   
name      19  0.7312 13.3002    199  0.005 **
yr.f      10  0.0650  2.2448    199  0.005 **
YrsOB      1  0.0067  2.3276    199  0.005 **
bison      1  0.0054  1.8640    199  0.005 **
YrsSLB     1  0.0062  2.1299    199  0.005 **
BP5Yrs     1  0.0044  1.5242    199  0.005 **
Residual 186  0.5382                         

anova(cca.full,by="margin",step=199,strata=yr.f);alarm()
Permutation test for cca under reduced model
Permutations stratified within `yr.f'
Marginal effects of terms
Model: cca(formula = (spdat.n) ~ name + yr.f + YrsOB + bison + YrsSLB + BP5Yrs)
          Df   Chisq       F N.Perm Pr(>F)   
name      19  0.7312 13.3002    199  0.005 **
yr.f      10  0.0650  2.2448    199  0.020 * 
YrsOB      1  0.0067  2.3276    199  0.005 **
bison      1  0.0054  1.8640    199  0.005 **
YrsSLB     1  0.0062  2.1299    199  0.005 **
BP5Yrs     1  0.0044  1.5242    199  0.005 **
Residual 186  0.5382                         

########################################################
###trade out site effects for soil topo effects#########
p.int<-function(object){
 summary(object)$constr.chi /object$tot.chi
}
cca.site<-cca(spdat.n~OM+slope+northness+eastness+logCa+logFe+pH+logMg+logK+logNa+logB+logMn+logCu+logZn+logAl+
          Condition(YrsOB)+Condition(YrsSLB)+Condition(BP5Yrs)+Condition(yr.f))
hold<-step(cca.site)

plot(cca.site,type='n')
text(cca.site, display = "cn", col = "blue",cex=.75) 
plot(cca.site)
anova(cca.site,by="margin",perm=199);alarm()

##subset of those variables 
ca.a<-rep(tapply(logCa,name,mean),each=11)

cca.site<-cca(spdat.n~logCa+slope+northness+Condition(YrsOB)+Condition(YrsSLB)+
                        Condition(BP5Yrs)+Condition(yr.f))
p.int(cca.site)
[1] 0.1299383
r2.adj(spdat.n,cbind(logCa,slope,northness),cbind(dist.m,yr.m),499)
[1] 0.1299383 0.1179565  ##w/100 reps
[1] 0.1299383 0.1179955  ##w/499 reps
cca.ca<-cca(spdat.n,logCa,cbind(slope,northness,dist.m,yr.m))
p.int(cca.ca)
r2.adj(spdat.n,logCa,cbind(slope,northness,dist.m,yr.m),reps=100)
[1] 0.07224288 0.06795699

cca.slope<-cca(spdat.n,slope,cbind(logCa,northness,dist.m,yr.m))
p.int(cca.slope)
r2.adj(spdat.n,slope,cbind(logCa,northness,dist.m,yr.m),reps=100)
[1] 0.02983976 0.02553035

cca.north<-cca(spdat.n,northness,cbind(logCa,slope,dist.m,yr.m))
p.int(cca.north)
r2.adj(spdat.n,northness,cbind(logCa,slope,dist.m,yr.m),reps=100)
[1] 0.01910410 0.01466706

cca.bison<-cca(spdat.n,YrsOB,cbind(BP5Yrs,YrsSLB,site.m,yr.m))
p.int(cca.bison)
[1] 0.004564033
r2.adj(spdat.n,YrsOB,cbind(BP5Yrs,YrsSLB,site.m,yr.m),reps=100)
[1]  4.564033e-03 -1.302320e-05
##0,0

cca.bp5<-cca(spdat.n,BP5Yrs,cbind(YrsOB,YrsSLB,site.m,yr.m))
p.int(cca.bp5)
[1] 0.003182903

cca.yrslb<-cca(spdat.n,YrsSLB,cbind(YrsOB,BP5Yrs,site.m,yr.m))
p.int(cca.yrslb)
[1] 0.004153339


pl<-plot(cca.site,type='n')
text(cca.site, display = "cn", col = "blue",cex=.75) 
points(cca.site, display = "sp",cex=.5)
identify(pl, "sp", labels = colnames(spdat.n),cex=.75)

anova(cca.site,by="margin",step=1000,strata=yr);alarm()
Permutation test for cca under reduced model
Permutations stratified within `yr'
Marginal effects of terms

Model: cca(formula = spdat.n ~ logCa + slope + northness + Condition(YrsOB) + Condition(YrsSLB) + Condition(BP5Yrs) + Condition(yr.f))
           Df   Chisq       F N.Perm Pr(>F)   
logCa       1  0.1066 19.7412    999  0.001 **
slope       1  0.0440  8.1541    999  0.001 **
northness   1  0.0282  5.2204    999  0.001 **
Residual  203  1.0960                      

cca.yr<-cca(spdat.n~PDSIavg+spr.rain.tot+sum.rain.tot+win.rain.tot+spr.temp+sum.temp+win.temp+Condition(YrsOB)+Condition(YrsSLB)+Condition(BP5Yrs)+Condition(name))
p.int(cca.yr)
[1] 0.03285512
r2.adj(spdat.n,cbind(PDSIavg,spr.rain.tot,sum.rain.tot,win.rain.tot,spr.temp,sum.temp,win.temp),cbind(dist.m,site.m),100)
[1] 0.032855125 0.001458442

test<-cca(spdat.n~PDSIavg+Condition(YrsOB)+Condition(YrsSLB)+Condition(BP5Yrs)+Condition(name))
r2.adj(test


plot(cca.yr)
anova(cca.yr,by='m')
Permutation test for cca under reduced model
Marginal effects of terms

Model: cca(formula = spdat.n ~ PDSIavg + spr.rain.tot + sum.rain.tot + win.rain.tot + spr.temp + sum.temp + win.temp + Condition(YrsOB) + Condition(YrsSLB) + Condition(BP5Yrs) + Condition(name))
              Df  Chisq      F N.Perm Pr(>F)   
PDSIavg        1 0.0080 2.7240    199  0.005 **
spr.rain.tot   1 0.0072 2.4263    199  0.005 **
sum.rain.tot   1 0.0061 2.0750    199  0.005 **
win.rain.tot   1 0.0064 2.1549    199  0.005 **
spr.temp       1 0.0050 1.6789    199  0.005 **
sum.temp       1 0.0082 2.7723    199  0.005 **
win.temp       1 0.0051 1.7239    199  0.005 **
Residual     190 0.5604                       
anova(cca.yr,by='m',strata=name)
##same result as above
cca.dist<-cca(spdat.n,dist.m,cbind(site.m,yr.m))
anova(cca.dist)
Model: cca(X = spdat.n, Y = dist.m, Z = cbind(site.m, yr.m))
          Df  Chisq      F N.Perm Pr(>F)   
Model      3 0.0187 2.1473    199  0.005 **
Residual 187 0.5436                        
cca.site<-cca(spdat.n,site.m,cbind(dist.m,yr.m))
anova(cca.site)
Model: cca(X = spdat.n, Y = site.m, Z = cbind(dist.m, yr.m))
          Df   Chisq       F N.Perm Pr(>F)   
Model     19  0.7441 13.4723    199  0.005 **
Residual 187  0.5436               
cca.yr<-cca(spdat.n,yr.m,cbind(dist.m,site.m))
anova(cca.yr)
Model: cca(X = spdat.n, Y = yr.m, Z = cbind(dist.m, site.m))
          Df  Chisq      F N.Perm Pr(>F)   
Model     10 0.0653 2.2465    199  0.005 **
Residual 187 0.5436   


############
##more sophisticated permuations
##torodial rotations
##this is the reduced model method
cca.perm<-function(y, x, z, strata = NA, torus = F,num.reps = 100) 
{
 y<-as.matrix(y)
 x<-as.matrix(x)
 z<-as.matrix(z)
 nsamples <- nrow(x)
 cca.emp <- cca(y,cbind(x,z))
 res.emp <- summary(cca.emp)$unconst.chi
 con.emp <- summary(cca.emp)$constr.chi
 F.val <- rep(NA,num.reps)
 F.val[1] <- (con.emp/(ncol(x)+ncol(z)))/(res.emp/(nsamples-ncol(x)-ncol(z)))
 cca.cov <- cca(y,z)
 y.fit<- fitted.values(cca.cov)
 y.res<- residuals(cca.cov)
 if(!is.na(strata[1])){ 
  uni.strata <- unique(strata) 
  nstrata <- length(uni.strata)
 }
 for (i in 2:num.reps) {
  if (is.na(strata[1])){
   o <- sample(nsamples)
   y.res <- y.res[o,]
  }
  else{
   for (j in 1:nstrata){ 
    srows <- (1:nsamples)[strata==uni.strata[j]]
    if(torus==T){#rows must be properly ordered already in dataset
     hold <- rep(srows,2)
     rstart <- sample(nsamples/nstrata,1)
     o <- hold[rstart:(rstart+((nsamples/nstrata)-1))]
     y.res[srows,] <- y.res[o,]
    }
    else{ 
     o <- sample(srows)
     y.res[ srows,] <- y.res[o,]
  }}}
  y.new <- y.fit + y.res
  cca.perm <- cca(y.new,cbind(x,z))
  res.perm <- summary(cca.perm)$unconst.chi
  con.perm <- summary(cca.perm)$constr.chi
  F.val[i]<- (con.perm/(ncol(x)+ncol(z)))/(res.perm/(nsamples-ncol(x)-ncol(z)))
  print(i)
 }
F.val
}

ca.tor<-cca.perm(spdat.n,ca.a,cbind(slope,northness,dist.m,yr.m),strata=yr,torus=T,num.reps=200);alarm()

ca.site<-cca.perm(spdat.n,ca.a,cbind(slope,northness,dist.m,yr.m),strata=name,num.reps=10);alarm()

ca.rand<-cca.perm(spdat.n,ca.a,cbind(slope,northness,dist.m,yr.m),num.reps=200);alarm()

sd(ca.tor[-1])
sd(ca.rand[-1])

hist(ca.tor[-1],breaks=10,col='red')
par(new=T)
hist(ca.rand[-1],breaks=10)

}
cca.perm(spdat.n,cbind(ca.a,slope,northness),cbind(dist.m,yr.m))

test<-cca(spdat.n,cbind(ca.a,slope,northness),cbind(dist.m,yr.m))
test
anova(test)$F

   res.emp <- summary(test)$unconst.chi
    con.emp <- summary(test)$constr.chi
    F.emp <- (con.emp/1)/(res.emp/(220-1-10-3-2-1))
nrow(as.matrix(ca.a))




pairs(~logCa+logFe+pH+OM)

site.pc<-prcomp(cbind(logFe,pH,logZn,northness,eastness,slope),scale=T)
plot(site.pc)
biplot(site.pc,cex=.5)
summary(site.pc)
site.pc1<-site.pc$x[,1]



site.m<-as.matrix(endatlv1[,50:66])
#rownames(site.m)<-name
princomp(site.m)
site.pc<-prcomp(site.m,scale=T)
plot(site.pc)
biplot(site.pc,cex=.5)
summary(site.pc)
site.pc1<-site.pc$x[,1]

cca.site.pc<-cca((spdat.n)~site.pc1+YrsOB+bison+YrsSLB+BP5Yrs)
cca.site.pc
anova(cca.site.pc,by="margin",perm=199,strata=yr.f);alarm()
site.pc1   1  0.1151 20.0104    199  0.005 **
YrsOB      1  0.0182  3.1600    199  0.005 **
bison      1  0.0185  3.2191    199  0.005 **
YrsSLB     1  0.0131  2.2847    199  0.005 **
BP5Yrs     1  0.0352  6.1194    199  0.005 **
Residual 214  1.2307
test<-step(cca((spdat.n)~logCa+logFe+pH+bison+YrsOB+YrsSLB+BP5Yrs+rain.tot+temp.avg))


cca.time<-cca((spdat.n)~spr.rain.tot+win.rain.tot+sum.rain.tot+spr.temp+win.temp+sum.temp+YrsOB+bison+YrsSLB+BP5Yrs+Condition(name))
pl<-plot(cca.time,display="sp",type='n')
totcov<-colSums(spdat.n)
##trys to make the sp names pretty
sel <- orditorp(cca.time, dis = "sp", lab = colnames(spdat.n), pcol = "gray", pch = "+")
sel <- orditorp(cca.time, dis = "sp", lab = colnames(spdat.n),priority = totcov, pcol = "gray", pch = "+")
identify(pl, "sp", labels = colnames(spdat.n),cex=.75)

anova(cca.time,by="margin");alarm()

###model selction###
soil.m<-as.data.frame(cbind(pH,endatlv1[,56:66]))
##
cca.base<-cca(spdat.n~name+yr.f)
resids<-residuals(cca.base)
cca.soils<-cca(resids,soil.m)
cca.soils.red<-cca(resids~logCa+logFe+pH)
cca.dist<-cca(resids,dist.m)
AIC.cca<-rbind(extractAIC(cca.base),extractAIC(cca.soils),extractAIC(cca.soils.red),extractAIC(cca.dist))
colnames(AIC.cca)<-c('rank','AIC')
rownames(AIC.cca)<-c('base','soils','soil red','dist')
AIC.cca
#######
soil.m<-as.matrix(cbind(pH,endatlv1[,56:66]))
soil.m.red<-as.matrix(cbind(pH,logCa,logFe))
topo.m<-as.matrix(endatlv1[,50:52])
clim.m<-as.matrix(cbind(spr.rain.tot,win.rain.tot,sum.rain.tot,spr.temp,win.temp,sum.temp))
pclim.m<-as.matrix(cbind(pspr.rain.tot,pwin.rain.tot,psum.rain.tot,pspr.temp,pwin.temp,psum.temp))
clim.m<-as.matrix(cbind(endatlv1[,96:98],endatlv1[,102:104]))
pclim.m<-as.matrix(cbind(pspr.rain.tot,pwin.rain.tot,psum.rain.tot,pspr.temp,pwin.temp,psum.temp))

cca1<-cca(spdat.n~dist.m+topo.m+clim.m+soil.m.red)
cca1
anova(cca1,by='m');alarm()
Permutation test for cca under reduced model
Marginal effects of terms
Model: cca(formula = spdat.n ~ dist.m + topo.m + clim.m + soil.m.red)
            Df   Chisq       F N.Perm Pr(>F)   
dist.m       4  0.0760  3.9434    199  0.005 **
topo.m       3  0.0926  6.4046    199  0.005 **
clim.m       6  0.0477  1.6503    199  0.005 **
soil.m.red   3  0.1853 12.8216    199  0.005 **
Residual   203  0.9780                    

cca.topo<-cca(spdat.n~slope+northness+eastness+Condition(dist.m)+Condition(clim.m)+Condition(soil.m.red))
cca.topo.cons<-cca(spdat.n~slope+northness+eastness+Condition(dist.m)+Condition(clim.m)+Condition(soil.m))
anova(cca.topo,by='m');alarm()
Model: cca(formula = spdat.n ~ slope + northness + eastness + Condition(dist.m) + Condition(clim.m) + Condition(soil.m.red))
           Df  Chisq      F N.Perm Pr(>F)   
slope       1 0.0293 6.0831    199  0.005 **
northness   1 0.0293 6.0753    199  0.005 **
eastness    1 0.0344 7.1337    199  0.005 **
Residual  203 0.9780                  
anova(cca.topo.cons,by='m');alarm()
slope       1 0.0243 5.8040    199  0.005 **
northness   1 0.0214 5.1073    199  0.005 **
eastness    1 0.0238 5.6947    199  0.005 **

pl<-plot(cca.topo,display=c('sp','bp'))
identify(pl,'sp',labels=colnames(spdat.n),cex=.75)
##no disernable patterns
cca.soil.red<-cca(spdat.n~soil.m.red+Condition(dist.m)+Condition(clim.m)+Condition(topo.m))
pl<-plot(cca.soil.red,display=c('sp','bp'))
identify(pl,'sp',labels=colnames(spdat.n),cex=.75)


####soil PCA
clim.m<-cbind(sum.rain.tot,win.rain.tot,spr.rain.tot,sum.temp,win.temp,spr.temp)
#clim.m<-cbind(endatlv1[,67:68],endatlv1[,72:95])
rownames(clim.m)<-name
princomp(clim.m)
clim.pc<-prcomp(clim.m,scale=T)
plot(clim.pc)
biplot(clim.pc,cex=.5)
summary(clim.pc)
clim.pc1<-clim.pc$x[,1]

site.m<-as.matrix(endatlv1[,53:66])
rownames(site.m)<-name
princomp(site.m)
site.pc<-prcomp(site.m,scale=T)
plot(site.pc)
biplot(site.pc,cex=.5)
summary(site.pc)
site.pc1<-site.pc$x[,1]


##############################
##############################
##Dissimilarity Analyses######
yr.m<-endatlv1[,28:38]
pl.m<-endatlv1[,8:27]

spdist<-vegdist(spdat.n,meth='jacc')
yrdist<-vegdist(yr,method='euclid')
pldist<-vegdist(pl.m,method='euclid')

par(mfrow=c(1,3))
plot(yrdist,spdist,type='n')
lines(lowess(yrdist,spdist),lwd=2,col='red')
plot(pldist,spdist,type='n')
lines(lowess(pldist,spdist),lwd=2,col='red')
plot(pldist,yrdist,type='n')
lines(lowess(pldist,yrdist),lwd=2,col='red')


mantel(spdist,yrdist);alarm()
Mantel statistic r: 0.08181 
      Significance: < 0.001

mantel.partial(spdist,yrdist,pldist);alarm()
Mantel statistic r: 0.1191 
      Significance: < 0.001 

mantel(spdist,pldist);alarm()
Mantel statistic r: 0.5644 
      Significance: < 0.001

mantel.partial(spdist,pldist,yrdist);alarm()
Mantel statistic r: 0.569 
      Significance: < 0.001 

mantel.partial(spdist,pldist,yrdist,strata=name);alarm()
Mantel statistic r: 0.569 
      Significance: < 0.001

test<-mantel.partial(spdist,pldist,yrdist,strata=c(yr,2009));alarm()
Mantel statistic r: 0.569 
      Significance: < 0.001


test<-matrix(NA,ncol=3,nrow=3)
test[,1]<-c(1,1,0)
test[,2]<-c(0,1,1)
test[,3]<-c(0,0,1)
print(test<-vegdist(test,meth='jacc'))
vegdist(test,meth='jacc')

mantel.cor(spdist,yrdist,breaks=0:11)
Error in if (mantel.result$statistic >= 0) { : 
  missing value where TRUE/FALSE needed
In addition: There were 50 or more warnings (use warnings() to see the first 50)
> warnings()
Warning messages:
1: In cor(x, y) ... : the standard deviation is zero
2: In cor(permvec, ydis, method = method) ... : the standard deviation is zero

##################################
##################################
library(ecodist)

envi.d<-distance(cbind(soil.m,topo.m,dist.m),"eucl")
time.d<-distance(yr,"eucl")
veg.d<-distance(spdat.n,"jacc")
v.yr.mgram<- mgram(veg.d,time.d,breaks=0:10)
v.yr.mgram
plot(v.yr.mgram,p=0.05/10)

v.en.mgram<- mgram(veg.d,envi.d)
v.en.mgram
plot(v.en.mgram,p=0.05/ncol(v.en.mgram$mgram))

v.yr.en<-pmgram(veg.d,time.d,envi.d,breaks=0:10)
plot(v.v.en,p=0.05/10)


