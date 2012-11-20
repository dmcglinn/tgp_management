endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
names(endat)
attach(endat)
library(nlme)
##work only with level 1
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)
yr.f<-factor(yr)

###Variation paritioning
##Palmer/McGlinn approach
##this is the way that makes the most sense to me but gives really large 
##negative shared varaiation
dist.m<-cbind(YrsOB,BP5Yrs,YrsSLB,bison)

##########################################
summary(lm(rich~-1+name+yr.f+dist.m))
mod<-gls(rich~-1+name+yr.f+dist.m,meth="ML")
anova(mod,type='m')##marginal F-tests
Denom. DF: 186 
       numDF  F-value p-value
name      20 57.81748  <.0001
yr.f      10  9.06527  <.0001
dist.m     4  6.96128  <.0001
plot(yr,resid(mod,resType="p"))
lines(lowess(yr,resid(mod,resType="n")),lwd=2)
plot(Variogram(mod,form=~yr|name))
plot(ACF(mod,form=~yr|name),alpha=0.05)
##need to consider autocorrelation
mod2<-update(mod,corr=corExp(form=~yr|name))
anova(mod,mod2)
mod2n<-update(mod,corr=corExp(form=~yr|name,nugget=T))
anova(mod2,mod2n)
plot(Variogram(mod2,form=~yr|name,robust=F),ylim=c(0,2))
plot(Variogram(mod2,form=~yr|name,robust=F,resType="n"),ylim=c(-2,2))
mod3<-update(mod,corr=corSpher(form=~yr|name,nugget=FALSE))
mod3n<-update(mod,corr=corSpher(form=~yr|name,nugget=TRUE))
anova(mod3,mod3n)
##nugget uncess with spher model
anova(mod2,mod3)
##expon much better than spher
mod4<-update(mod,corr=corRatio(form=~yr|name,nugget=FALSE))
mod4n<-update(mod,corr=corRatio(form=~yr|name,nugget=TRUE))
anova(mod4,mod4n)
anova(mod4n,mod2)
mod5<-update(mod,corr=corLin(form=~yr|name,nugget=FALSE))
mod5n<-update(mod,corr=corLin(form=~yr|name,nugget=TRUE))
anova(mod5,mod5n)
anova(mod5,mod2)
mod5<-update(mod,corr=corGaus(form=~yr|name,nugget=FALSE))
mod5n<-update(mod,corr=corGaus(form=~yr|name,nugget=TRUE))
anova(mod5,mod5n)
anova(mod5,mod2)
plot(Variogram(mod5n,form=~yr|name,robust=T),ylim=c(0,1.5))
plot(Variogram(mod5n,form=~yr|name,robust=F,resType="n"),ylim=c(-1.2,1.2))
##so expon still best
mod6<-update(mod,corr=corARMA(p=1,form=~yr|name))
mod6b<-update(mod,corr=corARMA(p=2,form=~yr|name))
mod6c<-update(mod,corr=corARMA(p=3,form=~yr|name))
anova(mod6,mod6b,mod6c)
anova(mod,mod6,mod2)
##exponential and AR1 are identical models ??
plot(ACF(mod6,resType="n",maxLag=5),alpha=0.05)
anova(mod2,type='m')##adjusted for autocorrelation
Denom. DF: 186 
       numDF  F-value p-value
name      20 48.67666  <.0001
yr.f      10 10.17113  <.0001
dist.m     4  5.36843   4e-04

plot(mod,main="not corrected")
plot(mod6)##standardized
plot(mod6,resid(.,type='n')~fitted(.))##normalized

mod7<-update(mod6,weights=varPower())
anova(mod6,mod7)
##heteroscadistic term not needed
intervals(mod6)
summary(mod6)

##ACF plots
library(lattice)
trellis.device(color=FALSE)
plot(ACF(mod6,resType='r',maxLag=5), alpha = 0.05)

ACF(mod6,resType='r',maxLag=5)
plot(ACF(mod6,resType='r',maxLag=5),alpha=0.05)
plot(ACF(mod6,resType='n',maxLag=5),alpha=0.05)

##begin var part with R2 calculations
##based on reg.ss/tot.ss or 1-err.ss/tot.ss
##experimentation shows that ssreg+sserr != sstot and therefore this 
##method should not be used
sstot<-sum((rich-mean(rich))^2)
sserr<-sum((rich-mod6$fit)^2)
1-sserr/sstot
[1] 0.7788224
##alternatively
ssreg<-sum((mod6$fit-mean(mod6$fit))^2)
ssreg/sstot
[1] 0.7745776 ##slightly smaller than above? (rounding errors?)
ssreg+sserr-sstot
[1] -152.3601
##indicates that this is not an appropriate technique
mod.full<-gls(rich~name+yr.f+dist.m,cor=corAR1(form=~yr|name),method="ML")
resids<-resid(update(mod.full,.~.-name))
mod.pl<-gls(resids~name,cor=corAR1(form=~yr|name),meth="ML")
resids<-resid(update(mod.full,.~.-yr.f))
mod.yr<-gls(resids~yr.f,cor=corAR1(form=~yr|name),meth="ML")
resids<-resid(update(mod.full,.~.-dist.m))
mod.ds<-gls(resids~dist.m,cor=corAR1(form=~yr|name),meth="ML")
resids<-resid(update(mod.full,.~.-name-yr.f))
mod.plyr<-gls(resids~name+yr.f,cor=corAR1(form=~yr|name),meth="ML")
resids<-resid(update(mod.full,.~.-name-dist.m))
mod.plds<-gls(resids~name+dist.m,cor=corAR1(form=~yr|name),meth="ML")
resids<-resid(update(mod.full,.~.-yr.f-dist.m))
mod.yrds<-gls(resids~yr.f+dist.m,cor=corAR1(form=~yr|name),meth="ML")

##explainable var
tot.ss<-sum((rich-mean(rich))^2)
##35893.71 is Sum Sq total (total variance in dataset)
full.ss<-sum((mod.full$fit-mean(mod.full$fit))^2)
pl.ss<-sum((mod.pl$fit-mean(mod.pl$fit))^2)
yr.ss<-sum((mod.yr$fit-mean(mod.yr$fit))^2)
ds.ss<-sum((mod.ds$fit-mean(mod.ds$fit))^2)
plyr.ss<-sum((mod.plyr$fit-mean(mod.plyr$fit))^2)-pl.ss-yr.ss
plds.ss<-sum((mod.plds$fit-mean(mod.plds$fit))^2)-pl.ss-ds.ss
yrds.ss<-sum((mod.yrds$fit-mean(mod.yrds$fit))^2)-ds.ss-yr.ss
plyrds.ss<-full.ss-pl.ss-yr.ss-ds.ss-plyr.ss-plds.ss-yrds.ss
c(full.ss/tot.ss,pl.ss/full.ss,yr.ss/full.ss,ds.ss/full.ss,plyr.ss/full.ss,plds.ss/full.ss,yrds.ss/full.ss,plyrds.ss/full.ss)
r2<-c(full.ss/tot.ss,pl.ss/tot.ss,yr.ss/tot.ss,ds.ss/tot.ss,plyr.ss/tot.ss,plds.ss/tot.ss,yrds.ss/tot.ss,plyrds.ss/tot.ss)
##function of r2 adj, Peres-Neto et al. 2006
r2adj<-function(reg.ss,p){1-((220-1)/(220-p-1)*(1-reg.ss/tot.ss))}
ss<-c(full.ss,pl.ss,yr.ss,ds.ss,plyr.ss,plds.ss,yrds.ss,plyrds.ss)
np<-c(19+10+4,19,10,4,19+10,19+4,10+4,0)

##function of r2 adj, Peres-Neto et al. 2006
r2adj<-function(reg.ss,p){1-((220-1)/(220-p-1)*(1-reg.ss/tot.ss))}
ss<-c(full.ss,pl.ss,yr.ss,ds.ss,plyr.ss,plds.ss,yrds.ss,plyrds.ss)
np<-c(19+10+4,19,10,4,19+10,19+4,10+4,0)
##calc r2adj
r2a<-mapply(r2adj,ss,np)
varexp<-cbind(r2,r2a)
rownames(varexp)<-c('expl var','plot','yr','dist','plot&yr','plot&dist','yr&dist','plot&yr&dist')
colnames(varexp)<-c('R2','R2adj')
round(varexp,3)
                 R2  R2adj
expl var      0.775  0.735
plot          0.546  0.503
yr            0.131  0.089
dist          0.024  0.006
plot&yr      -0.023 -0.179
plot&dist     0.011 -0.105
yr&dist       0.050 -0.015
plot&yr&dist  0.035  0.035


##begin var part with pseudo R2 calculations
mod.emp<-gls(rich~1,cor=corAR1(form=~yr|name),meth="ML")
mod.full<-update(mod.emp,.~.-1+name+yr.f+dist.m)
pR2.gls(mod.full)$r2ML
[1] 0.6532633
resids<-resid(update(mod.full,.~.-name),type='n')
mod.pl<-gls(resids~name,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.pl)$r2ML
[1] 0.286606
resids<-resid(update(mod.full,.~.-yr.f),type='n')
mod.yr<-gls(resids~yr.f,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.yr)$r2ML
[1] 0.3830418
resids<-resid(update(mod.full,.~.-dist.m),type='n')
mod.ds<-gls(resids~dist.m,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.ds)$r2ML
[1] 0.08355011
resids<-resid(update(mod.full,.~.-name-yr.f),type='n')
mod.plyr<-gls(resids~name+yr.f,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.plyr)$r2ML
[1] 0.5144256
resids<-resid(update(mod.full,.~.-name-dist.m),type='n')
mod.plds<-gls(resids~name+dist.m,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.plds)$r2ML
[1] 0.3486134
resids<-resid(update(mod.full,.~.-yr.f-dist.m),type='n')
mod.yrds<-gls(resids~yr.f+dist.m,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.yrds)$r2ML
[1] 0.4744509

pR2s<-rbind(pR2.gls(mod.full)$r2ML,
pR2.gls(mod.pl)$r2ML,
pR2.gls(mod.yr)$r2ML,
pR2.gls(mod.ds)$r2ML,
pR2.gls(mod.plyr)$r2ML-pR2.gls(mod.pl)$r2ML-pR2.gls(mod.yr)$r2ML,
pR2.gls(mod.plds)$r2ML-pR2.gls(mod.pl)$r2ML-pR2.gls(mod.ds)$r2ML,
pR2.gls(mod.yrds)$r2ML-pR2.gls(mod.yr)$r2ML-pR2.gls(mod.ds)$r2ML,
pR2.gls(mod.full)$r2ML-pR2.gls(mod.pl)$r2ML-pR2.gls(mod.yr)$r2ML-pR2.gls(mod.ds)$r2ML)
rownames(pR2s)<-c('full','plot','yr','dist','plotyr','plotdist','yrdist','plotyrdist')
colnames(pR2s)<-"pR2"
pR2s
                   pR2
full        0.64512447
plot        0.31857461
yr          0.37433136
dist        0.05375850
plotyr     -0.07145874
plotdist    0.01743712
yrdist      0.02572988
plotyrdist -0.10154001

##pseudo R2 formulas
##http://en.wikipedia.org/wiki/Coefficient_of_determination#Generalized_R.C2.A02
##N. Nagelkerke, “A Note on a General Definition of the Coefficient of Determination,” Biometrika, vol. 78, no. 3, pp. 691-692, 1991.
pR2Work <- function(llh,llhNull,n){
     McFadden <- 1 - llh/llhNull
     G2 <- -2*(llhNull-llh)
     r2ML <- 1 - exp(-G2/n)
     r2ML.max <- 1 - exp(llhNull*2/n)
     r2CU <- r2ML/r2ML.max
     out <- c(llh=llh,
              llhNull=llhNull,
              G2=G2,
              McFadden=McFadden,
              r2ML=r2ML,
              r2CU=r2CU)
     as.list(out)
}

pR2.glm <- function(object,...){
     llh <- logLik(object)
     objectNull <- update(object, ~ 1)
     llhNull <- logLik(objectNull)
     n <- dim(object$model)[1]
     pR2Work(llh,llhNull,n)
}

##dan's mod
pR2.gls <- function(object,...){
     llh <- logLik(object)
     objectNull <- update(object, ~ 1)
     llhNull <- logLik(objectNull)
     n <- as.numeric(object$dims[1])
     pR2Work(llh,llhNull,n)
}

###################################################3
###################################################
###Without Bison Binary Variable##################

###Variation paritioning
##Palmer/McGlinn approach
##this is the way that makes the most sense to me but gives really large 
##negative shared varaiation
dist.m<-cbind(YrsOB,BP5Yrs,YrsSLB)

##########################################
summary(lm(rich~-1+name+yr.f+dist.m))
mod<-gls(rich~-1+name+yr.f+dist.m,meth="ML")
anova(mod,type='m')##marginal F-tests
Denom. DF: 187 
       numDF  F-value p-value
name      20 65.43906  <.0001
yr.f      10  9.16105  <.0001
dist.m     3  9.37591  <.0001
plot(yr,resid(mod,resType="p"))
lines(lowess(yr,resid(mod,resType="n")),lwd=2)
plot(ACF(mod,form=~yr|name),alpha=0.05)
mod2<-update(mod,.~.,cor=corAR1(form=~yr|name))
mod3<-update(mod,.~.,cor=corARMA(p=2,form=~yr|name))
anova(mod,mod2,mod3)
plot(ACF(mod2,form=~yr|name,maxLag=5,resType='n'),alpha=0.05)
plot(mod2)
mod4<-update(mod2,weights=varPower())
anova(mod2,mod4)
##hetero term not necessary
###############################
##begin var part with pseudo R2 calculations
mod.emp<-gls(rich~1,cor=corAR1(form=~yr|name),meth="REML")
mod.full<-update(mod.emp,.~.+name+yr.f+dist.m)
mod.pl<-update(mod.full,.~.-yr.f-dist.m)
mod.yr<-update(mod.full,.~.-name-dist.m)
mod.ds<-update(mod.full,.~.-name-yr.f)##note 3 auto corr terms are more appr here
mod.plyr<-update(mod.full,.~.-dist.m)
mod.plds<-update(mod.full,.~.-yr.f)
mod.yrds<-update(mod.full,.~.-name)

pR2diff<-function(modfull,modsub){
 pR2.gls(modfull)$r2ML-pR2.gls(modsub)$r2ML
}
full.r2<-pR2.gls(mod.full)$r2ML
pl.r2<-pR2diff(mod.full,mod.yrds)
yr.r2<-pR2diff(mod.full,mod.plds)
ds.r2<-pR2diff(mod.full,mod.plyr)
##shared components
plyr.r2<-pR2diff(mod.full,mod.ds)-pl.r2-yr.r2
plds.r2<-pR2diff(mod.full,mod.yr)-pl.r2-ds.r2
yrds.r2<-pR2diff(mod.full,mod.pl)-yr.r2-ds.r2
plyrds.r2<-full.r2-pl.r2-yr.r2-ds.r2-plyr.r2-plds.r2-yrds.r2

##corrected pR2
pR2diff<-function(modfull,modsub){
 pR2.gls(modfull)$r2CU-pR2.gls(modsub)$r2CU
}
full.r2c<-pR2.gls(mod.full)$r2CU
pl.r2c<-pR2diff(mod.full,mod.yrds)
yr.r2c<-pR2diff(mod.full,mod.plds)
ds.r2c<-pR2diff(mod.full,mod.plyr)
##shared components
plyr.r2c<-pR2diff(mod.full,mod.ds)-pl.r2-yr.r2
plds.r2c<-pR2diff(mod.full,mod.yr)-pl.r2-ds.r2
yrds.r2c<-pR2diff(mod.full,mod.pl)-yr.r2-ds.r2
plyrds.r2c<-full.r2-pl.r2-yr.r2-ds.r2-plyr.r2-plds.r2-yrds.r2

parts<-c(full.r2,pl.r2,yr.r2,ds.r2,plyr.r2,yrds.r2,plds.r2,plyrds.r2)
partsc<-c(full.r2c,pl.r2c,yr.r2c,ds.r2c,plyr.r2c,yrds.r2c,plds.r2c,plyrds.r2c)
names<-c('full','pl','yr','ds','plyr','plds','yrds','plyrds')
data.frame(names=names,parts=parts,partsc=partsc)
   names        parts       partsc
1   full  0.770107878  0.770557560
2     pl  0.269132808  0.269289960
3     yr  0.184730544  0.184838412
4     ds  0.018956096  0.018967165
5   plyr  0.252689356  0.253101926
6   plds  0.049016543  0.049164102
7   yrds  0.001218682  0.001387614
8 plyrds -0.005636151 -0.005636151
##so in this case very little difference in Cox and Nagelkerke's def
##of a generalized R2 (i.e. max R2 of ratio def is near 1. 
anova(mod.full,type='m')

pR2.gls(mod.full)$r2ML
resids<-resid(update(mod.full,.~.-name),type='n')
mod.pl<-gls(resids~name,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.pl)$r2ML
resids<-resid(update(mod.full,.~.-yr.f),type='n')
mod.yr<-gls(resids~yr.f,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.yr)$r2ML
resids<-resid(update(mod.full,.~.-dist.m),type='n')
mod.ds<-gls(resids~dist.m,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.ds)$r2ML
resids<-resid(update(mod.full,.~.-name-yr.f),type='n')
mod.plyr<-gls(resids~name+yr.f,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.plyr)$r2ML
resids<-resid(update(mod.full,.~.-name-dist.m),type='n')
mod.plds<-gls(resids~name+dist.m,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.plds)$r2ML
[1] 0.3486134
resids<-resid(update(mod.full,.~.-yr.f-dist.m),type='n')
mod.yrds<-gls(resids~yr.f+dist.m,cor=corAR1(form=~yr|name),meth="ML")
pR2.gls(mod.yrds)$r2ML
[1] 0.4744509

pR2s<-rbind(pR2.gls(mod.full)$r2ML,
pR2.gls(mod.pl)$r2ML,
pR2.gls(mod.yr)$r2ML,
pR2.gls(mod.ds)$r2ML,
pR2.gls(mod.plyr)$r2ML-pR2.gls(mod.pl)$r2ML-pR2.gls(mod.yr)$r2ML,
pR2.gls(mod.plds)$r2ML-pR2.gls(mod.pl)$r2ML-pR2.gls(mod.ds)$r2ML,
pR2.gls(mod.yrds)$r2ML-pR2.gls(mod.yr)$r2ML-pR2.gls(mod.ds)$r2ML,
pR2.gls(mod.full)$r2ML-pR2.gls(mod.pl)$r2ML-pR2.gls(mod.yr)$r2ML-pR2.gls(mod.ds)$r2ML)
rownames(pR2s)<-c('full','plot','yr','dist','plotyr','plotdist','yrdist','plotyrdist')
colnames(pR2s)<-"pR2"
pR2s
                   pR2
full        0.64508830
plot        0.33261118
yr          0.37440794
dist        0.04933874
plotyr     -0.07802680
plotdist    0.00315778
yrdist      0.02658328
plotyrdist -0.11126956




