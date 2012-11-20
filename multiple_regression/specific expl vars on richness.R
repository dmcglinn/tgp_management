
##this analysis will mimic the pCCA analysis that attempts to 
##examine specific explanatory variables for site, yr, and management
endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
library(nlme)
library(car)
endatlv1<-endat[endat$lv==1&endat$cor==1,] ##to only work at the 100 m2 scale
attach(endatlv1)
##within and between plot variance can be quantified with aov or with lme
yr.f<-factor(yr)
srain<-scale(rain.tot)
stemp<-scale(temp.avg)
psrain<-scale(prain.tot)
pstemp<-scale(ptemp.avg)
sPDSI<-scale(PDSIavg)
ssum.rain<-scale(sum.rain.tot)
sspr.rain<-scale(spr.rain.tot)
swin.rain<-scale(win.rain.tot)
ssum.temp<-scale(sum.temp)
sspr.temp<-scale(spr.temp)
swin.temp<-scale(win.temp)
syob<-scale(YrsOB,scale=T)##b/c a continuous time var
sbp5<-scale(BP5Yrs,scale=T)##
syslb<-scale(YrsSLB,scale=T)##b/c a cont var
dist.m<-cbind(YrsOB,BP5Yrs,YrsSLB)
##################################
yr.m<-endatlv1[,28:37]
site.m<-endatlv1[,8:26]

library(vegan)
mod.clim<-gls(rich~psum.rain.tot+pwin.rain.tot+pspr.rain.tot+sum.rain.tot+win.rain.tot+spr.rain.tot+name+dist.m,meth="REML",cor=corAR1(form=~yr|name))
round(summary(mod.clim)$t,3)
test1<-gls(rich~psum.rain.tot+pwin.rain.tot+pspr.rain.tot+name+dist.m,cor=corAR1(form=~yr|name))
test2<-gls(rich~sum.rain.tot+win.rain.tot+spr.rain.tot+name+dist.m,cor=corAR1(form=~yr|name))
AIC(test1)
AIC(test2)
pR2.gls(test1)
pR2.gls(test2)
sum(spring)/sum(c(spring,summer,fall))

test<-rda(cbind(Forb.cov,G3.cov,G4.cov,Leg.cov,Wood.cov)~sum.rain.tot+win.rain.tot+spr.rain.tot+Condition(name)+Condition(dist.m))
plot(test,display=c("cn","sp"),scaling = 1)

sumr<-tapply(sum.rain.tot,yr,mean)
sprr<-tapply(spr.rain.tot,yr,mean)
winr<-tapply(win.rain.tot,yr,mean)
plot(1998:2008,sumr,type='l',col='brown')
points(1998:2008,sprr,type='l',col='green')
points(1998:2008,winr,type='l',col='grey')


spat.mods(rich,cbind(logCa,slope,northness,yr.m,dist.m),yr,name,"REML")
             [,1]
exp    19.8673793
gau    45.4263307
lin    47.0576265
ratio  26.2966336
sph    37.5875515
expN    0.0000000
gauN    1.4905260
linN    0.4936901
ratioN  0.1124233
sphN    0.2669581
mod.sit<-gls(rich~logCa+slope+northness+yr.f+dist.m,cor=corExp(form=~yr|name,n=F),meth="ML")
mod.sitn<-update(mod.sit,cor=corExp(form=~yr|name,n=T))
anova(mod.sit,mod.sitn)##nugget is worth it, now estimate with REML
mod.sit<-update(mod.sitn,meth="REML")
mod.sit<-gls(scale(rich)~scale(logCa)+scale(slope)+scale(northness)+yr.f+scale(dist.m),cor=corExp(form=~yr|name,n=T))

plot(Variogram(mod.sit),ylim=c(.2,1))
plot(Variogram(mod.sit,resType='n'),ylim=c(.2,1))
round(summary(mod.sit)$t,3)
               Value Std.Error t-value p-value
(Intercept)  139.577    27.177   5.136   0.000
logCa        -20.347     7.977  -2.551   0.011
slope         -0.689     0.817  -0.843   0.400
northness      4.030     2.092   1.926   0.055
                     Value Std.Error t-value p-value
(Intercept)         -0.632     0.176  -3.586   0.000
scale(logCa)        -0.275     0.108  -2.551   0.011
scale(slope)        -0.104     0.124  -0.843   0.400
scale(northness)     0.244     0.127   1.926   0.055

anova(mod.sit,type='m')
Denom. DF: 203 
            numDF   F-value p-value
(Intercept)     1 26.377059  <.0001
logCa           1  6.505766  0.0115
slope           1  0.710590  0.4002
northness       1  3.710967  0.0555
yr.f           10 12.811645  <.0001
dist.m          3  6.486342  0.0003

#R2s

mod.ca<-update(mod.sit,.~.-logCa)
pR2diff(mod.sit,mod.ca)
[1] 0.04699952
mod.slope<-update(mod.sit,.~.-slope)
pR2diff(mod.sit,mod.slope)
[1] 0.004875943
mod.north<-update(mod.sit,.~.-northness)
pR2diff(mod.sit,mod.north)
[1] 0.01043446

#spat.mods(rich,cbind(PDSIavg,PDSIavg^2,sum.rain.tot,win.rain.tot,spr.rain.tot,site.m,dist.m),yr,name,"REML")
spat.mods(rich,cbind(sum.rain.tot,win.rain.tot,spr.rain.tot,site.m,dist.m),yr,name,"REML")
            [,1]
exp     5.337600
gau     9.370104
lin     9.477604
ratio   5.467192
sph     9.477604
expN    2.673561
gauN    0.000000
linN   11.477604
ratioN  1.508481
sphN   11.477604
#mod.clim<-gls(rich~PDSIavg+I(PDSIavg^2)+sum.rain.tot+win.rain.tot+spr.rain.tot+name+dist.m,cor=corGaus(form=~yr|name,n=T))
mod.clim<-gls(rich~sum.rain.tot+win.rain.tot+spr.rain.tot+name+dist.m,meth="ML")
mod.clim1<-update(mod.clim,cor=corExp(form=~yr|name,n=F))
mod.climn<-update(mod.clim,cor=corGaus(form=~yr|name,n=T))
anova(mod.clim,mod.clim1,mod.climn)
##Exp w/o nugget is sig better but Guas has lowest AIC but not sig
##appears that extra parameter for nugget is not really worth it
##visually inspect
plot(Variogram(mod.clim))
plot(Variogram(mod.climn))
plot(ACF(mod.clim,resType='n'),alpha=0.05)
plot(ACF(mod.climn,resType='n'),alpha=0.05)
##the Exp model (w/o nugget seems to do best)##

mod.clim<-update(mod.clim1,meth="REML")
mod.clim<-gls(scale(rich)~scale(sum.rain.tot)+scale(win.rain.tot)+scale(spr.rain.tot)+name+scale(dist.m),cor=corExp(form=~yr|name))
round(summary(mod.clim)$t,3)
                    Value Std.Error t-value p-value
(Intercept)        89.330     5.198  17.185   0.000
sum.rain.tot       -0.192     0.037  -5.162   0.000
win.rain.tot        0.196     0.050   3.934   0.000
spr.rain.tot        0.182     0.054   3.381   0.001
                     Value Std.Error t-value p-value
(Intercept)          1.301     0.264   4.924   0.000
scale(sum.rain.tot) -0.191     0.037  -5.162   0.000
scale(win.rain.tot)  0.151     0.038   3.934   0.000
scale(spr.rain.tot)  0.137     0.041   3.381   0.001

anova(mod.clim,type='m')
Denom. DF: 194 
             numDF   F-value p-value
(Intercept)      1 295.32942  <.0001
sum.rain.tot     1  26.64880  <.0001
win.rain.tot     1  15.47409   1e-04
spr.rain.tot     1  11.43367   9e-04
name            19   8.69408  <.0001
dist.m           3  10.41317  <.0001

mod.sum<-update(mod.clim,.~.-sum.rain.tot)
pR2diff(mod.clim,mod.sum)
[1] 0.03443435
mod.win<-update(mod.clim,.~.-win.rain.tot)
pR2diff(mod.clim,mod.win)
[1] 0.01552391
mod.spr<-update(mod.clim,.~.-spr.rain.tot)
pR2diff(mod.clim,mod.spr)
[1] 0.01087412

spat.mods(rich,cbind(dist.m,site.m,yr.m),yr,name,"ML")
             [,1]
exp     6.8261730
gau    11.3022399
lin    11.4031473
ratio   6.6250050
sph    11.4031473
expN    1.3397267
gauN    0.0000000
linN   13.4031473
ratioN  0.3144198
sphN   13.4031473

mod.dist<-gls(rich~YrsOB+YrsSLB+BP5Yrs+name+yr.f,meth="ML")
mod.dist1<-update(mod.dist,cor=corExp(form=~yr|name))
mod.dist2<-update(mod.dist,cor=corGaus(form=~yr|name,n=T))
anova(mod.dist,mod.dist1,mod.dist2)
##temp autocorr needs to be corrected but somewhat abiguous if exp or guas is bets
##visual inspection
plot(Variogram(mod.dist1))
plot(Variogram(mod.dist2))
plot(ACF(mod.dist1,resType='n'),alpha=0.05)
plot(ACF(mod.dist2,resType='n'),alpha=0.05)
##neither model looks that good but Guas does remove some of the sig autocorr at short lags
mod.dist<-update(mod.dist2,meth="REML")
mod.dist<-gls(scale(rich)~scale(YrsOB)+scale(YrsSLB)+scale(BP5Yrs)+name+yr.f,cor=corGaus(form=~yr|name,n=T))
mod.dist<-gls(rich~dist.m+name+yr.f,cor=corGaus(form=~yr|name,n=T))
plot(Variogram(mod.dist))

round(summary(mod.dist)$t,3)
                    Value Std.Error t-value p-value
YrsOB               1.309     0.475   2.758   0.006
YrsSLB             -0.661     0.352  -1.877   0.062
BP5Yrs             -0.971     0.840  -1.156   0.249
                   Value Std.Error t-value p-value
(Intercept)        0.599     0.398   1.505   0.134
scale(YrsOB)       0.434     0.157   2.758   0.006
scale(YrsSLB)     -0.110     0.059  -1.877   0.062
scale(BP5Yrs)     -0.105     0.091  -1.156   0.249

anova(mod.dist,type='m')
Denom. DF: 187 
            numDF   F-value p-value
(Intercept)     1 246.19979  <.0001
YrsOB           1   7.60537  0.0064
YrsSLB          1   3.52374  0.0621
BP5Yrs          1   1.33740  0.2490
name           19   5.80048  <.0001
yr.f           10  12.87392  <.0001
Denom. DF: 187 
            numDF   F-value p-value
(Intercept)     1 246.19979  <.0001
dist.m          3   4.07308  0.0078
name           19   5.80048  <.0001
yr.f           10  12.87392  <.0001

#R2
mod.bis<-update(mod.dist,.~.-YrsOB)
pR2diff(mod.dist,mod.bis)
[1] 0.008609493
mod.bp5<-update(mod.dist,.~.-BP5Yrs)
pR2diff(mod.dist,mod.bp5)
[1] 0.003214671
mod.yrsb<-update(mod.dist,.~.-YrsSLB)
pR2diff(mod.dist,mod.yrsb)
[1] 0.003655319

##for table 2
mod.dist<-gls(rich~dist.m+name+yr.f,cor=corGaus(form=~yr|name,n=T))
anova(mod.dist,type='m')
Denom. DF: 187 
            numDF   F-value p-value
(Intercept)     1 246.19979  <.0001
dist.m          3   4.07308  0.0078
name           19   5.80048  <.0001
yr.f           10  12.87392  <.0001


mod.nam<-update(mod.dist,.~.-name)
pR2diff(mod.dist,mod.nam)
[1] 0.2501726
mod.yr<-update(mod.dist,.~.-yr.f)
pR2diff(mod.dist,mod.yr)
[1] 0.2101640
mod.dis3<-update(mod.dist,.~.-dist.m)
pR2diff(mod.dist,mod.dis3)
[1] 0.01417663



plot(Variogram(mod.dist),ylim=c(0,1))
plot(ACF(mod.dist,resType='n'),alpha=0.05)
test<-update(mod.dist,.~.-name)
spat.mods(rich,cbind(dist.m,yr.m),yr,name,"REML")#ExpN
spat.mods(rich,cbind(dist.m,site.m),yr,name,"REML")#gauN 
spat.mods(rich,cbind(yr.m,site.m),yr,name,"REML")#gauN
spat.mods(rich,site.m,yr,name,"REML")#ratioN 
spat.mods(rich,dist.m,yr,name,"REML",lin=F)#ExpN
spat.mods(rich,yr.m,yr,name,"REML")#ratioN

mod.dist<-update(mod.dist,cor=corARMA(q=2,form=~yr|name))


##mixed models for comparison##
mod.sit<-lme(rich~logCa+slope+northness+dist.m,random=~yr|name,cor=corExp(form=~yr|name,n=F))
anova(mod.sit,type='m')
            numDF denDF  F-value p-value
(Intercept)     1   196 36.78423  <.0001
logCa           1   196  8.96379  0.0031 ##196 is the listed dfs because this variable vars across years
slope           1    17  0.32056  0.5787
northness       1    17  3.58335  0.0755
dist.m          3   196 14.12575  <.0001
Ca.a<-rep(tapply(logCa,name,mean),each=11)
mod.sit<-lme(rich~Ca.a+slope+northness+dist.m,random=~yr|name,meth="ML")
mod.sit1<-lme(rich~Ca.a+slope+northness+dist.m,random=~yr|name,cor=corExp(form=~yr|name,n=F),meth="ML")
mod.sit2<-lme(rich~Ca.a+slope+northness+dist.m,random=~yr|name,cor=corExp(form=~yr|name,n=T),meth="ML")
anova(mod.sit,mod.sit1,mod.sit2)##spat mod not needed
plot(Variogram(mod.sit))
plot(Variogram(mod.sit,resType='n'))
anova(mod.sit,type='m')
            numDF denDF  F-value p-value
(Intercept)     1   197 54.22009  <.0001
Ca.a            1    16 20.28057  0.0004
slope           1    16  0.11639  0.7374
northness       1    16  2.58849  0.1272
dist.m          3   197 21.01870  <.0001


##########

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

pR2diff<-function(modfull,modsub){
 pR2.gls(modfull)$r2ML-pR2.gls(modsub)$r2ML
}

