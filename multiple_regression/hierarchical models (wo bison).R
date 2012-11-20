##Purpose: to perform multiple regression analyses on TGPP species richness
##(I): examination of each spatial grain and how variance is partitioned between sites and years
##(II): examinination of enviornmental variables on richness - includes several ways of looking at
## partial residuals and bootstrapped confidence intervals using the 'termplot' and 'boot' functions

endat<-read.table('C:/Users/dmcglinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
source('C:/Users/dmcglinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/multiple regression/bootstrap helper functions & permutation test for model coef.R')
#endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
library(nlme)
library(car)
names(endat)
attach(endat)

##examining within and between corner/plot variance through time
##begin working only at the 100 m2 scale so no corners
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)
##within and between plot variance can be quantified with aov or with lme
mod1.lm<-aov(rich~Error(name))
summary(mod1.lm)
Error: name
          Df  Sum Sq Mean Sq F value Pr(>F)
Residuals 19 19934.1  1049.2               

Error: Within
           Df  Sum Sq Mean Sq F value Pr(>F)
Residuals 200 15959.6    79.8    
##Within plot variance = 79.8
##Between plot variance is next line
## (MSE accounted by id - MSE not accounted by id)/n [n in our case is 11 years)
(1049.2 - 79.8)/11
#[1] 88.12727

##this will match the output of the mixed effect algo
mod1.lme<-lme(rich~1,random=~1|name,method="REML")
mod1.lme
 Fixed: rich ~ 1 
(Intercept) 
   76.23636 

Random effects:
 Formula: ~1 | name
        (Intercept) Residual
StdDev:    9.129898 8.932983
##square the stdevs to get the variances found with aov
9.387^2
8.933^2
##interclass correlation pg308 in DAAG
(9.1298^2)/(9.1298^2+8.933^2)
#51% of the variance in a particular plot is explained by plot to plot diffs
#Or 51% is the proporation of residual variance explained by diffs between sites
##note if method in lme was set to "ML" there would be some bias in these estimates
intervals(mod1.lme)
Approximate 95% confidence intervals
 Fixed effects:
               lower     est.    upper
(Intercept) 71.93017 76.23636 80.54256
attr(,"label")
[1] "Fixed effects:"
 Random Effects:
  Level: name 
                   lower     est.    upper
sd((Intercept)) 6.653638 9.387436 13.24448
 Within-group standard error:  ##not actually std error (i don't think)
   lower     est.    upper 
8.099095 8.932983 9.852728 
##Take home message, there is approximately equal variance within and between,
##their intervals overlap fairly well
detach(endatlv1)

##Function for extracting the proporation of variance due to plot and corner respectively 
var.ext<-function(lmemod){
  if(dim(VarCorr(lmemod))[1]>2){
   vplv1<-as.numeric(VarCorr(lmemod)[2,1])
   vplv2<-as.numeric(VarCorr(lmemod)[4,1])
   vpresid<-as.numeric(VarCorr(lmemod)[5,1])
   c("plot"=vplv1/(vplv1+vplv2+vpresid),"cor.in.plot"=vplv2/(vplv2+vpresid))
  }
  else{
   vplv1<-as.numeric(VarCorr(lmemod)[1,1])
   vresids<-as.numeric(VarCorr(lmemod)[1,2])
   c('plot'=vplv1/(vplv1+vresids))
}}
var.ext(mod1.lme)

##consider lv2 richness
attach(endat)
endatlv2<-endat[lv==2,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv2)
###now we have corners within plots
mod2.lm<-aov(rich~Error(name/cor))
summary(mod2.lm)
mod2.lme<-lme(rich~1,random=~1|name/cor)
summary(mod2.lme)
Random effects:
 Formula: ~1 | name
        (Intercept)
StdDev:    6.256476

 Formula: ~1 | cor %in% name
        (Intercept) Residual
StdDev:    2.339644 7.259241
##to get each of the above estimates look to summary(mod1.lm)
##see pg 321 in DAAG for explanation
sqrt(56)
sqrt((107.36-56)/11)
sqrt((1835-107.36)/(4*11))

var.ext(mod2.lme)
      plot cor.in.plot 
0.40223907  0.09410151 
##rest of lme output

##so at the 10m2 scale,
##within each corner (through time) sd of 7.25 species
##this is greater than the between corner sd of 2.33 species
##betwee the plot there is within sd of 6.25 sp

detach(endatlv2)

##consider lv3 richness
attach(endat)
endatlv3<-endat[lv==3,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv3)
###now we have corners within plots
mod3.lm<-aov(rich~Error(name/cor))
summary(mod3.lm)
mod3.lme<-lme(rich~1,random=~1|name/cor)
summary(mod3.lme)
Random effects:
 Formula: ~1 | name
        (Intercept)
StdDev:    3.220817

 Formula: ~1 | cor %in% name
        (Intercept) Residual
StdDev:    2.235098 4.792486
##to get each of the above estimates look to summary(mod1.lm)
##see pg 321 in DAAG for explanation
sqrt(25)
sqrt((102.8-25)/11)
sqrt((534.4-102.8)/(4*11))

var.ext(mod3.lme)
      plot cor.in.plot 
 0.2705896   0.1786489
##so at the 1m2 scale,
##within each corner (through time) sd of 4.79 species
##this is greater than the between corner sd of 2.235 species
##betwee the plot there is within sd of 3.22 sp

##how much of variance, for an individual yr, is explained by yr-to-yr diffs
##for yrs in the same corner ,pg325 in DAAG
4.792^2/(4.792^2+2.235^2)
##  0.8213342 so 82% of the variance is expl by yr-to-yr diffs if in diff corner same plot
##for diff corners in different plots
(4.792^2)/(4.792^2 + 2.235^2 + 3.221^2)
##0.5990417 so 60% of the variance is expl by yr-to-yr diffs in diff corners in diff plot

detach(endatlv3)
#################
##consider lv4 richness
attach(endat)
endatlv4<-endat[lv==4,] ##to only work at the .01 m2 scale
detach(endat)
attach(endatlv4)
###now we have corners within plots
mod4.lm<-aov(rich~Error(name/cor))
summary(mod4.lm)
mod4.lme<-lme(rich~1,random=~1|name/cor)
summary(mod4.lme)
Random effects:
 Formula: ~1 | name
        (Intercept)
StdDev:     1.44672

 Formula: ~1 | cor %in% name
        (Intercept) Residual
StdDev:    1.115273 3.057777
##
var.ext(mod4.lme)
     plot cor.in.plot 
0.1649740   0.1174112 

detach(endatlv4)
##consider lv5 richness
attach(endat)
endatlv5<-endat[lv==5,] ##to only work at the .01 m2 scale
detach(endat)
attach(endatlv5)
###now we have corners within plots
mod5.lm<-aov(rich~Error(name/cor))
summary(mod5.lm)
mod5.lme<-lme(rich~1,random=~1|name/cor)
summary(mod5.lme)
Random effects:
 Formula: ~1 | name
        (Intercept)
StdDev:   0.3735333

 Formula: ~1 | cor %in% name
        (Intercept) Residual
StdDev:   0.4619026 1.639152
##
var.ext(mod5.lme)
       plot cor.in.plot 
 0.04590161  0.07356598 
##Take home message is that the 
##for lv3 18%corner/37%plot = 
##for lv5 7%corner/4%plot
##So the take home is that the importance of a corner rel to a plot inc as you decrease grain
##Q: does this change occur regularly?

prop.vars<-rbind(c(var.ext(mod1.lme),NA),var.ext(mod2.lme),var.ext(mod3.lme),var.ext(mod4.lme),var.ext(mod5.lme))
colnames(prop.vars)<-c("plot","cor.in.plot")
rownames(prop.vars)<-c("lv1","lv2",'lv3','lv4','lv5')
prop.vars
          plot cor.in.plot
lv1 0.90372986          NA
lv2 0.40223907  0.09410151
lv3 0.27058964  0.17864887
lv4 0.16497403  0.11741121
lv5 0.04590161  0.07356598

plot(prop.vars,type='l')
##proportion of var expl
var.ext(mod1.lme)

plot(2:-2,prop.vars[,1],type='o',xlab='',ylab='')
points(2:-2,prop.vars[,2],col='red',type='o')
legend('topleft',c('between plot','between corner'),col=c('black','red'),lty=1,pch=1,bty='n')
##var compoents
VarCorr(mod2.lme)[5,1]
plot.var<-as.numeric(c(VarCorr(mod2.lme)[2,1],VarCorr(mod3.lme)[2,1],VarCorr(mod4.lme)[2,1],VarCorr(mod5.lme)[2,1]))
cor.var<-as.numeric(c(VarCorr(mod2.lme)[4,1],VarCorr(mod3.lme)[4,1],VarCorr(mod4.lme)[4,1],VarCorr(mod5.lme)[4,1]))
res.var<-as.numeric(c(VarCorr(mod2.lme)[5,1],VarCorr(mod3.lme)[5,1],VarCorr(mod4.lme)[5,1],VarCorr(mod5.lme)[5,1]))

plot(2:-2,(c(as.numeric(VarCorr(mod1.lme)[1,1]),plot.var)),type='o',ylim=c(0,100))
points(1:-2,(cor.var),type='o',col='red')
points(2:-2,(c(as.numeric(VarCorr(mod1.lme)[2,1]),res.var)),type='o',col='blue')
legend('topright',c('between plot','between corner','residual'),col=c('black','red','blue'),lty=1,pch=1,bty='n')

plot(2:-2,log(sqrt(c(as.numeric(VarCorr(mod1.lme)[1,1]),plot.var))),type='o',ylim=c(-2,4),xlab='',ylab='')
points(1:-2,log(sqrt(cor.var)),type='o',col='red')
points(2:-2,log(sqrt(c(as.numeric(VarCorr(mod1.lme)[2,1]),res.var))),type='o',col='blue')
legend('topright',c('between plot','between corner','residual'),col=c('black','red','blue'),lty=1,pch=1,bty='n')

#########################################################
#########################################################
##Now its time to build some linear fixed effect models to do var part with
##work only with level 1
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)
c1<-as.factor(cor)
p1<-as.factor(plot)
yr.f<-factor(yr)
syr<-scale(yr,scale=F)
syr2<-syr^2
syr3<-syr^3
srain<-scale(rain.tot)
stemp<-scale(temp.avg)
psrain<-scale(prain.tot)
pstemp<-scale(ptemp.avg)
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
#############################################################
##begin with gls analysis then move to mixed effects model
##1st Goal to address to examine partial regression plots for dist
##after factoring out plot and yr
gls.full<-gls(rich~-1+name+yr.f+YrsOB+BP5Yrs+YrsSLB,meth="ML")
plot(ACF(gls.full,maxLag=5),alpha=0.05)
gls.full2<-update(gls.full,cor=corAR1(form=~yr|name))
plot(ACF(gls.full2,maxLag=5,resType='n'),alpha=0.05)
anova(gls.full,gls.full2)
          Model df      AIC      BIC    logLik   Test  L.Ratio p-value
gls.full      1 34 1486.125 1601.508 -709.0625                        
gls.full2     2 35 1483.096 1601.873 -706.5481 1 vs 2 5.028723  0.0249
##to get sweet partial plots
coef(gls.full2)
mod.lm<-lm(rich~name+yr.f+YrsOB+BP5Yrs+YrsSLB)
summary(mod.lm)
YrsOB               1.4305     0.3153   4.537 1.02e-05 ***
BP5Yrs             -2.0217     0.7392  -2.735 0.006839 ** 
YrsSLB             -1.0145     0.3375  -3.006 0.003007 **
mod.lm<-lm(scale(rich)~name+yr.f+syob+sbp5+syslb)
summary(mod.lm)
syob               0.47369    0.10440   4.537 1.02e-05 ***
sbp5              -0.21869    0.07996  -2.735 0.006839 ** 
syslb             -0.16946    0.05637  -3.006 0.003007 **

library(hier.part)
#the general partitioning (site,yr,mang)
gen.part<-hier.part(rich,data.frame(name=name,yr=yr.f,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu')
##unconstrained test
gen.test<-rand.hp(rich,data.frame(name=name,yr=yr.f,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu',num.reps=499);alarm()
##constrained to obs rows
gen.test2<-nrand.hp(rich,data.frame(name=name,yr=yr.f,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu',num.reps=499);alarm()
pval<-matrix(NA,ncol=2,nrow=5)##non-parametric 1-sided p-values
for(i in 1:5){
  pval[i,1]<- sum(ifelse(gen.test$Irands[,i]>=gen.test$Irands[1,i],1,0))/500
  pval[i,2]<- sum(ifelse(gen.test2$Irands[,i]>=gen.test2$Irands[1,i],1,0))/500
}
rownames(pval)<-c('site','yr','yrsbison','BP5','YrsSLB')
pval
          [,1]  [,2]
site     0.002 0.002
yr       0.002 0.004
yrsbison 0.002 0.002
BO5      0.008 0.002
YrsSLB   0.116 0.074


##hierarchical plot with squared partial correlations
##for plot, yr, and mang vars
mod.lm<-lm(scale(rich)~name+yr.f+syob+sbp5+syslb)
partial.r2<-coef(mod.lm)[31:33]^2

holder<-rbind(gen.part$IJ[,3],c(NA,NA,partial.r2))
colnames(holder)<-NULL
barplot(holder,beside=T,ylim=c(0,.6),col=c('grey10','black'),density=c(NA,12),angle=c(NA,45),cex.axis=2)
par(new=T)
holder<-rbind(gen.part$IJ[,1],rep(NA,5))
barplot(holder,col='grey90',beside=T,ylim=c(0,.6),axes=F)

##the specific partitioning (soil,topo,clim,mang)
spe.part<-hier.part(rich,data.frame(logCa=logCa,slope=slope,north=northness,pdsi=PDSIavg,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu')
spe.part<-hier.part(rich,data.frame(logCa=logCa,slope=slope,north=northness,rain=rain.tot,prain=prain.tot,temp=temp.avg,ptemp=ptemp.avg,pdsi=PDSIavg,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu')
spe.part<-hier.part(rich,data.frame(logCa=logCa,slope=slope,north=northness,rain=rain.tot,prain=prain.tot,temp=temp.avg,ptemp=ptemp.avg,pdsi=PDSIavg,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu')
spe.part<-hier.part(rich,data.frame(logCa=logCa,slope=slope,north=northness,sumr=sum.rain.tot,winr=win.rain.tot,sprr=spr.rain.tot,sumt=sum.temp,wint=win.temp,sprt=spr.temp,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu')



system.time(spe.test<-rand.hp(rich,data.frame(logCa=logCa,slope=slope,north=northness,pdsi=PDSIavg,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu',num.reps=499));alarm()
  user  system elapsed 
 469.42    0.79  484.30 
##constrained to obs rows
system.time(spe.test2<-nrand.hp(rich,data.frame(logCa=logCa,slope=slope,north=northness,pdsi=PDSIavg,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu',num.reps=499));alarm()
  user  system elapsed 
 503.66    0.31  515.41 
pval<-matrix(NA,ncol=2,nrow=7)##non-parametric 1-sided p-values
for(i in 1:7){
  pval[i,1]<- sum(ifelse(spe.test$Irands[,i]>=spe.test$Irands[1,i],1,0))/nrow(spe.test2$Irands)
  pval[i,2]<- sum(ifelse(spe.test2$Irands[,i]>=spe.test2$Irands[1,i],1,0))/nrow(spe.test2$Irands)
}
rownames(pval)<-c('logCa','slope','north','pdsi','yrsbison','BP5','YrsSLB')
pval
          [,1]  [,2]
logCa    0.002 0.002
slope    0.028 0.032
north    0.006 0.006
pdsi     0.048 0.062
yrsbison 0.002 0.002
BP5      0.004 0.002
YrsSLB   0.012 0.008

##hierarchical plot with squared partial correlations
##for plot, yr, and mang vars
mod.lm<-lm(scale(rich)~scale(logCa)+scale(slope)+scale(northness)+scale(PDSIavg)+syob+sbp5+syslb)
cbind(round(coef(mod.lm),2))
cbind(round(coef(mod.lm),2))
                  [,1]
(Intercept)       0.00
scale(logCa)     -0.48
scale(slope)     -0.09
scale(northness)  0.17
scale(PDSIavg)   -0.04
syob              0.42
sbp5             -0.21
syslb            -0.29

partial.r2<-coef(mod.lm)[-1]^2
signs<-ifelse(coef(mod.lm)[-1]>0,"+","-")


holder<-rbind(spe.part$IJ[,3],partial.r2)
colnames(holder)<-NULL
##add directions
barplot(holder,beside=T,space=c(0,1),ylim=c(0,.25),col=c('grey10','black'),density=c(NA,12),angle=c(NA,45),cex.axis=2)
points(seq(2.5,20.5,3),holder[2,]+.01,pch=signs,cex=2)
par(new=T)
holder<-rbind(spe.part$IJ[,1],rep(NA,7))
barplot(holder,col='grey90',beside=T,space=c(0,1),ylim=c(0,.6),axes=F)


#######################################################
mod.lm<-gls(rich~-1+name+yr.f+YrsOB+BP5Yrs+YrsSLB)
anova(mod.lm,type='m')
Denom. DF: 186 
       numDF  F-value p-value
name      20 68.38627  <.0001
yr.f      10 10.72236  <.0001
YrsOB      1 19.97715  <.0001
BP5Yrs     1  7.31896  0.0075
YrsSLB     1  9.00689  0.0031
bison      1  0.01750  0.8949
anova(gls.full2,type='m')
Denom. DF: 186 
       numDF  F-value p-value
name      20 48.67666  <.0001
yr.f      10 10.17113  <.0001
YrsOB      1 12.78412  0.0004
bison      1  0.01636  0.8984
BP5Yrs     1  3.88856  0.0501
YrsSLB     1  6.04623  0.0149

round(summary(gls.full2)$tTable,3)
YrsOB              1.441     0.348   4.144   0.000
bison             -4.495     2.143  -2.097   0.037
BP5Yrs            -2.063     0.794  -2.598   0.010
YrsSLB            -0.988     0.348  -2.844   0.005
##so std errors are ajusted for autocorr

mod.lm<-lm(rich~name+yr.f+YrsOB+BP5Yrs+YrsSLB)
summary(mod.lm)
par(mfrow=c(1,3))
termplot(mod.lm,terms="YrsOB",xlab='',ylab='',partial=TRUE,se=T,smooth=panel.smooth,span.smth=1/3,lwd.term=1.5,lwd.se=2,col.term='blue',col.se='blue', col.res = 1,ylim=c(-28,28))
termplot(mod.lm,terms="YrsSLB",xlab='',ylab='',partial=TRUE,se=T,smooth=panel.smooth,span.smth=1/3,lwd.term=1.5,lwd.se=2,ylim=c(-28,28))
termplot(mod.lm,terms="BP5Yrs",xlab='',ylab='',partial=TRUE,se=T,smooth=panel.smooth,span.smth=1/3,lwd.term=1.5,lwd.se=2,ylim=c(-28,28))

termplot(mod.lm,terms="YrsOB",partial=TRUE,se=T,smooth=panel.smooth,span.smth=1/3,lwd.term=1.5,lwd.se=2,ylim=c(-28,28))

plot(yr,richexotic/rich)
mod<-lm(richexotic/rich~name+yr.f+YrsOB+BP5Yrs+YrsSLB)
summary(mod)
##mirrors results on just richness
termplot(mod)

##in a binomical modeling context, each species is treated as an observation
reper<-function(x,y){c(rep(1,x),rep(0,y-x))}
resp<-unlist(mapply(reper,richexotic,rich))
yr.rep<-rep(yr,rich)
yr.repf<-as.factor(rep(yr,rich))
site.rep<-as.factor(rep(name,rich))
YrsOB.rep<-rep(YrsOB,rich)
YrsSLB.rep<-rep(YrsSLB,rich)
BP5Yrs.rep<-rep(BP5Yrs,rich)

mod<-glm(resp~yr.repf+site.rep+YrsOB.rep+YrsSLB.rep+BP5Yrs.rep,family='binomial')
summary(mod) ##management variables are not significant




##the below recreates the termplot for one of the vars
##drop YrsSLB
resid.calc<-rich-as.matrix(cbind(endatlv1[,8:27],endatlv1[,29:38],YrsOB,BP5Yrs,bison))%*%coef(gls.full2)[-33]
#resid.calc<-rich-as.matrix(cbind(endatlv1[,8:27],endatlv1[,29:38],YrsOB,BP5Yrs,bison))%*%coef(mod.lm)[-33]
resid.calc<-resid.calc-mean(resid.calc)
##an alternative way is with the function model.matrix
##this function is sometimes glitchy though
## and a call to the function formula must first be used (see below for example)

####bootstrap confidence band
####fixed X values
##I'm not sure if I need to multipe the coef by a weights matrix with autocorr errors
library(boot)

###bison#############
resids<- rich - model.matrix(formula(gls.full2))[,-34]%*%coef(gls.full2)[-34]
resids.c<-resids-mean(resids)
mod<-gls(resids.c~bison,cor=corAR1(form=~yr|name),meth="ML")
par(mfrow=c(2,2));plot(lm(resids.c~bison))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
summary(mod)
fit<-predict(mod)
e<-residuals(mod)
e.c<- e-mean(e)
X <- bison
pbison <- boot(endatlv1, boot.pred, 1999);alarm()
plot(pbison,index=1)##shows bootstrap dist for just one fitted value

##my home brewed permutation test
cbison <- perm.coef(resp=resids.c,indep=bison,R=1999);alarm()
##shuffling only within plots
cbison.pl <- perm.coef(y=resids.c,x=bison,strata=name,R=1999);alarm()
##p.value
1-sum(c(cbison.pl,coef(mod)[2])>coef(mod)[2])/length(c(cbison.pl,coef(mod)[2]))
[1] 0.001
cbison.tor <- perm.coef(y=resids.c,x=bison,strata=name,R=1999,torus=TRUE);alarm()
1-sum(c(cbison.tor,coef(mod)[2])>coef(mod)[2])/length(c(cbison.tor,coef(mod)[2]))
[1] 0.004
par(mfrow=c(1,2))
hist(cbison.pl,breaks=50)
hist(cbison.tor,breaks=50)

boot.ci(cbison,index=2)
##generate con.intervals for normal, percent, and basic
bison.ci<-CI(pbison)
##confi band
plot(bison,resids.c,col='black',pch=19,ylab='',xlab='')
lines(lowess(bison,resids.c),col='red',lwd=2,lty=2)
newx<-data.frame(bison=pretty(bison,500))
lines(newx[,1],predict(mod,newdata=newx),lwd=3)
confi.band(bison,bison.ci,which='perc',lwd=3,col='darkgrey',lty=1)

##YrsSLB#############
resids<- rich - model.matrix(formula(gls.full2))[,-33]%*%coef(gls.full2)[-33]
resids.c<-resids-mean(resids)
mod<-gls(resids.c~YrsSLB,cor=corAR1(form=~yr|name),meth="ML")
par(mfrow=c(2,2));plot(lm(resids.c~YrsSLB))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
fit<-predict(mod)
e<-residuals(mod)
e.c<- e-mean(e)
X <- YrsSLB

pYrsSLB <- boot(endatlv1, boot.pred, 1999);alarm()
plot(pYrsSLB)

##generate con.intervals for normal, percent, and basic
YrsSLB.ci<-CI(pYrsSLB)
##confi band
plot(YrsSLB,resids.c,col='black',pch=19,ylab='',xlab='')
lines(lowess(YrsSLB,resids.c),col='red',lwd=2,lty=2)
newx<-data.frame(YrsSLB=pretty(YrsSLB,500))
lines(newx[,1],predict(mod,newdata=newx),lwd=3)
confi.band(YrsSLB,YrsSLB.ci,which='perc',lwd=3,col='darkgrey',lty=1)

##BP5Yrs#############
resids<- rich - model.matrix(formula(gls.full2))[,-32]%*%coef(gls.full2)[-32]
resids.c<-resids-mean(resids)
mod<-gls(resids.c~BP5Yrs,cor=corAR1(form=~yr|name),meth="ML")
par(mfrow=c(2,2));plot(lm(resids.c~BP5Yrs))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
summary(mod)
fit<-predict(mod)
e<-residuals(mod)
e.c<- e-mean(e)
X <- BP5Yrs

pBP5Yrs <- boot(endatlv1, boot.pred, 1999);alarm()
plot(pBP5Yrs)

##generate con.intervals for normal, percent, and basic
BP5Yrs.ci<-CI(pBP5Yrs)
##confi band
plot(BP5Yrs,resids.c,col='black',pch=19,ylab='',xlab='')
lines(lowess(BP5Yrs,resids.c),col='red',lwd=2,lty=2)
newx<-data.frame(BP5Yrs=pretty(BP5Yrs,500))
lines(newx[,1],predict(mod,newdata=newx),lwd=3)
confi.band(BP5Yrs,BP5Yrs.ci,which='perc',lwd=3,col='darkgrey',lty=1)

##YrsOB#############
resids<- rich - model.matrix(formula(gls.full2))[,-31]%*%coef(gls.full2)[-31]
resids.c<-resids-mean(resids)
mod<-gls(resids.c~YrsOB,cor=corAR1(form=~yr|name),meth="ML")
par(mfrow=c(2,2));plot(lm(resids.c~YrsOB))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
summary(mod)
fit<-predict(mod)
e<-residuals(mod)
e.c<- e-mean(e)
X <- YrsOB

pYrsOB <- boot(endatlv1, boot.pred, 1999);alarm()
plot(pYrsOB)

##examine case bootstrapping
boot.data<-data.frame(y=resids.c,x=YrsOB)
pYrsOB.c<- boot(boot.data, boot.pred.case, strata=name, 999);alarm()
cYrsOB.c<- boot(boot.data, boot.coef.case, strata=name, 999);alarm()
plot(pYrsOB.c,index=1)
plot(cYrsOB.c,index=2)
boot.ci(cYrsOB.c,index=2)
##case bootstrapping, results in wide predicted regression line, but consitent estimates of slope with fixed x shuffling
##it makes sense that case reshuffle is not appropriate for confid band
##you need the whole range of x values represented in each random step

##generate con.intervals for normal, percent, and basic
YrsOB.ci<-CI(pYrsOB)
YrsOB.ci.c<-CI(pYrsOB.c)

##confi band
plot(YrsOB,resids.c,col='black',pch=19,ylab='',xlab='')
lines(lowess(YrsOB,resids.c),col='red',lwd=2,lty=2)
newx<-data.frame(YrsOB=pretty(YrsOB,500))
lines(newx[,1],predict(mod,newdata=newx),lwd=3)
confi.band(YrsOB,YrsOB.ci,which='norm',col='darkgrey',lwd=3,lty=1)
##confi.band(YrsOB,YrsOB.ci.c,lwd=3,lty=1)
##the case bootstraps were 

#############################################
##now perform mixed model analysis to show similar results
gls.mod1<-gls(rich~-1+name+yr.f+dist.m,cor=corAR1(form=~yr|name),method="ML")
lme.mod1<-lme(rich~-1+yr.f+dist.m,random=~1|name,cor=corAR1(form=~yr|name),method="ML")
round(summary(lme.mod1)$tTable,3)
anova(lme.mod1,type='m')
anova(lme.mod1,gls.mod1)
##so gls with plot as a fixed effect is better but may be overfitted
plot(ACF(lme.mod1,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(lme.mod1,maxLag=5,resType='n'),alpha=0.05)
plot(lme.mod1)
qq.plot(lme.mod1$resid)
lme.mod2<-update(lme.mod1,weights=varPower())
anova(lme.mod1,lme.mod2)

##diagnostic plots
##assumption 1
mod.lme<-lme.mod1
##Are the within group errors centered around zero?
plot(mod.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
plot(predict(mod.lme),resid(mod.lme,type='p'))
lines(lowess(predict(mod.lme),resid(mod.lme,type='p')),lwd=2,col='red')
abline(h=0)
##look at same assumption on a group basis
plot(mod.lme,name~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod.lme,resid(.,type='p')~fitted(.)|name,id=0.05,adj=-0.3)
##examine difference in response and fitted 
plot(mod.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qq.plot(resid(mod.lme))
qqnorm(mod.lme,~resid(.)|name)
##assumption 2
##qqnorm() on ranef(.) for checking marginal normality and identifyin outliers
##pairs() on ranef(.) for id outliers and cheing assump of homogenity of ranef covariance matrix
qqnorm(mod.lme,~ranef(.,level=1),id=0.10,cex=.7)
qq.plot(ranef(mod.lme,level=1))
##diagnostic plots indicate
intervals(mod.lme)
##graphics!!#########################
##partial regression coeff plot with ranef
resids<- rich - model.matrix(formula(mod.lme))%*%fixef(mod.lme)
plot(residuals(mod.lme,level=0,type='r'),resids)
plot(residuals(mod.lme,level=0,type='n'),resids)
##bison##
resids<- rich - model.matrix(formula(mod.lme))[,-15]%*%fixef(mod.lme)[-15]
resids.c<-resids-mean(resids)
mod<-lme(resids.c~bison,random=~1|name,cor=corAR1(form=~yr|name),meth="ML")
par(mfrow=c(2,2));plot(lm(resids.c~bison))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
summary(mod)
fit<-predict(mod,level=0)
e<-residuals(mod,level=0)
e.c<- e-mean(e)
X <- bison
pbison.lme <- boot(endatlv1, boot.pred.lme, 499 ,strata=name);alarm()
plot(pbison.lme,index=1)##shows bootstrap dist for just one fitted value
##generate con.intervals for normal, percent, and basic
bison.ci<-CI(pbison.lme)

plot(bison,resids.c,pch=19,xlab='',ylab='')
lines(lowess(bison,resids.c),col='red',lty=2,lwd=3)
newx<-data.frame(bison=pretty(bison,500))
lines(newx[,1],predict(mod,newdata=newx,lev=0),lwd=3)
uni.name<-rownames(ranef(mod))
for(i in 1:20){
 lines(c(min(bison[name==uni.name[i]]),max(bison[name==uni.name[i]])),
       c(min(bison[name==uni.name[i]])*fixef(mod)[2]+fixef(mod)[1]+ranef(mod)[i,1],
         max(bison[name==uni.name[i]])*fixef(mod)[2]+fixef(mod)[1]+ranef(mod)[i,1]),
        type='l',col='grey',lwd=1)
}
confi.band(bison,bison.ci,which='perc',lwd=3,col='darkgrey',lty=2)

##YrsSLB##
resids<- rich - model.matrix(formula(mod.lme))[,-14]%*%fixef(mod.lme)[-14]
resids.c<-resids-mean(resids)
mod<-lme(resids.c~YrsSLB,random=~1|name,cor=corAR1(form=~yr|name),meth="ML")
par(mfrow=c(2,2));plot(lm(resids.c~YrsSLB))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
summary(mod)
fit<-predict(mod,level=0)
e<-residuals(mod,level=0)
e.c<- e-mean(e)
X <- YrsSLB
pYrsSLB.lme <- boot(endatlv1, boot.pred.lme, 499, strata=name);alarm()
plot(pYrsSLB.lme,index=1)##shows bootstrap dist for just one fitted value
##generate con.intervals for normal, percent, and basic
YrsSLB.ci<-CI(pYrsSLB.lme)

plot(YrsSLB,resids.c,pch=19,xlab='',ylab='')
lines(lowess(YrsSLB,resids.c),col='red',lty=2,lwd=3)
newx<-data.frame(YrsSLB=pretty(YrsSLB,500))
lines(newx[,1],predict(mod,newdata=newx,lev=0),lwd=3)
uni.name<-rownames(ranef(mod))
for(i in 1:20){
 lines(c(min(YrsSLB[name==uni.name[i]]),max(YrsSLB[name==uni.name[i]])),
       c(min(YrsSLB[name==uni.name[i]])*fixef(mod)[2]+fixef(mod)[1]+ranef(mod)[i,1],
         max(YrsSLB[name==uni.name[i]])*fixef(mod)[2]+fixef(mod)[1]+ranef(mod)[i,1]),
        type='l',col='grey',lwd=1)
}
confi.band(YrsSLB,YrsSLB.ci,which='perc',lwd=3,col='darkgrey',lty=2)

##BP5Yrs##
resids<- rich - model.matrix(formula(mod.lme))[,-13]%*%fixef(mod.lme)[-13]
resids.c<-resids-mean(resids)
mod<-lme(resids.c~BP5Yrs,random=~1|name,cor=corAR1(form=~yr|name),meth="ML")
par(mfrow=c(2,2));plot(lm(resids.c~BP5Yrs))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
summary(mod)
fit<-predict(mod,level=0)
e<-residuals(mod,level=0)
e.c<- e-mean(e)
X <- BP5Yrs
pBP5Yrs.lme <- boot(endatlv1, boot.pred.lme, 499, strata=name);alarm()
plot(pBP5Yrs.lme,index=1)##shows bootstrap dist for just one fitted value
##generate con.intervals for normal, percent, and basic
BP5Yrs.ci<-CI(pBP5Yrs.lme)

plot(BP5Yrs,resids.c,pch=19,xlab='',ylab='')
lines(lowess(BP5Yrs,resids.c),col='red',lty=2,lwd=3)
newx<-data.frame(BP5Yrs=pretty(BP5Yrs,500))
lines(newx[,1],predict(mod,newdata=newx,lev=0),lwd=3)
uni.name<-rownames(ranef(mod))
for(i in 1:20){
 lines(c(min(BP5Yrs[name==uni.name[i]]),max(BP5Yrs[name==uni.name[i]])),
       c(min(BP5Yrs[name==uni.name[i]])*fixef(mod)[2]+fixef(mod)[1]+ranef(mod)[i,1],
         max(BP5Yrs[name==uni.name[i]])*fixef(mod)[2]+fixef(mod)[1]+ranef(mod)[i,1]),
        type='l',col='grey',lwd=1)
}
confi.band(BP5Yrs,BP5Yrs.ci,which='perc',lwd=3,col='darkgrey',lty=2)

##YrsOB##
resids<- rich - model.matrix(formula(mod.lme))[,-12]%*%fixef(mod.lme)[-12]
resids.c<-resids-mean(resids)
mod<-lme(resids.c~YrsOB,random=~1|name,cor=corAR1(form=~yr|name),meth="ML")
par(mfrow=c(2,2));plot(lm(resids.c~YrsOB))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
summary(mod)
fit<-predict(mod,level=0)
e<-residuals(mod,level=0)
e.c<- e-mean(e)
X <- YrsOB
pYrsOB.lme <- boot(endatlv1, boot.pred.lme, 499, strata=name);alarm()
plot(pYrsOB.lme,index=1)##shows bootstrap dist for just one fitted value
##generate con.intervals for normal, percent, and basic
YrsOB.ci<-CI(pYrsOB.lme)

plot(YrsOB,resids.c,pch=19,xlab='',ylab='')
lines(lowess(YrsOB,resids.c),col='red',lty=2,lwd=3)
newx<-data.frame(YrsOB=pretty(YrsOB,500))
lines(newx[,1],predict(mod,newdata=newx,lev=0),lwd=3)
uni.name<-rownames(ranef(mod))
for(i in 1:20){
 lines(c(min(YrsOB[name==uni.name[i]]),max(YrsOB[name==uni.name[i]])),
       c(min(YrsOB[name==uni.name[i]])*fixef(mod)[2]+fixef(mod)[1]+ranef(mod)[i,1],
         max(YrsOB[name==uni.name[i]])*fixef(mod)[2]+fixef(mod)[1]+ranef(mod)[i,1]),
        type='l',col='grey',lwd=1)
}
confi.band(YrsOB,YrsOB.ci)
confi.band(YrsOB,YrsOB.ci,which='perc',lwd=3,col='darkgrey',lty=2)

##parabola for YrsOB?
mod.lme<-update(mod.lme,method="ML")
mod.lme2<-update(mod.lme,.~.+I(YrsOB^2))
anova(mod.lme,mod.lme2)
mod.gls<-gls(rich~-1+name+yr.f+dist.m,cor=corAR1(form=~yr|name),method="ML")
mod.gls2<-update(mod.gls,.~.+I(YrsOB^2))
anova(mod.gls,mod.gls2)
##a parabola is better but based on a lot of scatter and potentially few indep obs

###between plots gls analysis#########
##create avg ca variable
logCa.avg<-rep(tapply(logCa,name,mean),each=11)
gls.yr<-gls(rich~yr.f,cor=corAR1(form=~yr|name),meth="ML")
resids<-residuals(gls.yr,type='n')
gls.site<-gls(resids~YrsOB+YrsSLB+BP5Yrs+logCa.avg+northness,cor=corAR1(form=~yr|name),meth="ML")
round(summary(gls.site)$tTable,3)
plot(gls.site)
plot(ACF(gls.site,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(gls.site,maxLag=5,resType='n'),alpha=0.05)
gls.site2<-update(gls.site,cor=corARMA(p=2,form=~yr|name))
plot(ACF(gls.site2,maxLag=5,resType='n'),alpha=0.05)
anova(gls.site,gls.site2)
round(summary(gls.site2)$tTable,3)
gls.site3<-update(gls.site,.~.-bison-northness)
anova(gls.site,gls.site3)
round(summary(gls.site3)$tTable,3)
gls.site4<-update(gls.site3,.~.-YrsSLB)
anova(gls.site3,gls.site4)
round(summary(gls.site4)$tTable,3)
plot(ACF(gls.site4,maxLag=5,resType='n'),alpha=0.05)
gls.site5<-update(gls.site4,cor=corAR1(form=~yr|name))
anova(gls.site4,gls.site5)
##the 2 autoreg are still better
plot(gls.site4)
plot(predict(gls.site4),resids)
abline(0,1)
mod.lm<-lm(resids~logCa.avg+YrsOB)
par(mfrow=c(1,2))
termplot(mod.lm,partial=TRUE,se=T,smooth=panel.smooth,span.smth=1/3,lwd.term=1.5,lwd.se=2)
########
###within plots gls analysis#########
gls.name<-gls(rich~name,cor=corAR1(form=~yr|name),meth="ML")
resids<-residuals(gls.name,type='n')
#clim.mat<-cbind(spr.rain.tot,win.rain.tot,sum.rain.tot,spr.temp,win.temp,sum.temp)
#pclim.mat<-cbind(pspr.rain.tot,pwin.rain.tot,psum.rain.tot,pspr.temp,pwin.temp,psum.temp)
rain.mat<-cbind(spr.rain.tot,win.rain.tot,sum.rain.tot)
prain.mat<-cbind(pspr.rain.tot,pwin.rain.tot,psum.rain.tot)
temp.mat<-cbind(spr.temp,win.temp,sum.temp)
gls.time<-gls(resids~YrsOB+YrsSLB+BP5Yrs+rain.mat+prain.mat+temp.mat,cor=corAR1(form=~yr|name),meth="ML")
round(summary(gls.time)$tTable,3)
plot(ACF(gls.time,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(gls.time,maxLag=5,resType='n'),alpha=0.05)

mod.lm<-lm(resids~YrsOB+YrsSLB+BP5Yrs+spr.rain.tot+win.rain.tot+sum.rain.tot+pspr.rain.tot+pwin.rain.tot+psum.rain.tot)
par(mfrow=c(3,3))
termplot(mod.lm,partial=TRUE,se=T,smooth=panel.smooth,span.smth=1/3,lwd.term=1.5,lwd.se=2)

yr[pspr.rain.tot==min(pspr.rain.tot)]##1998
yr[spr.rain.tot==min(spr.rain.tot)]##2005
yr[pwin.rain.tot==max(pwin.rain.tot)]##1998
yr[win.rain.tot==max(win.rain.tot)]##1999
yr[psum.rain.tot==max(psum.rain.tot)]##1998
yr[sum.rain.tot==max(sum.rain.tot)]##1999
##1998 is exerting a lot of leverage on the previous year variable regressions
mod.lm<-lm(resids~YrsOB+YrsSLB+BP5Yrs+spr.rain.tot+win.rain.tot+sum.rain.tot+pspr.rain.tot+pwin.rain.tot+psum.rain.tot,subset=yr!=1998)
summary(mod.lm)
par(mfrow=c(3,3))
termplot(mod.lm,partial=TRUE,se=T,smooth=panel.smooth,span.smth=1/3,lwd.term=1.5,lwd.se=2)

mod.lm<-lm(resids~YrsOB+YrsSLB+BP5Yrs+spr.temp+win.temp+sum.temp+pspr.temp+pwin.temp+psum.temp)
summary(mod.lm)
par(mfrow=c(3,3))
termplot(mod.lm,partial=TRUE,se=T,smooth=panel.smooth,span.smth=1/3,lwd.term=1.5,lwd.se=2)

mod.lm<-lm(resids~dist.m+rain.tot*temp.avg)
anova(mod.lm)
##to simply interpretation just using total rain and avg temp
gls.name<-gls(rich~name,cor=corAR1(form=~yr|name),meth="ML")
resids<-residuals(gls.name,type='n')
gls.time<-gls(resids~YrsOB+YrsSLB+BP5Yrs+rain.tot*temp.avg,cor=corAR1(form=~yr|name),meth="ML")
round(summary(gls.time)$tTable,3)
plot(ACF(gls.time,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(gls.time,maxLag=5,resType='n'),alpha=0.05)
plot(gls.time)
anova(gls.time,gls.time2)
##for varpart
clim.m<-cbind(rain.tot,temp.avg,rain.tot*temp.avg)
lm.time<-lm(resids~dist.m+clim.m)
anova(lm.time)
################################################
################################################
##single mixed model with both within and between components
logCa.avg<-rep(tapply(logCa,name,mean),each=11)
cyr<-scale(yr,scale=F)
cyr2<-cyr^2
cyr3<-cyr^3
mod.lme<-lme(rich~cyr+logCa.avg+YrsOB+YrsSLB+BP5Yrs+rain.tot*temp.avg,random=~1|name,cor=corAR1(form=~yr|name),meth="ML")
mod2.lme<-lme(rich~cyr+cyr2+logCa.avg+YrsOB+YrsSLB+BP5Yrs+rain.tot*temp.avg,random=~1|name,cor=corAR1(form=~yr|name),meth="ML")
mod3.lme<-update(mod2.lme,.~.+cyr3)
plot(yr,resid(mod.lme,type='n'))
lines(lowess(yr,resid(mod.lme,type='n'),f=(1/3)))
lines(lowess(yr,resid(mod2.lme,type='n'),f=(1/3)),col='blue')
lines(lowess(yr,resid(mod3.lme,type='n'),f=(1/3)),col='red')
##shows no need for yr3 term
##mod3.lme is better but not by much and examination of resids suggests not worth it
anova(mod.lme,mod2.lme,mod3.lme)
summary(mod2.lme)
plot(ACF(mod2.lme,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(mod2.lme,maxLag=5,resType='n'),alpha=0.05)
plot(mod.lme)
mod.lme.var<-update(mod2.lme,weights=varPower())
anova(mod2.lme,mod.lme.var)
##heteroscatistic is good to consider
mod.lme<-mod.lme.var
##diagnostic plots
##assumption 1
##Are the within group errors centered around zero?
plot(mod.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
plot(predict(mod.lme),resid(mod.lme,type='p'))
lines(lowess(predict(mod.lme),resid(mod.lme,type='p')),lwd=2,col='red')
abline(h=0)
##look at same assumption on a group basis
plot(mod.lme,name~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod.lme,resid(.,type='p')~fitted(.)|name,id=0.05,adj=-0.3)
##examine difference in response and fitted 
plot(mod.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qq.plot(resid(mod.lme))
qqnorm(mod.lme,~resid(.)|name)
##assumption 2
##qqnorm() on ranef(.) for checking marginal normality and identifyin outliers
##pairs() on ranef(.) for id outliers and cheing assump of homogenity of ranef covariance matrix
qqnorm(mod.lme,~ranef(.,level=1),id=0.10,cex=.7)
qq.plot(ranef(mod.lme,level=1))
##diagnostic plots indicate
summary(mod.lme)
intervals(mod.lme)
pairs(formula(mod.lme),panel=panel.smooth,lwd=2)

###same analysis but on standardized expl vars###
logCa.avg<-rep(tapply(logCa,name,mean),each=11)
syr<-scale(yr)
syr2<-syr^2
slogca<-scale(logCa.avg)
mod.lme<-lme(rich~syr+syr2+slogca+syob+syslb+sbp5+srain*stemp,random=~1|name,cor=corAR1(form=~yr|name),meth="ML")
summary(mod.lme)
plot(ACF(mod.lme,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(mod.lme,maxLag=5,resType='n'),alpha=0.05)
plot(mod.lme)
mod.lme.var<-update(mod.lme,weights=varPower())
anova(mod.lme,mod.lme.var)
##heteroscatistic is good to consider
mod.lme<-mod.lme.var
##diagnostic plots
##assumption 1
##Are the within group errors centered around zero?
plot(mod.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
plot(predict(mod.lme),resid(mod.lme,type='p'))
lines(lowess(predict(mod.lme),resid(mod.lme,type='p')),lwd=2,col='red')
abline(h=0)
##look at same assumption on a group basis
plot(mod.lme,name~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod.lme,resid(.,type='p')~fitted(.)|name,id=0.05,adj=-0.3)
##examine difference in response and fitted 
plot(mod.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qq.plot(resid(mod.lme))
qqnorm(mod.lme,~resid(.)|name)
##assumption 2
##qqnorm() on ranef(.) for checking marginal normality and identifyin outliers
##pairs() on ranef(.) for id outliers and cheing assump of homogenity of ranef covariance matrix
qqnorm(mod.lme,~ranef(.,level=1),id=0.10,cex=.7)
qq.plot(ranef(mod.lme,level=1))
##diagnostic plots indicate
summary(mod.lme)
intervals(mod.lme)
pairs(formula(mod.lme),panel=panel.smooth,lwd=2)


##########END OF OFFICIAL ANALYSIS######################
########################################################
########################################################
########################################################
##Principle components Analysis attempts
##b/c deciding on different groups of variables is a pain
##summarize each with one PC1
dist.m<-cbind(YrsOB,BP5Yrs,YrsSLB,bison,spring,summer,fall)
rownames(dist.m)<-name
princomp(dist.m)
dist.pc<-prcomp(dist.m,scale=T)
plot(dist.pc)
biplot(dist.pc,cex=.5)
summary(dist.pc)
dist.pc1<-dist.pc$x[,1]


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

par(mfrow=c(1,3))
biplot(site.pc,cex=.5)
biplot(dist.pc,cex=.5)
biplot(clim.pc,cex=.5)


summary(mod.lm<-lm(rich~site.pc1+dist.pc1+clim.pc1))
#Multiple R-squared: 0.2944,     Adjusted R-squared: 0.2846 
##to get sweet partial plots
par(mfrow=c(1,3))
termplot(mod.lm,partial=TRUE,se=T,smooth=panel.smooth,span.smth=2/3,lwd.term=2,lwd.se=1.25)
##check output of termplot
##calc exp val with only first 3 coef then diffs will be resids
resids<-rich-cbind(rep(1,220),site.pc1,dist.pc1)%*%coef(mod.lm)[1:3]
plot(clim.pc1,resids)
abline(lm(resids~clim.pc1))
##create mod.gls just for anova purposes down the road
mod.gls<-gls(rich~site.pc1+dist.pc1+clim.pc1,method="ML")
summary(mod.gls)
anova(mod.gls)

yr.m<-cbind(syr,syr2,syr3)
pca.ob<-groupedData(rich~site.pc1+dist.pc1+clim.pc1|name)
mod.ls<-lmList(pca.ob)
plot(intervals(mod.ls))
##definately a strong plot to plot random influ on intercept
##possibily also pc1
mod.lme.int<-lme(rich~site.pc1+dist.pc1+clim.pc1,random=~1|name,data=pca.ob,method="ML")
anova(mod.lme.int,mod.gls)
mod.lme<-lme(rich~1,random=~1|name,data=pca.ob,method="ML")
anova(mod.lme,mod.lme.int)
mod.lme.sit<-update(mod.lme.int,random=~site.pc1|name)
anova(mod.lme.int,mod.lme.sit)
##additional pca.site var not necess as random effect
summary(mod.lme.int)
        (Intercept) Residual
StdDev:    6.973806 8.375835
Fixed effects: rich ~ site.pc1 + dist.pc1 + clim.pc1 
               Value Std.Error  DF  t-value p-value
(Intercept) 76.23636 1.6737742 197 45.54758  0.0000
site.pc1     2.26408 0.6181275 197  3.66280  0.0003
dist.pc1    -1.21348 0.5433549 197 -2.23331  0.0267
clim.pc1     1.50196 0.4268197 197  3.51896  0.0005
mod.lme<-mod.lme.int
##check temporal autocorr
plot(ACF(mod.lme,maxLag=5),alpha=0.05)
varo<-Variogram(mod.lme,maxDist=5,form=~yr|name)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')
mod.lme.acf<-update(mod.lme,correlation=corARMA(p=1,form=~yr|name))
mod.lme.acf2<-update(mod.lme,correlation=corARMA(p=2,form=~yr|name))
mod.lme.acf3<-update(mod.lme,correlation=corARMA(p=3,form=~yr|name))
anova(mod.lme,mod.lme.acf,mod.lme.acf2,mod.lme.acf3)
plot(ACF(mod.lme.acf2,maxLag=5),alpha=.05)
##this seems very strange that the autocorr is increasing with extra terms
##I think this is due to the fact that I did not include any year effects 
##therefore these autoregressive terms are picking up on the temporal trend
##better to not keep any autocorrelation values
##diagnostic plots
##assumption 1
##Are the within group errors centered around zero?
plot(mod.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
mod.lme.var<-update(mod.lme,weights=varPower())
anova(mod.lme.var,mod.lme)
##heteroscadistic term is good
mod.lme<-mod.lme.var
plot(mod.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
plot(predict(mod.lme),resid(mod.lme,type='p'))
lines(lowess(predict(mod.lme),resid(mod.lme,type='p')),lwd=2,col='red')
abline(h=0)
##look at same assumption on a group basis
plot(mod.lme,name~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod.lme,resid(.,type='p')~fitted(.)|name,id=0.05,adj=-0.3)
##examine difference in response and fitted 
plot(mod.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qq.plot(resid(mod.lme),label=name)
qqnorm(mod.lme,~resid(.)|name)
##assumption 2
##qqnorm() on ranef(.) for checking marginal normality and identifyin outliers
##pairs() on ranef(.) for id outliers and cheing assump of homogenity of ranef covariance matrix
qqnorm(mod.lme,~ranef(.,level=1),id=0.10,cex=.7)
qq.plot(ranef(mod.lme,level=1))
#pairs(mod3.lme,~ranef(.,level=1)|name)#inthis case does not make sense
##for multilevel models must assess random effects are each level of grouping
##diagnostic plots indicate
summary(mod.lme)
##graphics!!
plot(augPred(mod.lme,~yr,level=0:1),grid=T)
plot(augPred(mod.lme,~site.pc1,level=0:1),grid=T)
plot(augPred(mod.lme,~dist.pc1,level=0:1),grid=T)
plot(augPred(mod.lme,~clim.pc1,level=0:1),grid=T)

test.lme<-lme(rich~clim.pc1+dist.pc1,random=~site.pc1|name,data=pca.ob)

##partial regression coeff plot with ranef
par(mfrow=c(1,3))
uni.plot<-rownames(ranef(mod.lme)
for(i in 1:20){
 resids<-rich-cbind(rep(1,220),dist.pc1,clim.pc1)%*%fixef(mod.lme)[-2]
beta.cl<-coef(lm(resids~site.pc1))
plot(dist.pc1,resids,col=1,pch=19)
lines(c(min(site.pc1),max(site.pc1)),c(min(site.pc1)*beta.cl[2]+beta.cl[1],max(site.pc1)*beta.cl[2]+beta.cl[1]),type='l',col='red',lwd=2)
resids<-rich-cbind(rep(1,220),site.pc1,dist.pc1)%*%fixef(mod.lme)[-3]
beta.cl<-coef(lm(resids~dist.pc1))
plot(dist.pc1,resids,col=1,pch=19)
lines(c(min(dist.pc1),max(dist.pc1)),c(min(dist.pc1)*beta.cl[2]+beta.cl[1],max(dist.pc1)*beta.cl[2]+beta.cl[1]),type='l',col='red',lwd=2)
resids<-rich-cbind(rep(1,220),site.pc1,dist.pc1)%*%fixef(mod.lme)[-4]
beta.cl<-coef(lm(resids~clim.pc1))
plot(clim.pc1,resids,col=1,pch=19)
lines(c(min(clim.pc1),max(clim.pc1)),c(min(clim.pc1)*beta.cl[2]+beta.cl[1],max(clim.pc1)*beta.cl[2]+beta.cl[1]),type='l',col='red',lwd=2)
#############
##predicted lines
par(mfrow=c(1,3))
plot(site.pc1,rich,pch=19)
lines(c(min(site.pc1),max(site.pc1)),c(min(site.pc1)*fixef(mod.lme)[2]+fixef(mod.lme)[1],max(site.pc1)*fixef(mod.lme)[2]+fixef(mod.lme)[1]),col='red',lwd=2)
#lines(lowess(site.pc1,rich),col='blue',lwd=2)
##just for plots
cls<-rainbow(40)
uni.name<-rownames(ranef(mod.lme))
for(i in 1:20){
 lines(c(min(site.pc1[name==uni.name[i]]),max(site.pc1[name==uni.name[i]])),
       c(min(site.pc1[name==uni.name[i]])*fixef(mod.lme)[2]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1],
         max(site.pc1[name==uni.name[i]])*fixef(mod.lme)[2]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1]),
        type='l',col='grey',lwd=2)
}
for(i in 1:20){
 for(j in 1998:2008){
  xcord<-site.pc1[name==uni.name[i]&yr==j]
  ycord<-rich[name==uni.name[i]&yr==j]
  ynew<- xcord*fixef(mod.lme)[2]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1]
  arrows(xcord,ycord,xcord,ynew,length=0,col='red')
}}
legend("bottomright",c('fixed effect','plot ranef'),col=c("red","grey"),lty=1,lwd=2)

plot(dist.pc1,rich,pch=19)
lines(c(min(dist.pc1),max(dist.pc1)),c(min(dist.pc1)*fixef(mod.lme)[3]+fixef(mod.lme)[1],max(dist.pc1)*fixef(mod.lme)[3]+fixef(mod.lme)[1]),col='red',lwd=2)
#lines(lowess(dist.pc1,rich),col='blue',lwd=2)
##just for plots
uni.name<-rownames(ranef(mod.lme))
for(i in 1:20){
 lines(c(min(dist.pc1[name==uni.name[i]]),max(dist.pc1[name==uni.name[i]])),
       c(min(dist.pc1[name==uni.name[i]])*fixef(mod.lme)[3]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1],
         max(dist.pc1[name==uni.name[i]])*fixef(mod.lme)[3]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1]),
        type='l',col='grey',lwd=2)
}
for(i in 1:20){
 for(j in 1998:2008){
  xcord<-dist.pc1[name==uni.name[i]&yr==j]
  ycord<-rich[name==uni.name[i]&yr==j]
  ynew<- xcord*fixef(mod.lme)[3]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1]
  arrows(xcord,ycord,xcord,ynew,length=0,col='red')
}}
legend("bottomright",c('fixed effect','plot ranef'),col=c("red","grey"),lty=1,lwd=2)

plot(clim.pc1,rich,pch=19)
lines(c(min(clim.pc1),max(clim.pc1)),c(min(clim.pc1)*fixef(mod.lme)[4]+fixef(mod.lme)[1],max(clim.pc1)*fixef(mod.lme)[4]+fixef(mod.lme)[1]),col='red',lwd=2)
#lines(lowess(clim.pc1,rich),col='blue',lwd=2)
##just for plots
uni.name<-rownames(ranef(mod.lme))
for(i in 1:20){
 lines(c(min(clim.pc1[name==uni.name[i]]),max(clim.pc1[name==uni.name[i]])),
       c(min(clim.pc1[name==uni.name[i]])*fixef(mod.lme)[4]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1],
         max(clim.pc1[name==uni.name[i]])*fixef(mod.lme)[4]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1]),
        type='l',col='grey',lwd=2)
}
for(i in 1:20){
 for(j in 1998:2008){
  xcord<-clim.pc1[name==uni.name[i]&yr==j]
  ycord<-rich[name==uni.name[i]&yr==j]
  ynew<- xcord*fixef(mod.lme)[4]+fixef(mod.lme)[1]+ranef(mod.lme,lev=1)[i,1]
  arrows(xcord,ycord,xcord,ynew,length=0,col='red')
}}
legend("bottomright",c('fixed effect','plot ranef'),col=c("red","grey"),lty=1,lwd=2)

#for(i in 1:20){
# lines(lowess(c(clim.pc1[name==uni.name[i]],clim.pc1[name==uni.name[i]]),
#       c(rich[name==uni.name[i]],rich[name==uni.name[i]])),
#        type='l',col='skyblue',lwd=2)
#}
#legend("bottomright",c('fixed effect','plot ranef','lowess'),col=c("red","grey",'skyblue'),lty=1,lwd=2)

#################################################################
#################################################################
#################################################################

###interested in year to year (not long trend) and site to site var
dist.m<-cbind(YrsOB,BP5Yrs,YrsSLB,bison)
mod.yr<-lm(rich~syr+syr2+syr3)
resids<-resid(mod.yr)
dist.mod<-lm(resids~dist.m)
summary(dist.mod)
##drop YrsOB and grazer from reg
dist.m<-cbind(BP5Yrs,YrsSLB)
seas.mod<-lm(resids~sspr.rain+swin.rain+ssum.rain+sspr.temp+swin.temp+ssum.temp)
summary(seas.mod)
pairs(resids~sspr.rain+swin.rain+ssum.rain+sspr.temp+swin.temp+ssum.temp,panel=panel.smooth,lwd=2)
##keep swin.temp and ssum.rain
sea.m<-cbind(swin.temp,ssum.rain)
colnames(sea.m)<-c('win.temp','sum.rain')
pHavgs<-scale(tapply(pH,plot,mean))
Caavgs<-scale(tapply(logCa,plot,mean))
Feavgs<-scale(tapply(logFe,plot,mean))
ph<-rep(pHavgs,each=11)
lCa<-rep(Caavgs,each=11)
lFe<-rep(Feavgs,each=11)
snorth<-scale(northness)
seast<-scale(eastness)
site.mod<-lm(resids~lCa+ph+lFe+snorth+seast)
pairs(rich~lCa+northness+YrsOB+BP5Yrs+YrsSLB,panel=panel.smooth,lwd=2)
pairs(resids~lCa+northness+YrsOB+BP5Yrs+YrsSLB,panel=panel.smooth,lwd=2)
summary(site.mod)
##keep lCa and northness (although eastness outperforms this is just cause corr with Ca)
site.m<-cbind(lCa,snorth)


site.var<-groupedData(resids~lCa+snorth+syslb+sbp5|name)

mod.ls<-lmList(site.var)
plot(mod.ls)
plot(intervals(mod.ls))
mod.gls<-gls(resids~lCa+snorth+syslb+sbp5,method="ML")
summary(mod.gls)
mod.lme<-lme(resids~lCa+snorth+dist.m,random=~1|name,method="ML",data=site.var)
anova(mod.gls,mod.lme)
        Model df      AIC      BIC    logLik   Test  L.Ratio p-value
mod.gls     1  6 1640.487 1660.849 -814.2437                        
mod.lme     2  7 1576.760 1600.516 -781.3802 1 vs 2 65.72717  <.0001
##to test for the impor of a rand eff,consider impor of being close to boundary
.5*0 + .5*(1-pchisq(2*(logLik(mod.lme)-logLik(mod.gls)),1))
##pretty much zero so yes its better
mod2.lme<-update(mod.lme,random=~sbp5|name)
anova(mod.lme,mod2.lme)
##I checked each fixed effect as a corresponding random effect none looked good
plot(mod.lme)
mod2.lme<-update(mod.lme,weights=varPower())
anova(mod.lme,mod2.lme)
mod2.lme<-update(mod.lme,weights=varIdent(form=~1|name))
anova(mod.lme,mod2.lme)
plot(ACF(mod.lme,maxLag=5),alpha=0.01)
##doesn't look like a like of autocorr
mod2.lme<-update(mod.lme,correlation=corAR1(form=~yr|name))
intervals(mod2.lme)
plot(ACF(mod2.lme,maxLag=5),alpha=0.01)
anova(mod.lme,mod2.lme)
##model diagonsis conti
##diagnostic plots
##assumption 1
##Are the within group errors centered around zero?
plot(mod.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
plot(predict(mod.lme),resid(mod.lme,type='p'))
lines(lowess(predict(mod.lme),resid(mod.lme,type='p')),lwd=2,col='red')
abline(h=0)
##look at same assumption on a group basis
plot(mod.lme,name~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod.lme,resid(.,type='p')~fitted(.)|name,id=0.05,adj=-0.3)
##examine difference in response and fitted 
plot(mod.lme,resids~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qq.plot(resid(mod.lme),label=name)
qqnorm(mod.lme,~resid(.)|name)
##assumption 2
##qqnorm() on ranef(.) for checking marginal normality and identifyin outliers
##pairs() on ranef(.) for id outliers and cheing assump of homogenity of ranef covariance matrix
qqnorm(mod.lme,~ranef(.,level=1),id=0.10,cex=.7)
qq.plot(ranef(mod.lme,level=1))
#pairs(mod3.lme,~ranef(.,level=1)|name)#inthis case does not make sense
##for multilevel models must assess random effects are each level of grouping
##diagnostic plots indicate
###heteroscadisity
###residuals with strange irregular behavior
###randef effects not normally distrubited



#################################################
##################################################
##Try out lme4 for within and between plot var
##See DAAG chp10 for example (which is outdated and does not quite work in new version of package)
library(lme4)
mod1.lme<-lmer(rich~(1|name),REML=T)
Random effects:
 Groups   Name        Variance Std.Dev.
 name     (Intercept) 88.124   9.3874  
 Residual             79.798   8.9330 
Fixed effects:
            Estimate Std. Error t value
(Intercept)   76.236      1.617   47.14
##note that they are similar as calculated with aov
##use Markov Chain Monte Carlo approach to construct conf inter
plot.samp<-as.data.frame(mcmcsamp(mod1.lme,n=1000))
test<-as.matrix(mcmcsamp(mod1.lme,n=1000))
names(plot.samp)
[1] "(Intercept)" "ST1"         "sigma"  
CI95<-apply(plot.samp,2,function(x)quantile(x,prob=c(.025,.975)))
# 95% confidence limits of std.dev
CI95
      (Intercept)       ST1     sigma
2.5%     73.06225 0.5284636  8.527149
97.5%    79.34705 0.9534200 10.381144
##intercept is the overall mean
##not sure what ST1 is this method seems to be returning wacky info
##sigma is the residual std
