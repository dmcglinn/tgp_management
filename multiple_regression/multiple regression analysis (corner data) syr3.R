endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
library(nlme)
library(car)
names(endat)
 
##first examine multicollinarity in expl vars
season<-as.factor(ifelse(endat$spring==1,'spring','other'))
grazer<-as.factor(ifelse(endat$bison==1,'bison','cattle'))
mang.cat<-rep(NA,dim(endat)[1])
for(i in 1:length(mang.cat)){
 if(endat$OrBis[i]==1)
  mang.cat[i]<-'OBis'
 if(endat$ChBis[i]==1)
  mang.cat[i]<-'CBis'
 if(endat$OrBis[i]+endat$ChBis[i]==0)
  mang.cat[i]<-'OCat'
}
mang.cat<-as.factor(mang.cat)

site<-as.matrix(endat[,8:27])
dist<-as.matrix(endat[,43:45])
topo<-as.matrix(endat[,50:52])
mon.clim<-as.matrix(endat[,72:95])
sea.clim<-as.matrix(endat[,96:104])
avg.clim<-as.matrix(endat[,105:107])



attach(endat)

##create standardized variables for analysis##

##First explore climate effects on richness and z

yr.std<-scale(yr)
dist.std<-scale(dist)
topo.std<-scale(topo)
mon.clim.std<-scale(mon.clim)
sea.clim.std<-scale(sea.clim)

coef.var<-function(x){ sd(x)/mean(x)}
rain.cv<-apply(endat[,74:83],1,coef.var)
temp.cv<-apply(endat[,84:95],1,coef.var)

##begin with the full model of all climate vars ##
##fit residuals of site effects
##calc avg across 4 corners

plot<-as.factor(plot)
lv<-as.factor(lv)
cor<-as.factor(cor)

site.resids<-residuals(lm(rich1~plot+cor))
plot(yr,site.resids)

library(lattice)
splom(~cbind(site.resids,sea.clim.std)) ##scatterplot matrix
splom(~cbind(rich1,sea.clim.std))


sea.clim.mod<-lm(site.resids~-1+sea.clim.std)
summary(sea.clim.mod)

colnames(sea.clim.std)
[1] "sum.rain.tot" "win.rain.tot" "spr.rain.tot" "sum.rain.avg" "win.rain.avg"
[6] "spr.rain.avg" "sum.temp"     "win.temp"     "spr.temp"    

sea.clim.mod<-lm(site.resids~sea.clim.std[,1]*sea.clim.std[,2]*sea.clim.std[,3]*
 sea.clim.std[,7]*sea.clim.std[,8]*sea.clim.std[,9])

sea.clim.mod<-lm(site.resids~scale(sum.rain.tot)*scale(win.rain.tot)*scale(spr.rain.tot)*scale(sum.temp)*scale(win.temp)*scale(spr.temp))
sea.int.step<-step(sea.clim.mod)
sea.clim.mod<-lm(site.resids~scale(sum.rain.tot)+scale(win.rain.tot)+scale(spr.rain.tot)+scale(sum.temp)+scale(win.temp)+scale(spr.temp))
sea.add.step<-step(sea.clim.mod)
summary(sea.add.step)
##
scale(sum.rain.tot) -4.386e+00  1.570e-01  -27.936  < 2e-16 ***
scale(win.rain.tot)  1.076e+00  1.206e-01    8.926  < 2e-16 ***
scale(spr.rain.tot)  2.167e+00  1.292e-01   16.766  < 2e-16 ***
scale(sum.temp)     -1.405e+00  1.623e-01   -8.655  < 2e-16 ***
scale(win.temp)     -2.366e+00  1.442e-01  -16.412  < 2e-16 ***
scale(spr.temp)      5.834e-01  1.511e-01    3.862 0.000114 ***
Multiple R-squared: 0.2394,     Adjusted R-squared: 0.2384 
###
sea.clim.mod<-lm(site.resids~scale(sum.rain.tot)+scale(win.rain.tot)+scale(spr.rain.tot)+scale(sum.temp)+scale(win.temp)+scale(spr.temp)+
               scale(psum.rain.tot)+scale(pwin.rain.tot)+scale(pspr.rain.tot)+scale(psum.temp)+scale(pwin.temp)+scale(pspr.temp))
sea.lag1.step<-step(sea.clim.mod)
summary(sea.lag1.step)
##
scale(sum.rain.tot)  -4.129e+00  1.492e-01  -27.675  < 2e-16 ***
scale(win.rain.tot)  -3.591e-01  1.391e-01   -2.582  0.00986 ** 
scale(spr.rain.tot)   3.641e+00  3.143e-01   11.586  < 2e-16 ***
scale(sum.temp)      -4.134e+00  3.119e-01  -13.254  < 2e-16 ***
scale(win.temp)      -1.358e+00  1.864e-01   -7.287 3.75e-13 ***
scale(spr.temp)      -1.107e+00  2.236e-01   -4.951 7.66e-07 ***
scale(pspr.rain.tot)  3.044e+00  1.278e-01   23.819  < 2e-16 ***
scale(psum.temp)     -3.693e+00  3.225e-01  -11.451  < 2e-16 ***
Multiple R-squared: 0.4026,     Adjusted R-squared: 0.4016 
##
clim.mod<-lm(site.resids~rain.tot*temp.avg*prain.tot*ptemp.avg)
clim.lag1.step<-step(clim.mod)
summary(clim.lag1.step)
##
rain.tot            -1.757e+01  4.963e-01 -35.392  < 2e-16 ***
temp.avg            -1.279e+02  8.096e+00 -15.799  < 2e-16 ***
prain.tot            7.520e-01  2.182e-01   3.446 0.000574 ***
ptemp.avg           -9.178e+01  6.824e+00 -13.449  < 2e-16 ***
rain.tot:temp.avg    1.068e+00  2.652e-02  40.260  < 2e-16 ***
rain.tot:prain.tot  -5.228e-02  1.201e-03 -43.539  < 2e-16 ***
rain.tot:ptemp.avg   4.644e-01  1.874e-02  24.774  < 2e-16 ***
temp.avg:ptemp.avg   2.126e+00  4.578e-01   4.644 3.51e-06 ***
prain.tot:ptemp.avg  2.576e-01  1.578e-02  16.326  < 2e-16 ***
Multiple R-squared: 0.4026,     Adjusted R-squared: 0.4014
##
##from Adler and Levine 2007
mod1<-lm(site.resids~rain.tot)
mod2<-lm(site.resids~prain.tot)
mod3<-lm(site.resids~I(prain.tot^-1))
mod4<-lm(site.resids~rain.tot+I(prain.tot^-1))
mod5<-lm(site.resids~rain.tot*I(prain.tot^-1))##does appear to be best model
##all mods explain very little of the variance, low r2s
anova(mod1,mod5)
summary(mod5)

##drought indices
mod1<-lm(site.resids~PDSIavg)
mod2<-update(mod1,.~.+I(PDSIavg^2))
plot(PDSIavg,site.resids)
abline(mod1)
x<-seq(min(PDSIavg),max(PDSIavg),length.out=100)
lines(x,coef(mod2)[1]+coef(mod2)[2]*x+coef(mod2)[3]*x^2)
anova(mod1,mod2)
mod3<-lm(site.resids~SPI1+SPI12+SPI24)
mod4<-update(mod2,.~.+SPI1+SPI12+SPI24)
anova(mod2,mod3,mod4)
summary(mod4)
mod.full<-lm(site.resids~scale(PDSIavg)*scale(I(PDSIavg^2))*scale(SPI1)*scale(SPI12)*scale(SPI24))
summary(step(mod.full))

library(nlme)
?lme
rm(cor)
rm(lv)
rm(plot)
##first calculat the average of rich1 by each lv for each corner
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)
c1<-as.factor(cor)
p1<-as.factor(plot)
syr<-scale(yr,scale=F)
syr2<-syr^2
syr3<-syr^3
srain<-scale(rain.tot)
stemp<-scale(temp.avg)
psrain<-scale(prain.tot)
pstemp<-scale(ptemp.avg)

+srain*stemp
climgrp<-groupedData(rich~syr+syr2+syr3|name)

mod1.lis<-lmList(climgrp)
plot(intervals(mod1.lis))
plot(mod1.lis)
pairs(mod1.lis,id=0.01)
MSE<-summary(mod1.lis)$RSE/summary(mod1.lis)$df.residual


mod.pl<-lm(rich~name+syr+syr2+syr3)
plot(yr,resid(mod.pl))
lines(lowess(yr,resid(mod.pl)))
mod1.lm<-lm(resid(mod.pl)~srain*stemp)
summary(mod1.lm)
mod2.lm<-update(mod1.lm,.~.-srain:stemp)
summary(mod2.lm)
anova(mod1.lm,mod2.lm)##interaction is needed
plot(mod1.lm,labels.id=name)
plot(srain,resid(mod1.lm))
lines(lowess(srain,resid(mod1.lm)))
plot(stemp,resid(mod1.lm))
lines(lowess(stemp,resid(mod1.lm)))
##maybe a sq2 temp term would be good
mod3.lm<-update(mod1.lm,.~.+I(stemp^2))
anova(mod1.lm,mod3.lm)
##squared term does well
summary(mod3.lm)
plot(mod3.lm,labels.id=name)
plot(srain,resid(mod3.lm))
lines(lowess(srain,resid(mod3.lm)))
plot(stemp,resid(mod3.lm))
lines(lowess(stemp,resid(mod3.lm)))
plot(stemp*srain,mod3.lm$fit)
lines(lowess(stemp*srain,mod3.lm$fit))


##start with mixed model in which each plot has a diff intercept but common slope
mod1.gls<-gls(rich~syr+syr2+syr3+srain*stemp,data=climgrp)
mod1.lme<-lme(rich~syr+syr2+syr3+srain*stemp,random=~1|p1,data=climgrp)
anova(mod1.gls,mod1.lme)
##pays to consider plots as random effects
##consider a model with different incerpts for each plot and Beta for year effect
mod2.lme<-update(mod1.lme,random=~syr|p1) 
anova(mod1.lme,mod2.lme)
##mod2.lme is the superior model with the seperate year effects
##consider dropping uneven means for plots although this is proba really bad idea
mod2b.lme<-update(mod2.lme,random=~-1+syr|p1)
anova(mod2.lme,mod2b.lme)
##a seperate intercept for each plot is definately better
##next consider non common slope and intercept, pos diagonal corr matrix
##a pdDiag assumes independence between random effects (intercept and slope random independent)
##i.e. the correlation term is not needed
mod3.lme<-update(mod2.lme,random=pdDiag(~syr))
anova(mod2.lme,mod3.lme)
anova(mod
##mod3.lme is marginally superior
summary(mod3.lme)
mod2.lme<-mod3.lme
##diagnostic plots
##assumption 1
##Are the within group errors centered around zero?
plot(mod2.lme,abline=0,id=0.05,adj=-.3)
##look at same assumption on a group basis
plot(mod2.lme,p1~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod2.lme,resid(.,type='p')~fitted(.)|p1,id=0.05,adj=-0.3)
##because some plots seem to have slightly different variance
##we can consider a sep variance pred for each plot
mod3.lme<-update(mod2.lme,weights=varIdent(form=~1|p1))
plot(mod3.lme,p1~resid(.),abline=0)##centered around zero but plots seem to have diff variances

##b/c slight heteroscadisticity
mod4.lme<-update(mod2.lme,weights=varPower())
anova(mod2.lme,mod3.lme)
##one variance for all plots seems best
anova(mod2.lme,mod4.lme)
##is power function of var better
plot(mod4.lme,abline=0,id=0.05,adj=-.3)
##look at same assumption on a group basis
plot(mod4.lme,p1~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod4.lme,resid(.,type='p')~fitted(.)|p1,id=0.05,adj=-0.3)

##examine difference in response and fitted 
plot(mod4.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qq.plot(mod4.lme)
qq.plot(mod4.lme,~resid(.)|p1)
##all dianostic plots look reasonably good
random.effects(mod4.lme)
plot(augPred(mod4.lme,~srain,level=0:1),grid=T)
plot(augPred(mod4.lme,~stemp,level=0:1),grid=T)

##compare coef with lmList object
plot(compareFits(coef(mod1.lis),coef(mod4.lme)),mark=fixef(mod4.lme))
##compare pred with lmList object
plot(comparePred(mod1.lis,mod4.lme,primary=~syr),layout=c(5,4))
plot(comparePred(mod1.lis,mod4.lme,primary=~srain,level=0),layout=c(5,4))
plot(comparePred(mod1.lis,mod4.lme,primary=~srain,level=1),layout=c(5,4))
plot(comparePred(mod1.lis,mod4.lme,primary=~srain,level=1),layout=c(5,4))

plot(ACF(mod4.lme,maxLag=5),alpha=0.01)
varo<-Variogram(mod4.lme,maxDist=5,form=~syr)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')
##variogram indicates no autocorrelation

##there does not appear to be a large degree of autocorrelation
mod5.lme<-update(mod4.lme,corr=corARMA(p=1,form=~syr))
mod6.lme<-update(mod4.lme,corr=corAR1())
mod7.lme<-update(mod4.lme,corr=corARMA(p=0,q=2))#moving average 1 lag
anova(mod4.lme,mod5.lme)
anova(mod4.lme,mod7.lme)
anova(mod4.lme,mod6.lme)
mod7.lme<-update(mod2.lme,corr=corLin(form=~syr))
anova(mod2.lme,mod7.lme)
###mod4 still superior

##considering autocorrelation does not seem neccessary
summary(mod4.lme)
Linear mixed-effects model fit by REML
 Data: climgrp 
       AIC      BIC    logLik
  1560.249 1597.224 -769.1247

Random effects:
 Formula: ~syr | p1
 Structure: Diagonal
        (Intercept)      syr   Residual
StdDev:    9.427788 1.136930 0.01333510

Variance function:
 Structure: Power of variance covariate
 Formula: ~fitted(.) 
 Parameter estimates:
   power 
1.430094 
Fixed effects: rich ~ syr + syr2 + syr3 + srain * stemp 
               Value Std.Error  DF  t-value p-value
(Intercept) 76.76430 2.2416650 194 34.24432  0.0000
syr         -0.02142 0.4658798 194 -0.04598  0.9634
syr2        -0.10081 0.0586789 194 -1.71792  0.0874
syr3         0.07128 0.0186006 194  3.83212  0.0002
srain       -1.55049 0.6963887 194 -2.22647  0.0271
stemp       -1.65933 0.4491481 194 -3.69440  0.0003
srain:stemp  3.08628 0.7048169 194  4.37884  0.0000
 Correlation: 
            (Intr) syr    syr2   syr3   srain  stemp 
syr         -0.020                                   
syr2        -0.275  0.049                            
syr3         0.017 -0.774  0.004                     
srain        0.166 -0.133 -0.596  0.148              
stemp        0.014 -0.105 -0.067  0.058  0.199       
srain:stemp -0.118  0.367  0.412 -0.344 -0.702 -0.376

Standardized Within-Group Residuals:
         Min           Q1          Med           Q3          Max 
-2.813731019 -0.622997127  0.005382754  0.653559324  2.843193632 

Number of Observations: 220
Number of Groups: 20 
intervals(mod4.lme)


##attempt new data function of predict, 
##NOTE you need all covariables in equation to get correct predicted values
test<-cbind(p1,syr,stemp,srain)
test<-data.frame(test[order(stemp),])
colnames(test)<-c('p1','syr','stemp','srain')
plot(test$stemp,predict(mod4.lme,newdata=test,level=0),type='l',lwd=2,ylim=c(50,100))
##now its working but not really necessary unless new values will be used
##for example
newtemp<-rep(seq(min(stemp),max(stemp),length.out=11),each=20)
test<-cbind(p1,syr,newtemp,srain)
colnames(test)<-c('p1','syr','stemp','srain')
plot(newtemp,predict(mod4.lme,newdata=test,level=0),type='l',lwd=2,ylim=c(50,100))
##this example is not working correctly!!

##instead of using the newdata option of predict simply pull from fitted values
par(mfrow=c(1,3))
## temp effect
##if you want them connected with a line use the below
#plot(stemp[order(stemp)],predict(mod4.lme,level=0)[order(stemp)],type='l',lwd=2,col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
##b/c this does not make much sense (not seq in time) use
plot(stemp,predict(mod4.lme,l=0),col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(stemp,predict(mod4.lme,l=0)),col='red',lwd=2)
points(stemp,predict(mod4.lme,l=1),col='blue')##random eff
lines(lowess(stemp,predict(mod4.lme,l=1)),col='blue',lwd=2)
points(sort(unique(stemp)),tapply(rich,stemp,mean),type='b',pch=19)
lines(lowess(stemp,rich),lwd=2)
legend("bottomright",c('fixed effect','random effect'),col=c('red','blue'),pch=1)
## rain effect
plot(srain,mod4.lme$fit[,1],col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(srain,mod4.lme$fit[,1]),col='red')
points(srain,mod4.lme$fit[,2],col='blue')##random eff
lines(lowess(srain,mod4.lme$fit[,2]),col='blue')
points(sort(unique(srain)),tapply(rich,srain,mean),type='b',pch=19)
lines(lowess(srain,rich),lwd=2)
legend("bottomleft",c('fixed effect','random effect'),col=c('red','blue'),pch=1)
## interaction
plot(srain*stemp,mod4.lme$fit[,1],col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(srain*stemp,mod4.lme$fit[,1]),col='red')
points(srain*stemp,mod4.lme$fit[,2],col='blue')##random eff
lines(lowess(srain*stemp,mod4.lme$fit[,2]),col='blue')
lines(lowess(srain*stemp,rich),lwd=2)
legend("bottomleft",c('fixed effect','random effect'),col=c('red','blue'),pch=1)

###instead of looking at model output look at observed patts
##look at each covariate against centered richness
plrich<-tapply(rich,p1,mean)
crich<-rich-rep(plrich,11)
par(mfrow=c(1,3))
#plot(syr+2003,crich,xlab="Time",ylab="")
#lines(lowess(syr+2003,crich))
plot(stemp,crich,xlab="Rain standardized",ylab="")
lines(lowess(stemp,crich))
plot(srain,crich,xlab="Temp standardized",ylab="")
lines(lowess(srain,crich))
plot(stemp*srain,crich,xlab="Temp*Rain",ylab="")
lines(lowess(stemp*srain,crich))

##none of the above relationships look good so I will look only at within plot effects
par(mfrow=c(1,3))
pls<-unique(p1)
plot(stemp,crich,xlab="Temp",ylab="",type='n',ylim=c(-25,25))
for(i in 1:20){ #plot
 points(stemp[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(stemp[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}
plot(srain,crich,xlab="srain",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points(srain[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(srain[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}
plot(srain*stemp,crich,xlab="srain*temp",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points((srain*stemp)[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess((srain*stemp)[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}



##create hi low rain values
rain.f<-factor(ifelse(srain<median(srain),"low rain","high rain"))
interaction.plot(x.factor=factor(round(stemp,1)),trace.factor=rain.f,response=predict(mod4.lme,level=1))
#create a hot cold temp vals
temp.f<-factor(ifelse(stemp<median(stemp),"cold","hot"))
interaction.plot(x.factor=factor(round(srain,1)),trace.factor=temp.f,response=predict(mod4.lme,level=1))
interaction.plot(x.factor=factor(round(srain,1)),trace.factor=temp.f,response=predict(mod4.lme,level=1),type='p')
interaction.plot(x.factor=factor(round(srain,1)),trace.factor=factor(round(stemp,1)),response=rich,type='p')
interaction.plot(trace.factor=factor(round(srain,1)),x.factor=factor(round(stemp,1)),response=rich,type='p')
###these plots suggest that the interaction effect is based on a very sparse data
############################
##seasonal ananlysis
##begin by deciding on functional form of seasonal variables
ssum.rain<-scale(sum.rain.tot)
sspr.rain<-scale(spr.rain.tot)
swin.rain<-scale(win.rain.tot)
ssum.temp<-scale(sum.temp)
sspr.temp<-scale(spr.temp)
swin.temp<-scale(win.temp)


mod.pl1<-lm(rich~p1+syr)
mod.pl2<-update(mod.pl1,.~.+I(syr^2))
anova(mod.pl1,mod.pl2)
##quadratic term very important
mod.pl3<-update(mod.pl2,.~.+I(syr^3))
anova(mod.pl2,mod.pl3)
##cubic term better
mod.pl4<-update(mod.pl3,.~.+I(syr^4))
anova(mod.pl3,mod.pl4)
plot(syr,resid(mod.pl2))
lines(lowess(syr,resid(mod.pl2)),lwd=2)
plot(syr,resid(mod.pl3))
lines(lowess(syr,resid(mod.pl3)),lwd=2)
plot(syr,rich)
lines(unique(syr),tapply(predict(mod.pl3),syr,mean),lwd=2)
mod.pl<-mod.pl3

mod1.sea<-lm(resid(mod.pl)~ssum.rain+sspr.rain+swin.rain+ssum.temp+sspr.temp+swin.temp)
summary(mod1.sea)  
ssum.rain   -2.866e+00  6.370e-01    -4.499 1.12e-05 ***
sspr.rain    1.340e+00  5.244e-01     2.555   0.0113 *  
swin.rain    1.970e+00  4.893e-01     4.027 7.87e-05 ***
ssum.temp   -8.282e-01  6.586e-01    -1.257   0.2100    
sspr.temp    1.361e+00  6.130e-01     2.220   0.0275 *  
swin.temp   -2.811e+00  5.849e-01    -4.806 2.91e-06 ***
mod.sea<-lm(resid(mod.pl)~ssum.rain*sspr.rain*swin.rain*ssum.temp*sspr.temp*swin.temp)
summary(step(mod.sea))
###simplify model 1
mod2.sea<-update(mod1.sea,.~.-ssum.temp)
summary(mod2.sea)
anova(mod1.sea,mod2.sea)
##mod2 is superior
summary(mod2.sea)
Coefficients:
              Estimate Std. Error   t value Pr(>|t|)    
(Intercept) -1.724e-15  4.544e-01 -3.79e-15 1.000000    
ssum.rain   -2.396e+00  5.167e-01    -4.638 6.13e-06 ***
sspr.rain    1.020e+00  4.594e-01     2.221 0.027388 *  
swin.rain    1.861e+00  4.821e-01     3.859 0.000151 ***
sspr.temp    1.459e+00  6.088e-01     2.396 0.017456 *  
swin.temp   -3.002e+00  5.656e-01    -5.308 2.76e-07 ***
Residual standard error: 6.74 on 214 degrees of freedom
Multiple R-squared: 0.2017,     Adjusted R-squared: 0.183 
F-statistic: 10.81 on 5 and 214 DF,  p-value: 2.760e-09 
##build lmList model
clim.sea<-groupedData(rich~syr+syr2+syr3+ssum.rain+sspr.rain+swin.rain+swin.temp|p1)
mod1.lis<-lmList(clim.sea)
##runs out of degrees of freedom
plot(intervals(mod1.lis))##looks like common slopes will work except for yr
plot(mod1.lis)
pairs(mod1.lis,id=0.01)
summary(mod1.lis)$RSE
summary(mod1.lis)$df.residual
##build mixed effect model
mod1.lme<-lme(rich~syr+syr2+syr3+ssum.rain+sspr.rain+swin.rain+spr.temp+swin.temp,random=~1|p1,data=clim.sea,method="ML")
mod1b.lme<-update(mod1.lme,.~.-spr.temp)
anova(mod1.lme,mod1b.lme)
##spr.temp term not needed
mod1.lme<-lme(rich~syr+syr2+syr3+ssum.rain+sspr.rain+swin.rain+swin.temp,random=~1|p1,data=clim.sea)
summary(mod1.lme)
mod2.lme<-update(mod1.lme,random=~syr|p1)
anova(mod1.lme,mod2.lme)
##noncommon slopes for year effect superior
summary(mod2.lme)

##diagonstic
##assumption 1
##Are the within group errors centered around zero?
plot(mod2.lme,abline=0,id=0.05,adj=-.3)
##look at same assumption on a group basis
plot(mod2.lme,p1~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod2.lme,resid(.,type='p')~fitted(.)|p1,id=0.05,adj=-0.3)
##because some plots seem to have slightly different variance
##we can consider a sep variance pred for each plot
mod3.lme<-update(mod2.lme,weights=varIdent(form=~1|p1))
plot(mod3.lme,p1~resid(.),abline=0)##centered around zero but plots seem to have diff variances

##b/c slight heteroscadisticity
mod4.lme<-update(mod2.lme,weights=varPower())
anova(mod2.lme,mod3.lme)
##one variance for all plots seems best
anova(mod2.lme,mod4.lme)
##is power function of var better
plot(mod4.lme,abline=0,id=0.05,adj=-.3)
##look at same assumption on a group basis
plot(mod4.lme,p1~resid(.),abline=0)##centered around zero but plots seem to have diff variances
##Is there the same degree of variability in each group
plot(mod4.lme,resid(.,type='p')~fitted(.)|p1,id=0.05,adj=-0.3)

##examine difference in response and fitted 
plot(mod4.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qq.plot(mod4.lme)
qqnorm(mod4.lme,~resid(.)|p1)
##all dianostic plots look reasonably good
random.effects(mod4.lme)
plot(augPred(mod4.lme,~syr,level=0:1),grid=T)
plot(augPred(mod4.lme,~stemp,level=0:1),grid=T)
##not working

##compare coef with lmList object
plot(compareFits(coef(mod1.lis),coef(mod4.lme)),mark=fixef(mod4.lme))
##compare pred with lmList object
plot(comparePred(mod1.lis,mod4.lme,primary=~syr),layout=c(5,4))
plot(comparePred(mod1.lis,mod4.lme,primary=~srain,level=0),layout=c(5,4))
plot(comparePred(mod1.lis,mod4.lme,primary=~srain,level=1),layout=c(5,4))
plot(comparePred(mod1.lis,mod4.lme,primary=~srain,level=1),layout=c(5,4))

plot(ACF(mod4.lme,maxLag=5),alpha=0.01)
varo<-Variogram(mod4.lme,maxDist=5,form=~syr)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')
##variogram indicates no autocorrelation

##there does not appear to be a large degree of autocorrelation
mod5.lme<-update(mod4.lme,correlation=corAR1(form=~syr|p1))
anova(mod4.lme,mod5.lme)
###mod4 still superior

##considering autocorrelation does not seem neccessary
###compare season model to total rain and temp + interaction model
mod1.lme<-lme(rich~syr+syr2+syr3+srain*stemp,random=~syr|p1,corr=corAR1(form=~syr),method="ML")
mod2.lme<-lme(rich~syr+syr2+syr3ssum.rain+sspr.rain+swin.rain+swin.temp,random=~syr|p1,corr=corAR1(form=~syr),method="ML")
##mod1.lme is not quite what we found to be the best model above
##it has been fudged to make the comparison clearer
anova(mod1.lme,mod2.lme)
##mod2.lme with the seasonal breakdown does a much better job
##now lets check the mod1.lme that we actually know was the best
mod1.lme<-lme(rich~syr+srain*stemp,random=~syr|p1,method="ML")
anova(mod1.lme,mod2.lme)
##mod2.lme still the best even with the bigger diff in d.f.s

##Graphical Exploration
par(mfrow=c(1,4))
## sum rain effect
##if you want them connected with a line use the below
#plot(stemp[order(stemp)],predict(mod4.lme,level=0)[order(stemp)],type='l',lwd=2,col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
##b/c this does not make much sense (not seq in time) use
plot(ssum.rain,predict(mod4.lme,l=0),col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(ssum.rain,predict(mod4.lme,l=0)),col='red',lwd=2)
points(ssum.rain,predict(mod4.lme,l=1),col='blue')##random eff
lines(lowess(ssum.rain,predict(mod4.lme,l=1)),col='blue',lwd=2)
points(sort(unique(ssum.rain)),tapply(rich,ssum.rain,mean),type='b',pch=19)
lines(lowess(ssum.rain,rich),lwd=2)
legend("bottomright",c('fixed effect','random effect'),col=c('red','blue'),pch=1)
## spring rain effect
plot(sspr.rain,mod4.lme$fit[,1],col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(sspr.rain,mod4.lme$fit[,1]),col='red')
points(sspr.rain,mod4.lme$fit[,2],col='blue')##random eff
lines(lowess(sspr.rain,mod4.lme$fit[,2]),col='blue')
points(sort(unique(sspr.rain)),tapply(rich,sspr.rain,mean),type='b',pch=19)
lines(lowess(sspr.rain,rich),lwd=2)
legend("bottomleft",c('fixed effect','random effect'),col=c('red','blue'),pch=1)
## winter rain effect
plot(swin.rain,mod4.lme$fit[,1],col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(swin.rain,mod4.lme$fit[,1]),col='red')
points(swin.rain,mod4.lme$fit[,2],col='blue')##random eff
lines(lowess(swin.rain,mod4.lme$fit[,2]),col='blue')
points(sort(unique(swin.rain)),tapply(rich,swin.rain,mean),type='b',pch=19)
lines(lowess(swin.rain,rich),lwd=2)
legend("bottomleft",c('fixed effect','random effect'),col=c('red','blue'),pch=1)
## winter temp effect
plot(swin.temp,mod4.lme$fit[,1],col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(swin.temp,mod4.lme$fit[,1]),col='red')
points(swin.temp,mod4.lme$fit[,2],col='blue')##random eff
lines(lowess(swin.temp,mod4.lme$fit[,2]),col='blue')
points(sort(unique(swin.temp)),tapply(rich,swin.temp,mean),type='b',pch=19)
lines(lowess(swin.temp,rich),lwd=2)
legend("bottomleft",c('fixed effect','random effect'),col=c('red','blue'),pch=1)

###instead of looking at model output look at observed patts
##look at each covariate against centered richness
plrich<-tapply(rich,p1,mean)
crich<-rich-rep(plrich,11)
par(mfrow=c(1,4))
#plot(syr+2003,crich,xlab="Time",ylab="")
#lines(lowess(syr+2003,crich))
plot(ssum.rain,crich,xlab="sum rain",ylab="")
lines(lowess(ssum.rain,crich))
plot(sspr.rain,crich,xlab="spr rain",ylab="")
lines(lowess(sspr.rain,crich))
plot(swin.rain,crich,xlab="win rain",ylab="")
lines(lowess(swin.rain,crich))
plot(swin.temp,crich,xlab="win temp",ylab="")
lines(lowess(swin.temp,crich))

##none of the above relationships look good so I will look only at within plot effects
par(mfrow=c(1,4))
pls<-unique(p1)
plot(ssum.rain,crich,xlab="sum rain",ylab="",type='n',ylim=c(-25,25))
for(i in 1:20){ #plot
 points(ssum.rain[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(ssum.rain[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}
plot(sspr.rain,crich,xlab="spr rain",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points(sspr.rain[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(sspr.rain[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}
plot(swin.rain,crich,xlab="win rain",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points(swin.rain[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(swin.rain[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}
plot(swin.temp,crich,xlab="win temp",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points(swin.temp[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(swin.temp[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}

detach(endatalv1)
#######################################################################
##############disturbance variables##############
##begin looking at 100m2 scale

syob<-scale(YrsOB,scale=T)##b/c a continuous time var
sbp5<-scale(BP5Yrs,scale=FALSE)##
syslb<-scale(YrsSLB,scale=T)##b/c a cont var
grazer<-factor(ifelse(bison==1,'bison','cattle'))
##Mike's approach (plot avg year avg)
yr.f<-factor(yr)
mod1.lm<-lm(rich~p1+yr.f)
resids<-resid(mod1.lm)
mod1.gls<-gls(resids~syob+sbp5+syslb+grazer,method="ML")
summary(mod1.gls)
(Intercept)  -2.150932 0.7676804 -2.801859  0.0055
syob          2.010551 0.6574761  3.057983  0.0025
sbp5         -1.244453 0.4981085 -2.498357  0.0132
syslb        -1.934519 0.6062993 -3.190701  0.0016
grazercattle  4.594224 1.3602285  3.377539  0.0009
#residuals
par(mfrow=c(1,4))
plot(syr,resid(mod1.gls))
lines(lowess(syr,resid(mod1.gls)))
plot(syslb,resid(mod1.gls))
lines(lowess(syslb,resid(mod1.gls)))
plot(sbp5,resid(mod1.gls))
lines(lowess(sbp5,resid(mod1.gls)))
##predicted rel
par(mfrow=c(1,4))
plot(syslb,resid(mod1.lm),ylim=c(-20,20))
abline(coef=c(coef(mod1.gls)[1],coef(mod1.gls)[4]),col='red')
plot(sbp5,resid(mod1.lm),ylim=c(-20,20))
abline(coef=c(coef(mod1.gls)[1],coef(mod1.gls)[3]),col='red')
plot(syob,resid(mod1.lm),ylim=c(-20,20))
abline(coef=c(coef(mod1.gls)[1],coef(mod1.gls)[2]),col='red')
plot(grazer,resid(mod1.lm),ylim=c(-20,20))
abline(coef=c(coef(mod1.gls)[1],coef(mod1.gls)[5]),col='red')
##check autocorr within plots
plot(ACF(mod1.gls,form=~syr|p1,maxLag=5),alpha=0.01)
mod2.gls<-update(mod1.gls,correlation=corAR1(form=~syr|p1))
plot(ACF(mod2.gls,maxLag=5),alpha=0.01)
#no apparent visual change
anova(mod1.gls,mod2.gls)
#but model is sig better
###this model indicates that dist vars can explain
#a sig amt of remaining var after considering plot and yr effects
#but only a tiny amount of the variance
#also the relationships between richness and the explanatory variables is very messy
#the goal now is to consider mixed effect to see if improvement over gls model
mod1.lme<-lme(resids~syob+sbp5+syslb+grazer,random=~1|p1,method="ML")
summary(mod1.lme)
anova(mod1.gls,mod1.lme)
##given that plot means were already factored out the grouping variable does not help at all

##build a stardard lm model to look at functional form
mod.p1<-lm(rich~syr+p1)
mod.p2<-update(mod.p1,.~.+I(syr^2))
mod.p3<-update(mod.p2,.~.+I(syr^3))
anova(mod.p1,mod.p2,mod.p3)
mod.pl<-mod.p3
mod1.lm<-lm(resid(mod.pl)~-1+grazer)
summary(mod1.lm)
##no grazer effect
##fit a stepwise for fun
mod.full<-update(mod1.lm,.~.+grazer*syob*sbp5*syslb)
summary(step(mod.full))
##indicates that bison,syob may be the only important main effects
mod2.lm<-update(mod1.lm,.~.+syob+sbp5+syslb)
summary(mod2.lm)
anova(mod1.lm,mod2.lm)
##mod2 is sig better
mod3.lm<-update(mod2.lm,.~.+syob:sbp5)
summary(mod3.lm)
anova(mod2.lm,mod3.lm)
##mod3.lm sig better, but interaction may be diff to interpt

syr2<-syr^2
##now fit lmList mode
dobject<-groupedData(rich~syr+syr2+syr3+syob+sbp5+syslb|name)
mod1.lis<-lmList(dobject)
plot(intervals(mod1.lis))
plot(mod1.lis)##heteroscadistic
pairs(mod1.lis,id=0)
summary(mod1.lis)$RSE
summary(mod1.lis)$df.residual
##build 4 models
#basic ols
mod.ols<-gls(rich~syr+syr2+syob+grazer+syob+sbp5+syslb,method="ML")
#uncond means mod
mod.unc<-lme(rich~1,random=~1|p1,data=dobject,method="ML")
VarCorr(mod.unc)
p1 = pdLogChol(1) 
            Variance StdDev  
(Intercept) 88.12396 9.387436
Residual    79.79818 8.932983
tau.sq<-as.numeric(VarCorr(mod.unc)[1,1])
sigma.sq<-as.numeric(VarCorr(mod.unc)[2,1] )
tau.sq/(tau.sq+sigma.sq)
[1] 0.5109004 #is the correlation of observations in same group
##between and within plot var almost the same, but w/in slightly greater
#random intercept model
mod1.lme<-lme(rich~syr+syr2+syr3+syob+grazer+syob+sbp5+syslb,random=~1|name,data=dobject,method="ML")
summary(mod1.lme)
##to test for importance of random effects
##creat ols model
mod.OLS<-lm(rich~syr+syr2+syr3+syob+grazer+syob+sbp5+syslb)
##chi sqr test
1-pchisq(2*(logLik(mod1.lme)-logLik(mod.OLS)),1)
##not strictly correct (Verbeke and Molenberghs 2000, p. 64–73
##http://www.unc.edu/courses/2006spring/ecol/145/001/docs/lectures/lecture43.htm
.5*0 + .5*(1-pchisq(2*(logLik(mod1.lme)-logLik(mod.OLS)),1))
##pseudo R2
(as.numeric(VarCorr(mod.unc)[2,1]) - as.numeric(VarCorr(mod1.lme)[2,1]))/ as.numeric(VarCorr(mod.unc)[2,1])
[1] 0.3415035 ## so 34% of var in rich is explained by disturbance variables
anova(mod.unc,mod1.lme)
#random inter and slope model
mod2.lme<-update(mod1.lme,random=~syr|name)
anova(mod1.lme,mod2.lme)
         Model df      AIC      BIC    logLik   Test L.Ratio p-value
mod1.lme     1 10 1573.366 1607.302 -776.6828                       
mod2.lme     2 12 1571.240 1611.964 -773.6201 1 vs 2 6.12541  0.0468
#random slope is better but only slightly and test is flawed
 .5*(1-pchisq(6.12541,1)) + .5*(1-pchisq(6.12541,2))
[1] 0.03004314
#still only marginal support

mod3.lme<-update(mod1.lme,weights=varPower(),control=lmeControl(maxIter=200,msMaxIter=200,niterEM=100,msVerbose=TRUE))
# maxIter=200 increases the number of iterations for the lme optimization from 50 to 200.
# msMaxIter=200 increases the maximum number of iterations for the nlm optimization step from 50 to 200.
# niterEM=100 increases the number of iterations used by the EM algorithm in finding initial estimates for the random effects covariance matrix from 25 to 100.
# msVerbose=TRUE
## no convergence if attempting to update mod2.lme
anova(mod1.lme,mod3.lme)
##considering heteroscadisticity is helpful
summary(mod3.lme)
Linear mixed-effects model fit by maximum likelihood
 Data: dobject 
       AIC      BIC    logLik
  1564.681 1602.011 -771.3406

Random effects:
 Formula: ~1 | name
        (Intercept)   Residual
StdDev:    8.741482 0.01856811

Variance function:
 Structure: Power of variance covariate
 Formula: ~fitted(.) 
 Parameter estimates:
   power 
1.374707 
Fixed effects: rich ~ syr + syr2 + syr3 + syob + grazer + syob + sbp5 + syslb 
                Value Std.Error  DF  t-value p-value
(Intercept)  76.25414 2.2792505 193 33.45580  0.0000
syr          -1.23035 0.4139506 193 -2.97222  0.0033
syr2         -0.25192 0.0535585 193 -4.70360  0.0000
syr3          0.08453 0.0197241 193  4.28555  0.0000
syob          5.42744 1.2323749 193  4.40405  0.0000
grazercattle  4.98371 1.9293798 193  2.58306  0.0105
sbp5         -2.91665 0.7083873 193 -4.11731  0.0001
syslb        -2.51818 0.7602566 193 -3.31228  0.0011

##simple noncommon intercpet for plots and varPower superior

##diagnostic plots
##assumption 1
##Are the within group errors centered around zero?
plot(mod3.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
##look at same assumption on a group basis
plot(mod3.lme,name~resid(.),abline=0)##centered around zero and roughly all same var
##Is there the same degree of variability in each group
plot(mod3.lme,resid(.,type='p')~fitted(.)|name,adj=-0.3)
##because some plots seem to have slightly different variance
##we can consider a sep variance pred for each plot
mod4.lme<-update(mod2.lme,weights=varIdent(form=~1|name))
anova(mod3.lme,mod4.lme)
##examine difference in response and fitted 
plot(mod3.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qq.plot(resid(mod1.lme),label=name)
qqnorm(mod3.lme,~resid(.)|name)
##assumption 2
##qqnorm() on ranef(.) for checking marginal normality and identifyin outliers
##pairs() on ranef(.) for id outliers and cheing assump of homogenity of ranef covariance matrix
qqnorm(mod3.lme,~ranef(.,level=1),id=0.10,cex=.7)
qq.plot(ranef(mod3.lme,level=1))
#pairs(mod3.lme,~ranef(.,level=1)|name)#inthis case does not make sense
##for multilevel models must assess random effects are each level of grouping
##all dianostic plots look reasonably good

random.effects(mod3.lme)

plot(ACF(mod3.lme,maxLag=5),alpha=0.01)
varo<-Variogram(mod2.lme,maxDist=5,form=~syr)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')
##variogram indicates no autocorrelation

##there does not appear to be a large degree of autocorrelation
mod4.lme<-update(mod2.lme,corr=corARMA(p=1,form=~syr|name))
plot(ACF(mod4.lme,maxLag=5),alpha=0.01)
anova(mod3.lme,mod4.lme)
##autoregressive not necessary

###exploring disturbance variables graphically
##year effect
plot(syr,mod3.lme$fit[,1],col='red',ylim=c(50,100))
points(syr,mod3.lme$fit[,2],col='blue')
points(-5:5,tapply(mod3.lme$fit[,1],syr,mean),type='l')
points(-5:5,tapply(rich,syr,mean),type='o')
legend("bottomright",c('fixed effect','random effect'),col=c('red','blue'),pch=1)

##why do the random predicts and fixed predicted have the same mean?

par(mfrow=c(1,3))
##years of bison effect
plot(syob,mod3.lme$fit[,1],col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(syob,mod3.lme$fit[,1]),col='red')
points(syob,mod3.lme$fit[,2],col='blue')##random eff
lines(lowess(syob,mod3.lme$fit[,2]),col='blue')
lines(lowess(syob,rich),lwd=2)
legend("bottomright",c('fixed effect','random effect'),col=c('red','blue'),pch=1)

##sbp5
plot(sbp5,mod3.lme$fit[,1],col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(sbp5,mod3.lme$fit[,1]),col='red')
points(sbp5,mod3.lme$fit[,2],col='blue')##random eff
lines(lowess(sbp5,mod3.lme$fit[,2]),col='blue')
points(sort(unique(sbp5)),tapply(rich,BP5Yrs,mean),type='b',pch=19)
#lines(lowess(sbp5,rich),lwd=2)
legend("bottomleft",c('fixed effect','random effect'),col=c('red','blue'),pch=1)

##yrsb
plot(syslb,mod3.lme$fit[,1],col='red',ylim=c(50,100),ylab='',xlab='')##fixed effe
lines(lowess(syslb,mod3.lme$fit[,1]),col='red')
points(syslb,mod3.lme$fit[,2],col='blue')##random eff
lines(lowess(syslb,mod3.lme$fit[,2]),col='blue')
lines(lowess(syslb,rich),lwd=2)
legend("bottomleft",c('fixed effect','random effect'),col=c('red','blue'),pch=1)

gyr.avgs<-tapply(rich,list(syr,grazer),mean)
gyr.sd<-tapply(rich,list(syr,grazer),sd)
gyr.n<-tapply(rich,list(syr,grazer),length)
gyr.se<-gyr.sd/sqrt(gyr.n)
plot(1998:2008,gyr.avgs[,1],type='l',col='red',ylim=c(50,100))
arrows(1998:2008,gyr.avgs[,1],1998:2008,gyr.avgs[,1]+gyr.se[,1],angle=90,length=.05)
arrows(1998:2008,gyr.avgs[,1],1998:2008,gyr.avgs[,1]-gyr.se[,1],angle=90,length=.05)
points(1998:2008,gyr.avgs[,2],type='l',col='blue')
arrows(1998:2008,gyr.avgs[,2],1998:2008,gyr.avgs[,2]+gyr.se[,2],angle=90,length=.05)
arrows(1998:2008,gyr.avgs[,2],1998:2008,gyr.avgs[,2]-gyr.se[,2],angle=90,length=.05)
points(1998:2008,tapply(rich,syr,mean),type='l')
legend("bottomright",c('bison','cattle'),col=c('red','blue'),lty=1)


##with richness centered on plot means
plrich<-tapply(rich,name,mean)
crich<-rich-rep(plrich,11)
gyr.avgs<-tapply(crich,list(syr,grazer),mean)
gyr.sd<-tapply(crich,list(syr,grazer),sd)
gyr.n<-tapply(crich,list(syr,grazer),length)
gyr.se<-gyr.sd/sqrt(gyr.n)
plot(1998:2008,gyr.avgs[,1],type='l',col='red',ylim=c(-20,15),xlab="",ylab='',lwd=2)
arrows(1998:2008,gyr.avgs[,1],1998:2008,gyr.avgs[,1]+gyr.se[,1],angle=90,length=.05)
arrows(1998:2008,gyr.avgs[,1],1998:2008,gyr.avgs[,1]-gyr.se[,1],angle=90,length=.05)
points(1998:2008,gyr.avgs[,2],type='l',col='blue',lwd=2)
arrows(1998:2008,gyr.avgs[,2],1998:2008,gyr.avgs[,2]+gyr.se[,2],angle=90,length=.05)
arrows(1998:2008,gyr.avgs[,2],1998:2008,gyr.avgs[,2]-gyr.se[,2],angle=90,length=.05)
legend("bottomright",c('bison','cattle'),col=c('red','blue'),lty=1)

par(mfrow=c(1,3))
plot(syob,crich,xlab="Yrs of Bison",ylab="")
lines(lowess(syob,crich))
plot(sbp5,crich,xlab="# Burns in Past 5 Years",ylab="")
lines(lowess(sbp5,crich))
plot(syslb,crich,xlab="Years Since Last Burn",ylab="")
lines(lowess(syslb,crich))

##none of the above relationships look good so I will look only at within plot effects
par(mfrow=c(1,3))
pls<-unique(name)
plot(syob,crich,xlab="Yrs of Bison",ylab="",type='n',ylim=c(-25,25))
for(i in 1:20){ #plot
 points(syob[name==pls[i]],crich[name==pls[i]],col=i,cex=.75)
 lines(lowess(syob[name==pls[i]],crich[name==pls[i]]),col=i,lwd=2)
}
plot(sbp5,crich,xlab="# Burns in Past 5 Years",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points(sbp5[name==pls[i]],crich[name==pls[i]],col=i,cex=.75)
 lines(lowess(sbp5[name==pls[i]],crich[name==pls[i]]),col=i,lwd=2)
}
plot(syslb,crich,xlab="Years Since Last Burn",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points(syslb[name==pls[i]],crich[name==pls[i]],col=i,cex=.75)
 lines(lowess(syslb[name==pls[i]],crich[name==pls[i]]),col=i,lwd=2)
}



#############################################

##compare the lme model with a model which considers plots in a standard way
mod.pl<-lm(rich~syr+I(syr^2)+name)
mod.clim<-lm(resid(mod.pl)~srain*stemp)
summary(mod.clim)
syr           4.1487     0.5618   7.384 3.33e-12 ***
srain        -2.7245     0.6448  -4.226 3.52e-05 ***
stemp        -1.9807     0.5204  -3.806 0.000184 ***
srain:stemp   4.3452     0.7205   6.031 7.04e-09 ***
Multiple R-squared: 0.07768,    Adjusted R-squared: 0.06486 
plot(mod.clim)
mod.clim2<-update(mod.clim,.~.+I(syr^2))
anova(mod.clim,mod.clim2)
##quadratic yr effect not visible
mod.clim3<-update(mod.clim,.~.+psrain+pstemp)
anova(mod.clim,mod.clim3)
summary(mod.clim3)
mod.clim4<-lm(resid(mod.pl)~syr+srain*psrain)
mod.clim4<-lm(resid(mod.pl)~syr+stemp*pstemp)
mod.clim5<-lm(resid(mod.pl)~syr+srain+psrain+stemp*pstemp+srain:stemp)
summary(mod.clim5)

step.clim<-(lm(resid(mod.pl)~srain+psrain+stemp+pstemp+I(psrain^-1)+I(pstemp^-1)))

summary(step.clim)

##Adler and Levine mods
mod1<-lm(resid(mod.pl)~srain)
mod2<-lm(resid(mod.pl)~psrain)
mod3<-lm(resid(mod.pl)~I(psrain^-1))
anova(mod1,mod2,mod3)
plot(mod3)
mod4<-update(mod3,.~.+srain)
anova(mod3,mod4)
plot(I(psrain^-1),resid(mod.pl))
lines(lowess(I(psrain^-1),resid(mod.pl)))

plot(I(prain.tot^-1),resid(mod.pl))
lines(lowess(I(prain.tot^-1),resid(mod.pl)))

mod1<-lm(resid(mod.pl)~stemp)
mod2<-lm(resid(mod.pl)~pstemp)
mod3<-lm(resid(mod.pl)~I(pstemp^-1))
anova(mod1,mod2,mod3)
plot(mod3)
mod4<-update(mod3,.~.+stemp)
anova(mod3,mod4)
plot(I(pstemp^-1),resid(mod.pl))
lines(lowess(I(pstemp^-1),resid(mod.pl)))
summary(mod3)

##interaction model
mod1<-lm(resid(mod.pl)~psrain*pstemp)
mod2<-lm(resid(mod.pl)~I(psrain^-1)*I(pstemp^-1))
mod3<-lm(resid(mod.pl)~(psrain)*I(pstemp^-1))
mod4<-lm(resid(mod.pl)~psrain+I(pstemp^-1))
anova(mod1,mod2,mod3,mod4)
#mod4 is best
summary(mod4)
scatterplot3d(psrain,I(pstemp^-1),resid(mod.pl))

####
summary(lm(rich~-1+name))




##################looking at z######################
###
mod1.pl<-lm(z~name+syr)
mod2.pl<-update(mod1.pl,.~.+syr2)
mod3.pl<-update(mod2.pl,.~.+I(syr^3))
anova(mod2.pl,mod3.pl)
mod.pl<-mod2.pl
mod1.sea<-lm(resid(mod.pl)~srain+stemp+ssum.rain+sspr.rain+swin.rain+ssum.temp+sspr.temp+swin.temp)
summary(mod1.sea)
ssum.rain   -2.569e-03  3.018e-03    -0.851  0.39555   
sspr.rain    7.284e-05  2.484e-03     0.029  0.97664   
swin.rain   -8.095e-04  2.318e-03    -0.349  0.72729   
ssum.temp    1.086e-03  3.120e-03     0.348  0.72816   
sspr.temp    5.462e-03  2.904e-03     1.881  0.06135 . 
swin.temp   -7.779e-03  2.771e-03    -2.807  0.00546 **
###simplify model 1
mod2.sea<-lm(resid(mod.pl)~swin.temp)
summary(mod2.sea)
anova(mod1.sea,mod2.sea)
##mod2 is superior

