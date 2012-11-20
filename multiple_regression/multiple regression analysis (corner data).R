endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
names(endat)
  [1] "plot"           "yr"             "lv"             "cor"           
  [5] "plot.yr.lv.cor" "plotnum"        "yr.plot"        "X205"          
  [9] "X206"           "X208"           "X220"           "X222"          
 [13] "X226"           "X238"           "X244"           "X254"          
 [17] "X259"           "X303"           "X307"           "X308"          
 [21] "X309"           "X317"           "X319"           "X331"          
 [25] "X343"           "X346"           "X350"           "X1998"         
 [29] "X1999"          "X2000"          "X2001"          "X2002"         
 [33] "X2003"          "X2004"          "X2005"          "X2006"         
 [37] "X2007"          "X2008"          "OrBis"          "ChBis"         
 [41] "bison"          "cattle"         "YrsOB"          "BP5Yrs"        
 [45] "YrsSLB"         "burn"           "spring"         "summer"        
 [49] "fall"           "slope"          "northness"      "eastness"      
 [53] "OM"             "pH"             "SolS"           "logP"          
 [57] "logCa"          "logFe"          "logMg"          "logK"          
 [61] "logNa"          "logB"           "logMn"          "logCu"         
 [65] "logZn"          "logAl"          "PDSIavg"        "SodTemp"       
 [69] "SPI1"           "SPI12"          "SPI24"          "rain1"         
 [73] "rain2"          "rain3"          "rain4"          "rain5"         
 [77] "rain6"          "rain7"          "rain8"          "rain9"         
 [81] "rain10"         "rain11"         "rain12"         "temp1"         
 [85] "temp2"          "temp3"          "temp4"          "temp5"         
 [89] "temp6"          "temp7"          "temp8"          "temp9"         
 [93] "temp10"         "temp11"         "temp12"         "sum.rain.tot"  
 [97] "win.rain.tot"   "spr.rain.tot"   "sum.rain.avg"   "win.rain.avg"  
[101] "spr.rain.avg"   "sum.temp"       "win.temp"       "spr.temp"      
[105] "rain.tot"       "rain.avg"       "temp.avg"       "rich1"         
[109] "rich2"          "rich3"          "rich4"          "rich5"         
[113] "A.cov"          "P.cov"          "S.cov"          "T.cov"         
[117] "A.p"            "P.p"            "S.p"            "T.p"           
[121] "ForbLeg.cov"    "Forb.cov"       "Gsum.cov"       "G3.cov"        
[125] "G4.cov"         "Leg.cov"        "Wood.cov"       "F.p"           
[129] "Gsum.p"         "G3.p"           "G4.p"           "Leg.p"         
[133] "Wood.p"         "rich"           "z"              "c"             
[137] "prain1"         "prain2"         "prain3"         "prain4"        
[141] "prain5"         "prain6"         "prain7"         "prain8"        
[145] "prain9"         "prain10"        "prain11"        "prain12"       
[149] "ptemp1"         "ptemp2"         "ptemp3"         "ptemp4"        
[153] "ptemp5"         "ptemp6"         "ptemp7"         "ptemp8"        
[157] "ptemp9"         "ptemp10"        "ptemp11"        "ptemp12"       
[161] "psum.rain.tot"  "pwin.rain.tot"  "pspr.rain.tot"  "psum.rain.avg" 
[165] "pwin.rain.avg"  "pspr.rain.avg"  "psum.temp"      "pwin.temp"     
[169] "pspr.temp"      "prain.tot"      "prain.avg"      "ptemp.avg"

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
srain<-scale(rain.tot)
stemp<-scale(temp.avg)
psrain<-scale(prain.tot)
pstemp<-scale(ptemp.avg)

climgrp<-groupedData(rich~syr+srain*stemp|p1)

mod1.lis<-lmList(climgrp)
plot(intervals(mod1.lis))
plot(mod1.lis)
pairs(mod1.lis,id=0.01)
MSE<-summary(mod1.lis)$RSE/summary(mod1.lis)$df.residual

mod1.lm<-lm(rich~syr+srain*stemp)
summary(mod1.lm)
mod2.lm<-update(mod1.lm,.~.+I(syr^2))
summary(mod2.lm)
anova(mod1.lm,mod2.lm)##drop quadratic time term

##start with mixed model in which each plot has a diff intercept but common slope
mod1.gls<-gls(rich~syr+srain*stemp,data=climgrp)
mod1.lme<-lme(rich~syr+srain*stemp,random=~1|p1,data=climgrp)
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
##a positive diagonal correlation assumes ...
mod3.lme<-update(mod2.lme,random=pdDiag(~syr))
anova(mod2.lme,mod3.lme)
##mod2.lme is superior
summary(mod2.lme)
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
qqnorm(mod4.lme)
qqnorm(mod4.lme,~resid(.)|p1)
##all dianostic plots look reasonably good
random.effects(mod4.lme)
plot(augPred(mod4.lme,~syr,level=0:1),grid=T)
plot(augPred(mod4.lme,~stemp,level=0:1),grid=T)

##compare coef with lmList object
plot(compareFits(coef(mod1.lis),coef(mod4.lme)),mark=fixef(mod4.lme))
##compare pred with lmList object
plot(comparePred(mod1.lis,mod4.lme,length.out=2),layout=c(5,4))
##this last plot does not seem to do what I want

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
anova(mod2.lme,mod6.lme)
mod7.lme<-update(mod2.lme,corr=corLin(form=~syr))
anova(mod2.lme,mod7.lme)
###mod4 still superior

##considering autocorrelation does not seem neccessary
summary(mod4.lme)
Linear mixed-effects model fit by REML
 Data: climgrp 
       AIC      BIC    logLik
  1561.745 1595.451 -770.8723

Random effects:
 Formula: ~syr | p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev     Corr  
(Intercept) 9.48983640 (Intr)
syr         1.10311121 0.516 
Residual    0.03082096       

Variance function:
 Structure: Power of variance covariate
 Formula: ~fitted(.) 
 Parameter estimates:
   power 
1.246848 
Fixed effects: rich ~ syr + srain * stemp 
               Value Std.Error  DF  t-value p-value
(Intercept) 75.69062 2.1722788 196 34.84388   0e+00
syr          1.40138 0.2922498 196  4.79514   0e+00
srain       -2.69874 0.5733595 196 -4.70690   0e+00
stemp       -1.84615 0.4723388 196 -3.90852   1e-04
srain:stemp  4.48344 0.6402459 196  7.00269   0e+00
 Correlation: 
            (Intr) syr    srain  stemp 
syr          0.437                     
srain        0.000  0.000              
stemp       -0.004 -0.080  0.185       
srain:stemp -0.001  0.178 -0.618 -0.366

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-2.29722772 -0.66964437  0.01507455  0.61204947  2.87233656 

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
mod.pl2<-update(mod.pl,.~.+I(syr^2))
anova(mod.pl1,mod.pl2)
##quadratic term very important
mod.pl<-mod.pl2
mod1.sea<-lm(resid(mod.pl)~ssum.rain+sspr.rain+swin.rain+ssum.temp+sspr.temp+swin.temp)
summary(mod1.sea)
ssum.rain   -2.850e+00  6.440e-01    -4.425 1.54e-05 ***
sspr.rain    2.429e+00  5.302e-01     4.582 7.85e-06 ***
swin.rain    2.224e+00  4.947e-01     4.496 1.14e-05 ***
ssum.temp   -6.298e-01  6.659e-01    -0.946 0.345320    
sspr.temp    1.469e-01  6.197e-01     0.237 0.812842    
swin.temp   -2.325e+00  5.914e-01    -3.932 0.000114 ***
mod.sea<-lm(resid(mod.pl)~syr*ssum.rain*sspr.rain*swin.rain*ssum.temp*sspr.temp*swin.temp)
summary(step(mod.sea))
###simplify model 1
mod2.sea<-update(mod1.sea,.~.-ssum.temp-sspr.temp)
summary(mod2.sea)
anova(mod1.sea,mod2.sea)
##mod2 is superior
summary(mod2.sea)
Coefficients:
              Estimate Std. Error   t value Pr(>|t|)    
(Intercept) -6.853e-16  4.578e-01 -1.50e-15        1    
ssum.rain   -2.419e+00  4.786e-01    -5.054 9.25e-07 ***
sspr.rain    2.202e+00  4.608e-01     4.778 3.28e-06 ***
swin.rain    2.169e+00  4.794e-01     4.524 1.00e-05 ***
swin.temp   -2.351e+00  4.629e-01    -5.079 8.23e-07 ***
Residual standard error: 6.79 on 215 degrees of freedom
Multiple R-squared: 0.2518,     Adjusted R-squared: 0.2379 
F-statistic: 18.09 on 4 and 215 DF,  p-value: 8.023e-13

##build lmList model
syr2<-syr^2
clim.sea<-groupedData(rich~syr+syr2+ssum.rain+sspr.rain+swin.rain+swin.temp|p1)
mod1.lis<-lmList(clim.sea)
plot(intervals(mod1.lis))##looks like common slopes will work except for yr
plot(mod1.lis)
pairs(mod1.lis,id=0.01)
summary(mod1.lis)$RSE
summary(mod1.lis)$df.residual
##build mixed effect model
mod1.lme<-lme(rich~syr+syr2+ssum.rain+sspr.rain+swin.rain+swin.temp,random=~1|p1,data=clim.sea)
summary(mod1.lme)
mod2.lme<-update(mod1.lme,random=~syr|p1)
anova(mod1.lme,mod2.lme)
##noncommon slopes for year effect superior
summary(mod2.lme)
Random effects:
 Formula: ~syr | p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev   Corr  
(Intercept) 9.598010 (Intr)
syr         1.094839 0.501 
Residual    6.149255   

Fixed effects: rich ~ syr + ssum.rain + sspr.rain + swin.rain + swin.temp 
               Value Std.Error  DF  t-value p-value
(Intercept) 79.16273 2.2541774 194 35.11823   0e+00
syr          0.96695 0.2826551 194  3.42096   8e-04
syr2        -0.29355 0.0550769 194 -5.32974   0e+00
ssum.rain   -2.22677 0.4581761 194 -4.86008   0e+00
sspr.rain    2.44521 0.4549883 194  5.37423   0e+00
swin.rain    2.32250 0.4632685 194  5.01329   0e+00
swin.temp   -2.34636 0.4193867 194 -5.59475   0e+00

##diagnostic plots
plot(mod2.lme)
plot(mod2.lme,resid(.)~syr,abline=0)
plot(mod2.lme,resid(.)~ssum.rain,abline=0)
plot(mod2.lme,resid(.)~sspr.rain,abline=0)
plot(mod2.lme,resid(.)~sspr.rain,abline=0)

qqnorm(mod2.lme)
##all dianostic plots look reasonably good
plot(mod2.lme,resid(.)~fitted(.)|p1) ##plots do not appear to have diff variances
plot(mod2.lme,resid(.,type='n')~fitted(.)|p1) ##plots do not appear to have diff variances
random.effects(mod2.lme)
plot(augPred(mod2.lme,~syr),aspect='xy',grid=T)
plot(augPred(mod2.lme,~ssum.rain),grid=T)
plot(augPred(mod2.lme,~swin.rain),grid=T)
plot(augPred(mod2.lme,~sspr.rain),grid=T)

##examine a model allowing each plot to have a different residual variance
mod4.lme<-update(mod2.lme,weights=varIdent(form=~1|p1))
anova(mod2.lme,mod4.lme)
##mod2 with one variance for all plots is superior

plot(ACF(mod2.lme,maxLag=5),alpha=0.01)
varo<-Variogram(mod2.lme,maxDist=5,form=~syr)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')
##variogram indicates no autocorrelation

##there does not appear to be a large degree of autocorrelation
mod4.lme<-update(mod2.lme,corr=corARMA(p=1,form=~syr))
mod5.lme<-update(mod2.lme,corr=corARMA(p=2,form=~syr))##does not converge
plot(ACF(mod4.lme,maxLag=5),alpha=0.01)
anova(mod2.lme,mod4.lme)
##model 4 which has a single autoregressive term is superior
mod6.lme<-update(mod2.lme,corr=corLin(form=~syr))
anova(mod4.lme,mod6.lme)
###mod4 still superior
summary(mod4.lme)
Random effects:
 Formula: ~syr | p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev    Corr  
(Intercept) 9.4300710 (Intr)
syr         0.9410423 0.567 
Residual    6.5020154       

Correlation Structure: AR(1)
 Formula: ~syr | p1 
 Parameter estimate(s):
      Phi 
0.2106036 
Fixed effects: rich ~ syr + syr2 + ssum.rain + sspr.rain + swin.rain + swin.temp 
               Value Std.Error  DF  t-value p-value
(Intercept) 79.12972 2.2682430 194 34.88591   0e+00
syr          0.99418 0.2692941 194  3.69181   3e-04
syr2        -0.28968 0.0609784 194 -4.75057   0e+00
ssum.rain   -2.13522 0.4441266 194 -4.80769   0e+00
sspr.rain    2.33685 0.4570053 194  5.11339   0e+00
swin.rain    2.29969 0.4425371 194  5.19660   0e+00
swin.temp   -2.30413 0.3888834 194 -5.92498   0e+00

###compare season model to total rain and temp + interaction model
mod1.lme<-lme(rich~syr+syr2+srain*stemp,random=~syr|p1,corr=corAR1(form=~syr),method="ML")
mod2.lme<-lme(rich~syr+syr2+ssum.rain+sspr.rain+swin.rain+swin.temp,random=~syr|p1,corr=corAR1(form=~syr),method="ML")
##mod1.lme is not quite what we found to be the best model above
##it has been fudged to make the comparison clearer
anova(mod1.lme,mod2.lme)
##mod2.lme with the seasonal breakdown does a much better job
##now lets check the mod1.lme that we actually know was the best
mod1.lme<-lme(rich~syr+srain*stemp,random=~syr|p1,method="ML")
anova(mod1.lme,mod2.lme)
##mod2.lme still the best even with the bigger diff in d.f.s
detach(endatalv1)
#######################################################################
##clim vars at lv 3
attach(endat)
endatlv3<-endat[endat$lv==3,] ##to only work at the 1 m2 scale
detach(endat)
attach(endatlv3)
c1<-as.factor(cor)
p1<-as.factor(plot)
syr<-scale(yr,scale=F)
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


mod.pl1<-lm(rich~p1+syr)
mod.pl2<-update(mod.pl1,.~.+I(syr^2))
mod.pl3<-update(mod.pl2,.~.+I(syr^3))
mod.pl4<-update(mod.pl3,.~.+I(syr^4))
anova(mod.pl1,mod.pl2,mod.pl3,mod.pl4)
##3 degree polynomial term very important
mod.pl<-mod.pl3
mod1.sea<-lm(resid(mod.pl)~ssum.rain+sspr.rain+swin.rain+ssum.temp+sspr.temp+swin.temp)
summary(mod1.sea)
ssum.rain   -1.168e+00  2.314e-01    -5.047 5.46e-07 ***
sspr.rain    7.106e-01  1.905e-01     3.730 0.000203 ***
swin.rain    8.760e-01  1.777e-01     4.929 9.88e-07 ***
ssum.temp   -2.006e-01  2.392e-01    -0.838 0.402024    
sspr.temp    6.025e-01  2.227e-01     2.706 0.006942 ** 
swin.temp   -1.067e+00  2.125e-01    -5.022 6.20e-07 ***
Residual standard error: 6.645 on 873 degrees of freedom
Multiple R-squared: 0.2094,     Adjusted R-squared: 0.2039 
F-statistic: 38.53 on 6 and 873 DF,  p-value: < 2.2e-16 
###simplify model 1
mod2.sea<-update(mod1.sea,.~.-ssum.temp)
summary(mod2.sea)
anova(mod1.sea,mod2.sea)
##mod2 is superior
##build lmList model
syr2<-syr^2
syr3<-syr^3
syr4<-syr^4
clim.sea<-groupedData(rich~syr+syr2+ssum.rain+sspr.rain+swin.rain+swin.temp|p1/c1)
mod1.lis<-lmList(clim.sea)
plot(intervals(mod1.lis))##looks like common slopes will work except for yr
plot(mod1.lis)
pairs(mod1.lis,id=0.01)
summary(mod1.lis)$RSE
summary(mod1.lis)$df.residual
##build mixed effect model
mod1.lme<-lme(rich3~syr+syr2+ssum.rain+sspr.rain+swin.rain+swin.temp,random=~1|p1,data=clim.sea)
mod2.lme<-update(mod1.lme,random=~1|p1/c1)
anova(mod1.lme,mod2.lme)
summary(mod1.lme)
mod2.lme<-update(mod1.lme,random=~syr|p1/c1)
anova(mod1.lme,mod2.lme)
##noncommon slopes for year effect superior
summary(mod2.lme)
Random effects:
 Formula: ~syr | p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev    Corr  
(Intercept) 3.4629076 (Intr)
syr         0.5234711 0.156 

 Formula: ~syr | c1 %in% p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev       Corr  
(Intercept) 6.881739e-05 (Intr)
syr         1.604078e-05 0     
Residual    2.593479e+00       

Fixed effects: rich3 ~ syr + syr2 + ssum.rain + sspr.rain + swin.rain + swin.temp 
                Value Std.Error  DF   t-value p-value
(Intercept) 24.226798 0.7878575 794 30.750230  0.0000
syr          0.066960 0.1207838 794  0.554375  0.5795
syr2        -0.042225 0.0116145 794 -3.635572  0.0003
ssum.rain   -0.864814 0.0964540 794 -8.966079  0.0000
sspr.rain    1.047895 0.0957829 794 10.940310  0.0000
swin.rain    1.036926 0.0975260 794 10.632305  0.0000
swin.temp   -0.802554 0.0882882 794 -9.090162  0.0000

##diagnostic plots
plot(mod2.lme)
plot(mod2.lme,resid(.)~syr,abline=0)
plot(mod2.lme,resid(.)~ssum.rain,abline=0)
plot(mod2.lme,resid(.)~sspr.rain,abline=0)
plot(mod2.lme,resid(.)~sspr.rain,abline=0)

qqnorm(mod2.lme)
##all dianostic plots look reasonably good
plot(mod2.lme,resid(.)~fitted(.)|p1) ##plots do not appear to have diff variances
plot(mod2.lme,resid(.,type='n')~fitted(.)|p1) ##plots do not appear to have diff variances
random.effects(mod2.lme)
plot(augPred(mod2.lme,~syr),aspect='xy',grid=T)
plot(augPred(mod2.lme,~ssum.rain),grid=T)
plot(augPred(mod2.lme,~swin.rain),grid=T)
plot(augPred(mod2.lme,~sspr.rain),grid=T)

##examine a model allowing each plot to have a different residual variance
mod4.lme<-update(mod2.lme,weights=varIdent(form=~1|p1))
anova(mod2.lme,mod4.lme)
##mod2 with one variance for all plots is superior

plot(ACF(mod2.lme,maxLag=5),alpha=0.01)
varo<-Variogram(mod2.lme,maxDist=5,form=~syr)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')
##variogram indicates no autocorrelation

##there does not appear to be a large degree of autocorrelation
mod4.lme<-update(mod2.lme,corr=corARMA(p=1,form=~syr))
mod5.lme<-update(mod2.lme,corr=corARMA(p=2,form=~syr))##does not converge
plot(ACF(mod4.lme,maxLag=5),alpha=0.01)
anova(mod2.lme,mod4.lme)
##model 4 which has a single autoregressive term is superior
mod6.lme<-update(mod2.lme,corr=corLin(form=~syr))
anova(mod4.lme,mod6.lme)
###mod4 still superior
summary(mod4.lme)
Random effects:
 Formula: ~syr | p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev    Corr  
(Intercept) 3.4911522 (Intr)
syr         0.5359877 0.156 

 Formula: ~syr | c1 %in% p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev       Corr  
(Intercept) 6.526067e-05 (Intr)
syr         1.961388e-05 0     
Residual    2.581380e+00       

Correlation Structure: AR(1)
 Formula: ~syr | p1/c1 
 Parameter estimate(s):
       Phi 
-0.1399773 
Fixed effects: rich3 ~ syr + syr2 + ssum.rain + sspr.rain + swin.rain + swin.temp 
                Value Std.Error  DF   t-value p-value
(Intercept) 24.222167 0.7916647 794 30.596498  0.0000
syr          0.071524 0.1228900 794  0.582014  0.5607
syr2        -0.041469 0.0110349 794 -3.758018  0.0002
ssum.rain   -0.847939 0.0989468 794 -8.569647  0.0000
sspr.rain    1.035686 0.0948957 794 10.913940  0.0000
swin.rain    1.036815 0.1020597 794 10.158898  0.0000
swin.temp   -0.780436 0.0932838 794 -8.366251  0.0000

#####lv5###################
detach(endatlv3)
attach(endat)
endatlv5<-endat[endat$lv==5,] ##to only work at the .01 m2 scale
detach(endat)
attach(endatlv3)
c1<-as.factor(cor)
p1<-as.factor(plot)
syr<-scale(yr,scale=F)
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


mod.pl1<-lm(rich~p1+syr)
mod.pl2<-update(mod.pl1,.~.+I(syr^2))
mod.pl3<-update(mod.pl2,.~.+I(syr^3))
mod.pl4<-update(mod.pl3,.~.+I(syr^4))
anova(mod.pl1,mod.pl2,mod.pl3,mod.pl4)
anova(mod.pl1,mod.pl3,mod.pl4)
##no trend through time
mod.pl<-lm(rich~p1)
anova(mod.pl1,mod.pl)
mod1.sea<-lm(resid(mod.pl)~ssum.rain+sspr.rain+swin.rain+ssum.temp+sspr.temp+swin.temp)
summary(mod1.sea)
ssum.rain   -1.208e+00  2.310e-01    -5.228 2.14e-07 ***
sspr.rain    9.670e-01  1.902e-01     5.085 4.50e-07 ***
swin.rain    8.898e-01  1.774e-01     5.015 6.43e-07 ***
ssum.temp   -1.773e-01  2.389e-01    -0.742    0.458    
sspr.temp    3.227e-01  2.223e-01     1.452    0.147    
swin.temp   -9.460e-01  2.121e-01    -4.460 9.27e-06 ***
Residual standard error: 4.891 on 873 degrees of freedom
Multiple R-squared: 0.0941,     Adjusted R-squared: 0.08787 
F-statistic: 15.11 on 6 and 873 DF,  p-value: < 2.2e-16
###simplify model 1
mod2.sea<-update(mod1.sea,.~.-ssum.temp-sspr.temp)
summary(mod2.sea)
anova(mod1.sea,mod2.sea)
##mod2 is superior
##build lmList model
syr2<-syr^2
syr3<-syr^3
syr4<-syr^4
clim.sea<-groupedData(rich~ssum.rain+sspr.rain+swin.rain+swin.temp|p1/c1)
mod1.lis<-lmList(clim.sea)
plot(intervals(mod1.lis))##looks like common slopes will work except for yr
plot(mod1.lis)
pairs(mod1.lis,id=0.01)
summary(mod1.lis)$RSE
summary(mod1.lis)$df.residual
##build mixed effect model
mod1.lme<-lme(rich~ssum.rain+sspr.rain+swin.rain+swin.temp,random=~1|p1,data=clim.sea)
summary(mod1.lme)
mod2.lme<-update(mod1.lme,random=~1|p1/c1)
anova(mod1.lme,mod2.lme)
##considering the corner random effects is superior
summary(mod2.lme)
mod3.lme<-update(mod2.lme,random=~swin.rain|p1/c1)
anova(mod1.lme,mod2.lme,mod3.lme)
##simple noncommon intercpet for corners within plots superior
summary(mod2.lme)
Linear mixed-effects model fit by REML
 Data: clim.sea 
       AIC      BIC    logLik
  5315.863 5354.056 -2649.931

Random effects:
 Formula: ~1 | p1
        (Intercept)
StdDev:    3.220822

 Formula: ~1 | c1 %in% p1
        (Intercept) Residual
StdDev:    2.285751 4.521957

Fixed effects: rich ~ ssum.rain + sspr.rain + swin.rain + swin.temp 
                Value Std.Error  DF   t-value p-value
(Intercept) 23.804545 0.7792493 796 30.548049       0
ssum.rain   -0.992631 0.1590863 796 -6.239574       0
sspr.rain    0.922860 0.1531770 796  6.024793       0
swin.rain    0.909955 0.1593590 796  5.710098       0
swin.temp   -0.800842 0.1538609 796 -5.204971       0
##diagnostic plots
##assumption 1
##Are the within group errors centered around zero?
plot(mod2.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
##look at same assumption on a group basis
plot(mod2.lme,p1~resid(.),abline=0)##centered around zero and roughly all same var
##Is there the same degree of variability in each group
plot(mod2.lme,resid(.,type='p')~fitted(.)|p1,id=0.05,adj=-0.3)
##because some plots seem to have slightly different variance
##we can consider a sep variance pred for each plot
mod3.lme<-update(mod2.lme,weights=varPower())##herto follows power funct
mod4.lme<-update(mod2.lme,weights=varIdent(form=~1|p1))
anova(mod2.lme,mod3.lme,mod4.lme)
##seperate variance for each plot seems best
plot(mod4.lme,abline=0,id=0.05,adj=-.3) ##looks better
plot(mod4.lme,p1~resid(.),abline=0)##centered around zero and roughly all same var
plot(mod4.lme,resid(.,type='p')~fitted(.)|p1,id=0.05,adj=-0.3)

##examine difference in response and fitted 
plot(mod4.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qqnorm(mod4.lme)
qqnorm(mod4.lme,~resid(.)|p1)
##assumption 2
##qqnorm() on ranef(.) for checking marginal normality and identifyin outliers
##pairs() on ranef(.) for id outliers and cheing assump of homogenity of ranef covariance matrix
qqnorm(mod4.lme,~ranef(.,level=1),id=0.10,cex=.7)
qqnorm(mod4.lme,~ranef(.,level=2),id=0.10,cex=.7)
#pairs(mod4.lme,~ranef(.,level=1)|p1)#inthis case does not make sense
##for multilevel models must assess random effects are each level of grouping
##all dianostic plots look reasonably good

random.effects(mod2.lme)
plot(augPred(mod2.lme,~syr),aspect='xy',grid=T)
plot(augPred(mod2.lme,~syr),grid=T)
##compare coef with lmList object
plot(compareFits(coef(mod1.lis),coef(mod2.lme)),mark=fixef(mod2.lme))
##compare pred with lmList object
plot(comparePred(mod1.lis,mod2.lme,length.out=2),layout=c(5,4))
##this last plot does not seem to do what I want

plot(mod2.lme,resid(.)~fitted(.)|p1) ##plots do not appear to have diff variances
plot(mod2.lme,resid(.,type='n')~fitted(.)|p1) ##plots do not appear to have diff variances
random.effects(mod2.lme)
plot(augPred(mod2.lme,~syr),aspect='xy',grid=T)
plot(augPred(mod2.lme,~ssum.rain),grid=T)
plot(augPred(mod2.lme,~swin.rain),grid=T)
plot(augPred(mod2.lme,~sspr.rain),grid=T)

##examine a model allowing each plot to have a different residual variance
mod4.lme<-update(mod2.lme,weights=varIdent(form=~1|p1))
anova(mod2.lme,mod4.lme)
##mod2 with one variance for all plots is superior

plot(ACF(mod2.lme,maxLag=5),alpha=0.01)
varo<-Variogram(mod2.lme,maxDist=5,form=~syr)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')
##variogram indicates no autocorrelation

##there does not appear to be a large degree of autocorrelation
mod4.lme<-update(mod2.lme,corr=corARMA(p=1,form=~syr))
mod5.lme<-update(mod2.lme,corr=corARMA(p=2,form=~syr))##does not converge
plot(ACF(mod4.lme,maxLag=5),alpha=0.01)
anova(mod2.lme,mod4.lme)
##model 4 which has a single autoregressive term is superior
mod6.lme<-update(mod2.lme,corr=corLin(form=~syr))
anova(mod4.lme,mod6.lme)
###mod4 still superior
summary(mod4.lme)
Random effects:
 Formula: ~syr | p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev    Corr  
(Intercept) 3.4911522 (Intr)
syr         0.5359877 0.156 

 Formula: ~syr | c1 %in% p1
 Structure: General positive-definite, Log-Cholesky parametrization
            StdDev       Corr  
(Intercept) 6.526067e-05 (Intr)
syr         1.961388e-05 0     
Residual    2.581380e+00       

Correlation Structure: AR(1)
 Formula: ~syr | p1/c1 
 Parameter estimate(s):
       Phi 
-0.1399773 
Fixed effects: rich3 ~ syr + syr2 + ssum.rain + sspr.rain + swin.rain + swin.temp 
                Value Std.Error  DF   t-value p-value
(Intercept) 24.222167 0.7916647 794 30.596498  0.0000
syr          0.071524 0.1228900 794  0.582014  0.5607
syr2        -0.041469 0.0110349 794 -3.758018  0.0002
ssum.rain   -0.847939 0.0989468 794 -8.569647  0.0000
sspr.rain    1.035686 0.0948957 794 10.913940  0.0000
swin.rain    1.036815 0.1020597 794 10.158898  0.0000
swin.temp   -0.780436 0.0932838 794 -8.366251  0.0000

##############disturbance variables##############
##begin looking at 100m2 scale
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)
c1<-as.factor(cor)
p1<-as.factor(plot)
syr<-scale(yr,scale=F)##center's year
syob<-scale(YrsOB,scale=T)##b/c a continuous time var
sbp5<-scale(BP5Yrs,scale=FALSE)##
syslb<-scale(YrsSLB,scale=T)##b/c a cont var
grazer<-factor(ifelse(bison==1,'bison','cattle'))

##build a stardard lm model to look at functional form
mod.p1<-lm(rich~syr+p1)
mod.p2<-update(mod.p1,.~.+I(syr^2))
mod.p3<-update(mod.p2,.~.+I(syr^3))
anova(mod.p1,mod.p2,mod.p3)
mod.pl<-mod.p2
mod1.lm<-lm(resid(mod.pl)~bison)
summary(mod1.lm)
##no bison effect
##fit a stepwise for fun
mod.full<-update(mod1.lm,.~.+bison*syob*sbp5*syslb)
summary(step(mod.full))
##indicates that bison,syob may be the only important main effects
mod2.lm<-update(mod1.lm,.~.+syob+sbp5+syslb)
summary(mod2.lm)
anova(mod1.lm,mod2.lm)
##mod2 is sig better
mod3.lm<-update(mod2.lm,.~.+syob:sbp5)
summary(mod3.lm)
anova(mod2.lm,mod3.lm)


mod1.lm<-lm(resid(mod.pl)~-1+grazer)
summary(mod1.lm)
##no bison effect
##fit a stepwise for fun
mod.full<-update(mod1.lm,.~.+grazer*syob*sbp5*syslb)
summary(step(mod.full))
##indicates that bison,syob may be the only important main effects
mod2.lm<-update(mod1.lm,.~.+syob+sbp5+syslb)
summary(mod2.lm)
             Estimate Std. Error t value Pr(>|t|)   
grazerbison   -2.5273     0.9078  -2.784  0.00585 **
grazercattle   2.8708     0.9944   2.887  0.00429 **
syob           2.1459     0.7775   2.760  0.00628 **
sbp5          -1.5411     0.5890  -2.616  0.00952 **
syslb         -2.1174     0.7170  -2.953  0.00350 **
anova(mod1.lm,mod2.lm)
##mod2 is sig better
mod3.lm<-update(mod2.lm,.~.+syob:sbp5)
summary(mod3.lm)
anova(mod2.lm,mod3.lm)
##marginally better,and interaction term is diff to interpt

syr2<-syr^2
##now fit lmList mode
dobject<-groupedData(rich1~syr+syr2+syob+sbp5+syslb|p1)
mod1.lis<-lmList(dobject)
plot(intervals(mod1.lis))
plot(mod1.lis)##heteroscadistic
pairs(mod1.lis,id=0)
summary(mod1.lis)$RSE
summary(mod1.lis)$df.residual
##build mixed effect model
mod1.lme<-lme(rich~syr+syr2+syob+grazer+syob+sbp5+syslb,random=~1|p1,data=dobject)
summary(mod1.lme)
mod2.lme<-update(mod1.lme,random=~syr|p1)
anova(mod1.lme,mod2.lme)
##just considering intercepts is enough
mod3.lme<-update(mod1.lme,weights=varPower())
anova(mod1.lme,mod3.lme)
##considering heteroscadisticity is helpful
summary(mod3.lme)
Linear mixed-effects model fit by REML
 Data: dobject 
       AIC      BIC    logLik
  1573.773 1607.386 -776.8865
Random effects:
 Formula: ~1 | p1
        (Intercept)   Residual
StdDev:    9.130343 0.03152518

Variance function:
 Structure: Power of variance covariate
 Formula: ~fitted(.) 
 Parameter estimates:
   power 
1.264495 
Fixed effects: rich ~ syr + syr2 + syob + grazer + syob + sbp5 + syslb 
                Value Std.Error  DF  t-value p-value
(Intercept)  76.00287 2.3449860 194 32.41080  0.0000
syr           0.32932 0.2106197 194  1.56358  0.1195
syr2         -0.27766 0.0551973 194 -5.03029  0.0000
syob          5.47503 1.2807922 194  4.27472  0.0000
grazercattle  6.10237 2.0100388 194  3.03595  0.0027
sbp5         -3.11036 0.7428572 194 -4.18702  0.0000
syslb        -2.64988 0.7865268 194 -3.36910  0.0009
 Correlation: 
             (Intr) syr    syr2   syob   grzrct sbp5  
syr          -0.074                                   
syr2         -0.203  0.104                            
syob         -0.052 -0.596 -0.086                     
grazercattle -0.363  0.136 -0.123  0.218              
sbp5          0.035  0.100  0.135  0.037 -0.202       
syslb        -0.017 -0.102  0.070  0.120 -0.030  0.583
##simple noncommon intercpet for corners within plots and varPower superior

##diagnostic plots
##assumption 1
##Are the within group errors centered around zero?
plot(mod3.lme,abline=0,id=0.05,adj=-.3) ##yes but appear heteroscadistic
##look at same assumption on a group basis
plot(mod3.lme,p1~resid(.),abline=0)##centered around zero and roughly all same var
##Is there the same degree of variability in each group
plot(mod3.lme,resid(.,type='p')~fitted(.)|p1,id=0.05,adj=-0.3)
##because some plots seem to have slightly different variance
##we can consider a sep variance pred for each plot
mod4.lme<-update(mod2.lme,weights=varIdent(form=~1|p1))
anova(mod3.lme,mod4.lme)
##examine difference in response and fitted 
plot(mod3.lme,rich~fitted(.),id=0.05,adj=-0.3,abline=c(0,1))
##Normality within group error
qqnorm(mod3.lme)
qqnorm(mod3.lme,~resid(.)|p1)
##assumption 2
##qqnorm() on ranef(.) for checking marginal normality and identifyin outliers
##pairs() on ranef(.) for id outliers and cheing assump of homogenity of ranef covariance matrix
qqnorm(mod3.lme,~ranef(.,level=1),id=0.10,cex=.7)
#pairs(mod3.lme,~ranef(.,level=1)|p1)#inthis case does not make sense
##for multilevel models must assess random effects are each level of grouping
##all dianostic plots look reasonably good

random.effects(mod3.lme)

plot(ACF(mod3.lme,maxLag=5),alpha=0.01)
varo<-Variogram(mod2.lme,maxDist=5,form=~syr)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')
##variogram indicates no autocorrelation

##there does not appear to be a large degree of autocorrelation
mod4.lme<-update(mod2.lme,corr=corARMA(p=1,form=~syr))
mod5.lme<-update(mod2.lme,corr=corARMA(p=2,form=~syr))##does not converge
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
plrich<-tapply(rich,p1,mean)
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
pls<-unique(p1)
plot(syob,crich,xlab="Yrs of Bison",ylab="",type='n',ylim=c(-25,25))
for(i in 1:20){ #plot
 points(syob[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(syob[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}
plot(sbp5,crich,xlab="# Burns in Past 5 Years",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points(sbp5[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(sbp5[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}
plot(syslb,crich,xlab="Years Since Last Burn",ylab="",ylim=c(-25,25))
for(i in 1:20){ #plot
 points(syslb[p1==pls[i]],crich[p1==pls[i]],col=i,cex=.75)
 lines(lowess(syslb[p1==pls[i]],crich[p1==pls[i]]),col=i,lwd=2)
}



#############################################

##compare the lme model with a model which considers plots in a standard way
mod.pl<-lm(rich~syr+I(syr^2)+p1)
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
summary(lm(rich~-1+p1))




##################looking at z######################
###
mod1.pl<-lm(z~p1+syr)
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

