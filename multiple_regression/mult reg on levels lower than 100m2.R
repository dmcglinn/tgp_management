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

