endat<-read.table('C:/Users/dmcglinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
library(nlme)
library(car)
endat1<-endat[endat$lv==1&endat$cor==1,]
endat3<-endat[endat$lv==3,]
endat5<-endat[endat$lv==5,]

mod1<-gls(rich~plot+sum.rain.tot + win.rain.tot + spr.rain.tot,cor=corAR1(form=~yr|plot),data=endat1)
mod3<-gls(rich~plot+cor+sum.rain.tot + win.rain.tot + spr.rain.tot,cor=corAR1(form=~yr|plot/cor),data=endat3)
mod5<-gls(rich~plot+cor+sum.rain.tot + win.rain.tot + spr.rain.tot,cor=corAR1(form=~yr|plot/cor),data=endat5)
coef(mod1)
coef(mod3)
coef(mod5)

par(mfrow=c(2,2))
termplot(lm(rich~name+sum.rain.tot + win.rain.tot + spr.rain.tot,data=endat1),partial.resid=TRUE)
windows()
par(mfrow=c(2,3))
termplot(lm(rich~name+cor+sum.rain.tot + win.rain.tot + spr.rain.tot,data=endat3),partial.resid=TRUE)
windows()
par(mfrow=c(2,3))
termplot(lm(rich~name+cor+sum.rain.tot + win.rain.tot + spr.rain.tot,data=endat5),partial.resid=TRUE)



detach(endat)
endat.no1<-endat[endat$lv!=1,]
attach(endat.no1)
dist.m<-cbind(YrsOB,BP5Yrs,YrsSLB)
full.mod<-lme(rich~dist.m+lv*sum.rain.tot + lv* win.rain.tot +lv*spr.rain.tot,cor=corAR1(form=~yr|plot/cor/lv),random=~1|plot/cor/lv,data=endat.no1)
round(summary(full.mod)$tTable,3)
                  Value Std.Error   DF t-value p-value
(Intercept)      66.314     1.696 3191  39.104   0.000
dist.mYrsOB       0.354     0.043 3191   8.315   0.000
dist.mBP5Yrs     -0.571     0.135 3191  -4.241   0.000
dist.mYrsSLB     -0.574     0.059 3191  -9.761   0.000
lv              -13.318     0.414  239 -32.205   0.000
sum.rain.tot     -0.234     0.018 3191 -13.349   0.000
win.rain.tot      0.142     0.023 3191   6.102   0.000
spr.rain.tot      0.226     0.025 3191   9.060   0.000
lv:sum.rain.tot   0.051     0.005 3191  10.708   0.000
lv:win.rain.tot  -0.020     0.006 3191  -3.105   0.002
lv:spr.rain.tot  -0.043     0.007 3191  -6.262   0.000
summary(full.mod)
Linear mixed-effects model fit by REML
 Data: NULL 
       AIC      BIC    logLik
  21025.76 21124.36 -10496.88

Random effects:
 Formula: ~1 | plot
        (Intercept)
StdDev:    2.854503

 Formula: ~1 | cor %in% plot
         (Intercept)
StdDev: 0.0009010382

 Formula: ~1 | lv %in% cor %in% plot
        (Intercept) Residual
StdDev:    4.294188 4.531262

Correlation Structure: AR(1)
 Formula: ~yr | plot/cor/lv 
 Parameter estimate(s):
      Phi 
0.3113578 
Fixed effects: rich ~ dist.m + lv * sum.rain.tot + lv * win.rain.tot + lv *      spr.rain.tot 
                    Value Std.Error   DF   t-value p-value
(Intercept)      66.31373 1.6958257 3191  39.10409  0.0000
dist.mYrsOB       0.35417 0.0425964 3191   8.31458  0.0000
dist.mBP5Yrs     -0.57087 0.1346222 3191  -4.24056  0.0000
dist.mYrsSLB     -0.57403 0.0588116 3191  -9.76053  0.0000
lv              -13.31766 0.4135279  239 -32.20498  0.0000
sum.rain.tot     -0.23374 0.0175099 3191 -13.34918  0.0000
win.rain.tot      0.14234 0.0233267 3191   6.10218  0.0000
spr.rain.tot      0.22648 0.0249967 3191   9.06040  0.0000
lv:sum.rain.tot   0.05101 0.0047634 3191  10.70809  0.0000
lv:win.rain.tot  -0.01970 0.0063434 3191  -3.10546  0.0019
lv:spr.rain.tot  -0.04251 0.0067894 3191  -6.26190  0.0000
 Correlation: 
                (Intr) ds.YOB d.BP5Y d.YSLB lv     sm.rn. wn.rn. spr.r. lv:sm.. lv:w..
dist.mYrsOB     -0.127                                                                
dist.mBP5Yrs    -0.216  0.177                                                         
dist.mYrsSLB    -0.150  0.094  0.580                                                  
lv              -0.853  0.000  0.000  0.000                                           
sum.rain.tot    -0.386  0.024 -0.013 -0.013  0.411                                    
win.rain.tot    -0.334  0.037 -0.006  0.005  0.353 -0.019                             
spr.rain.tot    -0.556 -0.025  0.022 -0.030  0.590  0.005  0.068                      
lv:sum.rain.tot  0.369  0.000  0.000  0.000 -0.432 -0.952  0.019 -0.005               
lv:win.rain.tot  0.316  0.000  0.000  0.000 -0.370  0.019 -0.952 -0.067 -0.020        
lv:spr.rain.tot  0.530  0.000  0.000  0.000 -0.621 -0.005 -0.067 -0.951  0.006   0.070

Standardized Within-Group Residuals:
         Min           Q1          Med           Q3          Max 
-5.086806102 -0.497554023  0.003558986  0.479292842  5.011359797 

Number of Observations: 3520
Number of Groups: 
                 plot         cor %in% plot lv %in% cor %in% plot 
                   20                    80                   320 
