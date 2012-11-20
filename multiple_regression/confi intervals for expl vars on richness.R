endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
library(nlme)
endatlv1<-endat[endat$lv==1&endat$cor==1,] ##to only work at the 100 m2 scale
attach(endatlv1)

test<-lm(rich~PDSIavg+name+dist.m)
test2<-lm(rich~PDSIavg+I(PDSIavg^2)+name+dist.m)
anova(test,test2)
newdata<-data.frame(PDSIavg=sort(unique(PDSIavg)),name=rep("Cowboy",11))
plot(PDSIavg,rich)
lines(newdata$PDSIavg,predict(test2,newdata=newdata),type='l')

test<-gls(rich~PDSIavg+I(PDSIavg^2)+name+dist.m)
test2<-gls(rich~PDSIavg+I(PDSIavg^2)+name+dist.m,cor=corAR1(form=~yr|name))
anova(test,test2)
summary(test2)$t
newdata<-data.frame(PDSIavg=sort(unique(PDSIavg)),name=rep("Cowboy",11))
plot(PDSIavg,rich)

lines(newdata$PDSIavg,predict(test2,newdata=newdata),type='l')

rich.m<-gls(rich~1,cor=corRatio(form=~yr|name,n=T))

dist.m<-cbind(YrsOB,YrsSLB,BP5Yrs)
site.m<-endatlv1[,8:26]
yr.m<-endatlv1[,28:37]

spat.mods(logCa,site.m,yr,name)
ca.m<-gls(logCa~1,cor=corRatio(form=~yr|name))
plot(ACF(ca.m,resType='n'),alpha=0.05)
summary(lm(residuals(rich.m,type='n')~residuals(ca.m,type='n')))


holder<-matrix(NA,ncol=3,nrow=5)
for(i in 1:5){
 y<-endat$rich[endat$lv==i]
 x<-endat$logCa[endat$lv==i]
 mod<-lm(y~x)
 holder[i,1]<-coef(mod)[2]
 holder[i,2]<-summary(mod)$r.s
 holder[i,3]<-anova(mod)$P[1]
}
round(holder,4)
        [,1]   [,2]   [,3]
[1,] -31.6807 0.1829 0.0000
[2,] -21.6390 0.1462 0.0000
[3,] -10.6542 0.0894 0.0000
[4,]  -3.1320 0.0232 0.0000
[5,]  -0.1154 0.0001 0.7348

test<-rda(srs[,-1]~logCa+slope+northness+Condition(dist.m)+Condition(yr.f))
anova(test,by='m')
plot(test,display = c("sp", "cn"),scaling = 1)

mod<-gls(rich~logCa+slope+northness+dist.m+PDSIavg+I(PDSIavg^2),cor=corAR1(form=~yr|name))
spat.mods(rich,cbind(logCa,slope,northness,dist.m,PDSIavg,PDSIavg^2),yr,name)
mod2<-update(mod,cor=corRatio(form=~yr|name,nugget=T))
anova(mod,mod2)
round(summary(mod2)$t,3)
               Value Std.Error t-value p-value
(Intercept)  148.249    24.872   5.961   0.000
logCa        -19.559     7.430  -2.632   0.009
slope         -0.523     0.796  -0.656   0.512
northness      4.061     2.025   2.006   0.046
dist.mYrsOB    0.979     0.300   3.261   0.001
dist.mYrsSLB  -0.770     0.404  -1.904   0.058
dist.mBP5Yrs  -1.658     0.908  -1.826   0.069
PDSIavg        1.658     0.598   2.774   0.006
I(PDSIavg^2)  -1.292     0.277  -4.666   0.000
##standardized variables
mod3<-gls(scale(rich)~scale(logCa)+scale(slope)+scale(northness)+scale(dist.m)+scale(PDSIavg)+scale(PDSIavg^2),cor=corRatio(form=~yr|name,nugget=T))
intervals(mod3)
Approximate 95% confidence intervals
 Coefficients:
                           lower        est.       upper
(Intercept)         -0.271174568 -0.04797443  0.17522572
scale(logCa)        -0.461757302 -0.26403021 -0.06630311
scale(slope)        -0.317170784 -0.07922320  0.15872439
scale(northness)     0.004220176  0.24613090  0.48804162
scale(dist.m)YrsOB   0.128141815  0.32407593  0.52001004
scale(dist.m)YrsSLB -0.261672231 -0.12857843  0.00451538
scale(dist.m)BP5Yrs -0.372897205 -0.17931402  0.01426916
scale(PDSIavg)       0.056344492  0.19477217  0.33319986
scale(PDSIavg^2)    -0.486528839 -0.34202719 -0.19752554
attr(,"label")
[1] "Coefficients:"

 Correlation structure:
           lower      est.      upper
range  3.2063509 6.9443484 15.0401427
nugget 0.3760211 0.4995376  0.6231107
attr(,"label")
[1] "Correlation structure:"

 Residual standard error:
    lower      est.     upper 
0.6787118 0.7796957 0.8957047 

confi<-intervals(mod3)
betas<-confi$coef[-1,2]
low<-confi$coef[-1,1]
up<-confi$coef[-1,3]
plot(c(-1,1),c(1,10),type='n',axes=F,frame.plot=F,xlab='',ylab='',xlim=c(-2,1))
axis(side=1,cex.axis=1.5,lwd=2,at=c(-1,-.5,0,.5,1))
abline(v=0,lty=2,col='grey',lwd=2)
points(betas,9:2,pch=19,cex=1.5)
arrows(betas,9:2,low,9:2,angle=90,len=.1,lwd=2)
arrows(betas,9:2,up,9:2,angle=90,len=.1,lwd=2)
labs<-c('log Ca','slope','northness','years of bison','years since last burn','# burns in 5 years','PDSI',expression(PSDI^2))
text(x=rep(-2,8),y=8.8:1.8,labels=labs,adj=c(0,0),cex=1.5)


plot(ACF(mod2,resType='n'),alpha=0.05)
plot(Variogram(mod2,resTyp='n'),ylim=c(0,2),xlim=c(0,5))

####
spat.mods(rich,cbind(logCa,slope,northness,yr.f,dist.m))
sitmod<-gls(rich~logCa+slope+northness+yr.f+dist.m,cor=corRatio(form=~yr|name,n=T))
anova(sitmod,type='m')
Denom. DF: 203 
            numDF   F-value p-value
(Intercept)     1 29.643582  <.0001
logCa           1  7.573431  0.0065
slope           1  0.761062  0.3840
northness       1  4.037643  0.0458
yr.f           10 12.572956  <.0001
dist.m          3  5.999731  0.0006
summary(sitmod)$t
pR2.gls(sitmod)$r2ML
[1] 0.5434781
pR2diff<-function(modfull,modsub){
 pR2.gls(modfull)$r2ML-pR2.gls(modsub)$r2ML
}
mod.ca<-update(sitmod,.~.-logCa)
mod.slope<-update(sitmod,.~.-slope)
mod.north<-update(sitmod,.~.-northness)

ca.r2<-pR2diff(sitmod,mod.ca)
sl.r2<-pR2diff(sitmod,mod.slope)
no.r2<-pR2diff(sitmod,mod.north)


#################################
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
