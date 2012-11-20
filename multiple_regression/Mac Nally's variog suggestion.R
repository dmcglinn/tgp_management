endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
attach(endat)
##begin working only at the 100 m2 scale so no corners
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)

library(nlme)

mod.ar1<-gls(rich~1,cor=corAR1(form=~yr|name))
mod.exp<-update(mod.ar1,.~.,cor=corExp(form=~yr|name,nugget=TRUE))
mod.lin<-update(mod.ar1,.~.,cor=corLin(form=~yr|name,nugget=TRUE))
mod.gau<-update(mod.ar1,.~.,cor=corGaus(form=~yr|name,nugget=TRUE))
mod.sph<-update(mod.ar1,.~.,cor=corSpher(form=~yr|name,nugget=TRUE))
mod.rat<-update(mod.ar1,.~.,cor=corRatio(form=~yr|name,nugget=T))

anova(mod.ar1,mod.rat)

plot(Variogram(mod.rat))
plot(Variogram(mod.sph))
plot(Variogram(mod.gau))

plot(ACF(mod.rat,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(mod.rat,maxLag=5,resType='n'),alpha=0.05)

resids<-as.vector(residuals(mod.rat,type='n'))
plot(yr,resids)
lines(lowess(yr,resids),lwd=2,col='red')
library(hier.part)

rich.part<-hier.part(rich,data.frame(logCa=logCa,slope=slope,north=northness,pdsi=PDSIavg,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu')
resid.part<-hier.part(resids,data.frame(logCa=logCa,slope=slope,north=northness,pdsi=PDSIavg,YrsOB=YrsOB,BP5Yrs=BP5Yrs,YrsSLB=YrsSLB),gof='Rsqu')

round(rich.part$IJ,3)
round(resid.part$IJ,3)

##and lastly
library(geoR)
x<-rep(1:10,each=10)
y<-rep(1:10,times=10)
z<-rnorm(100)+.25*x+.5*y
dat<-as.geodata(cbind(x,y,z))
plot(variog(dat,uv=1:50))
grid()

semi<-matrix(NA,ncol=6,nrow=20)
uni.name<-unique(name)
plot(NA,xlim=c(0,6),ylim=c(0,300),type='n')
for(i in 1:20){ 
 dat<-as.geodata(cbind(1:11,rep(0,11),rich[name==uni.name[i]]))
 semi[i,]<-variog(dat,max.dist=6,uv=1:6)$v
 points(1:6,semi[i,],type='l')
}

plot(apply(semi,2,mean))
lines(lowess(1:6,apply(semi,2,mean)))
avg.semi<-apply(semi,2,mean)
mod<-lm(avg.semi~c(1:6))
abline(mod)
##not sure how I would calculate residuals at this point


