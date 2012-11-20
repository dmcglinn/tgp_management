endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
spdat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
sp<-spdat[,-1]
library(nlme)
library(car)
names(endat)
attach(endat)

endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)

yr.f<-as.factor(yr)
mod<-lm(rich~name+yr.f+bison+YrsSLB+BP5Yrs)
mod<-lm(rich~name+yr.f+YrsOB+YrsSLB+BP5Yrs)

summary(mod)
par(mfrow=c(2,3))
termplot(mod,partial.resid=T)

ca.coef<-matrix(NA,nrow=11,ncol=2)
uni.yr<-1998:2008
for(i in 1:11){
 mod<-lm(rich[yr==uni.yr[i]]~logCa[yr==uni.yr[i]])
 ca.coef[i,1]<-coef(mod)[2]
 ca.coef[i,2]<-summary(mod)$coef[2,2]
}

slope.coef<-matrix(NA,nrow=11,ncol=2)
uni.yr<-1998:2008
for(i in 1:11){
 mod<-lm(rich[yr==uni.yr[i]]~slope[yr==uni.yr[i]])
 slope.coef[i,1]<-coef(mod)[2]
 slope.coef[i,2]<-summary(mod)$coef[2,2]
}

plot(1998:2008,ca.coef[,1],type='o',ylim=c(-100,20),ylab='slope of richness on logCa',xlab='year')
abline(h=0)
arrows(1998:2008,ca.coef[,1],1998:2008,ca.coef[,1]+1.96*ca.coef[,2],angle=90,len=.1)
arrows(1998:2008,ca.coef[,1],1998:2008,ca.coef[,1]-1.96*ca.coef[,2],angle=90,len=.1)
par(new=T)
#plot(1998:2008,tapply(rich,yr.f,mean),axes=F,xlab='',ylab='',type='o',pch=19)
plot(1998:2008,tapply(logCa,yr.f,mean),axes=F,xlab='',ylab='',type='o',pch=19)
axis(side=4)

ca.coef<-rep(NA,11)
uni.yr<-1998:2008
for(i in 1:11){
 ca.std<-scale(logCa[yr==uni.yr[i]])
 ca.coef[i]<-coef(lm(rich[yr==uni.yr[i]]~ca.std))[2]
}

plot(1998:2008,ca.coef,type='o')

##cumulative richness
crich<-rep(NA,20)
uni.name<-unique(name)
for(i in 1:20){
 crich[i]<-sum(ifelse(apply(sp[name==uni.name[i],],2,sum)>0,1,0))
}
avg.ca<-tapply(logCa,name,mean)

avg.rich<-tapply(rich,name,mean)

##Palmer et al. 2003 Folia Fig4
plot(avg.ca,avg.rich,ylim=c(45,160),ylab='richness',xlab='average log Ca')
mod<-lm(avg.rich~avg.ca)
newx<-data.frame(avg.ca=seq(min(avg.ca),max(avg.ca),length.out=3))
points(newx$avg.ca,predict(mod,newdata=newx),type='l')
points(avg.ca,crich,pch=19)
mod<-lm(crich~avg.ca)
points(newx$avg.ca,predict(mod,newdata=newx),type='l')
legend('bottomleft',c('cumulative richness','average richness'),pch=c(19,1),lty=1,bty='n')


plot(avg.ca,crich)
abline(lm(crich~avg.ca))
summary(lm(crich~avg.ca))
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)   340.07      62.61   5.432 3.68e-05 ***
avg.ca        -62.98      18.31  -3.440  0.00292 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

Residual standard error: 13.66 on 18 degrees of freedom
Multiple R-squared: 0.3967,     Adjusted R-squared: 0.3632

rich98<-rich[yr==1998]
srich<-crich-rich98

num.burns<-tapply(burn,name,sum)
par(mfrow=c(1,2))
plot(num.burns,srich)
lines(lowess(num.burns,srich))
num.burns<-tapply(burn,name,sum) + BP5Yrs[yr==1998]
plot(num.burns,srich)
lines(lowess(num.burns,srich))

srich.98<-crich-rich98
srich.avg<-(crich)-(avg.rich)
srich.lavg<-log10(crich)-log10(avg.rich)
mod<-lm(srich~num.burns)
mod2<-update(mod,.~.+I(num.burns^2))
anova(mod,mod2)
summary(mod)
abline(mod)

yrs.sd<-tapply(YrsSLB,name,sd)
plot(yrs.sd,srich)
lines(lowess(yrs.sd,srich))
mod2<-lm(srich~yrs.sd+I(yrs.sd^2))
newx<-data.frame(yrs.sd=seq(min(yrs.sd),max(yrs.sd),length.out=100))
points(newx$yrs.sd,predict(mod2,newdata=newx),col='red',type='l')
#Multiple R-squared: 0.3365,     Adjusted R-squared: 0.2584 

yrs.avg<-tapply(YrsSLB,name,mean)
bis.avg<-tapply(bison,name,mean)
bis.avg<-tapply(YrsOB,name,mean)
plot(yrs.avg,srich,pch=19,ylab='',xlab='')
#points(yrs.avg,srich,cex=sqrt(bis.avg))
lines(lowess(yrs.avg,srich))
mod2<-lm(srich~yrs.avg+I(yrs.avg^2))
newx<-data.frame(yrs.avg=seq(min(yrs.avg),max(yrs.avg),length.out=100))
points(newx$yrs.avg,predict(mod2,newdata=newx),col='red',type='l')
summary(mod2)
#Multiple R-squared: 0.3365,     Adjusted R-squared: 0.2584 
par(new=T)
bp5.avg<-tapply(BP5Yrs,name,mean)
plot(yrs.avg[order(yrs.avg)],bp5.avg[order(yrs.avg)],type='o',xlab='',ylab='',axes=F)
axis(side=4)

par(mfrow=c(1,3))
plot(yrs.avg,srich.98,ylab='',xlab='',pch=19)
lines(lowess(yrs.avg,srich.98))
mod<-lm(srich.98~yrs.avg)
mod2<-update(mod,.~.+I(yrs.avg^2))
anova(mod,mod2)
newx<-data.frame(yrs.avg=seq(min(yrs.avg),max(yrs.avg),length.out=100))
points(newx$yrs.avg,predict(mod2,newdata=newx),col='red',type='l')
plot(yrs.avg,srich.avg,ylab='',xlab='',pch=19)
lines(lowess(yrs.avg,srich.avg))
mod<-lm(srich.avg~yrs.avg)
mod2<-update(mod,.~.+I(yrs.avg^2))
anova(mod,mod2)
newx<-data.frame(yrs.avg=seq(min(yrs.avg),max(yrs.avg),length.out=100))
points(newx$yrs.avg,predict(mod2,newdata=newx),col='red',type='l')
plot(yrs.avg,srich.lavg,ylab='',xlab='',pch=19)
lines(lowess(yrs.avg,srich.lavg))
mod<-lm(srich.lavg~yrs.avg)
mod2<-update(mod,.~.+I(yrs.avg^2))
anova(mod,mod2)
newx<-data.frame(yrs.avg=seq(min(yrs.avg),max(yrs.avg),length.out=100))
points(newx$yrs.avg,predict(mod,newdata=newx),col='red',type='l')

##with bison plots visible
par(mfrow=c(1,1))
plot(yrs.avg,srich.avg,ylab='',xlab='',pch=19)
points(yrs.avg,srich.avg,cex=1+sqrt(bis.avg))
lines(lowess(yrs.avg,srich.avg))
mod<-lm(srich.avg~yrs.avg)
mod2<-update(mod,.~.+I(yrs.avg^2))
anova(mod,mod2)
newx<-data.frame(yrs.avg=seq(min(yrs.avg),max(yrs.avg),length.out=100))
points(newx$yrs.avg,predict(mod2,newdata=newx),col='red',type='l')



##
plot(bp5.avg,srich)
lines(lowess(bp5.avg,srich))

rich.avg<-tapply(rich,name,mean)
plot(yrs.avg,rich.avg)
lines(lowess(yrs.avg,rich.avg))


bp5.sd<-tapply(BP5Yrs,name,sd)
plot(bp5.sd,srich)
lines(lowess(bp5.sd,srich))

################################
################################
YrsOC<-rep(0,length(YrsOB))
for(i in 1:length(YrsOC)){
 if(bison[i]==0){
  YrsOC[i]<-yr[i]-1997
}}
summary(YrsOC)

par(mfrow=c(1,2))
mod.full<-lm(rich~name+YrsOB)
termplot(mod.full,partial.resid = T,se=T,term='YrsOB',ylim=c(-25,25))
mod.full<-lm(rich~name+YrsOC)
termplot(mod.full,partial.resid = T,se=T,term='YrsOC',ylim=c(-25,25))

par(mfrow=c(1,2))
termplot(mod.full,partial.resid = T,se=T,term='YrsOB',col.term=1,ylim=c(-25,25),ylab='',xlab='')
termplot(mod.full,partial.resid = T,se=T,term='YrsOC',ylim=c(-25,25))



par(mfrow=c(1,2))
mod.pl<-lm(rich~name)
resids<-resid(mod.pl)
mod.bis<-lm(resids~YrsOB)
mod.bis2<-update(mod.bis,.~.+I(YrsOB^2))
anova(mod.bis,mod.bis2)
plot(YrsOB,resids)
newx<-data.frame(YrsOB=seq(min(YrsOB),max(YrsOB),length.out=100))
preds<-predict(mod.bis,newdata=newx)
confid<-predict(mod.bis,newdata=newx,interval="confidence")
points(newx$YrsOB,preds,type='l')
lines(newx$YrsOB,confid[,2],lty=2)
lines(newx$YrsOB,confid[,3],lty=3)

mod.cat<-lm(resids~YrsOC)
plot(YrsOC,resids)
summary(mod.cat)

####
gls.b<-gls(resids~YrsOB)
gls.bar1<-update(gls.b,cor=corAR1(form=~yr|name))
gls.bar2<-update(gls.bar1,cor=corARMA(p=2))
anova(gls.b,gls.bar1,gls.bar2)
##marginally better on ar1
summary(gls.bar1)##still sig
summary(gls.bar2)
plot(ACF(gls.bar1,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(gls.bar1,maxLag=5,resType='n'),alpha=0.05)

gls.c<-gls(resids~YrsOC)
gls.car1<-update(gls.c,cor=corAR1(form=~yr|name))
gls.car2<-update(gls.car1,cor=corARMA(p=2))
anova(gls.c,gls.car1,gls.car2)
##marginally better on ar1
summary(gls.car1)##still sig
summary(gls.car2)
plot(ACF(gls.car1,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(gls.car1,maxLag=5,resType='n'),alpha=0.05)


###
mod.bis<-lm(resids~YrsOB)
mod.bis.n<-lm(resids~YrsOB,subset=bison==1)
summary(mod.bis)
plot(YrsOB[bison==1],resids[bison==1])
abline(mod.bis)
abline(mod.bis.n,lty=2)


##########
par(mfrow=c(1,2))
mod.bis<-lm(resids~YrsOB)
mod.bis2<-update(mod.bis,.~.+I(YrsOB^2))
anova(mod.bis,mod.bis2)
plot(YrsOB,resids)
newx<-data.frame(YrsOB=seq(min(YrsOB),max(YrsOB),length.out=100))
preds<-predict(mod.bis,newdata=newx)
confid<-predict(mod.bis,newdata=newx,interval="confidence")
points(newx$YrsOB,preds,type='l')
lines(newx$YrsOB,confid[,2],lty=2)
lines(newx$YrsOB,confid[,3],lty=3)

full.mod<-gls(rich~YrsOB+YrsSLB+BP5Yrs+name+yr.f,cor=corRatio(form=~yr|name,nugget=T),meth="ML")
resids<- rich - model.matrix(formula(full.mod))[,-2]%*%coef(full.mod)[-2]
resids.c<-resids-mean(resids)
spat.mods(resids.c,YrsOB,yr,name,meth="REML")
mod.bis<-gls(resids.c~YrsOB,cor=corAR1(form=~yr|name),meth="ML")
mod.bis<-lm(resids.c~YrsOB)
plot(Variogram(mod.bis))
plot(ACF(mod.bis,resType='n'),alpha=0.05)
plot(YrsOB,resids)
newx<-data.frame(YrsOB=seq(min(YrsOB),max(YrsOB),length.out=100))
preds<-predict(mod.bis,newdata=newx)
confid<-predict(mod.bis,newdata=newx,interval="confidence")
points(newx$YrsOB,preds,type='l')
lines(newx$YrsOB,confid[,2],lty=2)
lines(newx$YrsOB,confid[,3],lty=3)



par(mfrow=c(2,2));plot(lm(resids.c~YrsOB))
plot(ACF(mod,maxLag=5,resType='n'),alpha=0.05)
summary(mod)
fit<-predict(mod)
e<-residuals(mod)
e.c<- e-mean(e)
X <- YrsOB
