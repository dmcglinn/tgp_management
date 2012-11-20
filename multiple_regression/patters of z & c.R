endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
library(nlme)
library(car)
names(endat)
attach(endat)
endatlv2<-endat[lv==2,]
detach(endat)
attach(endatlv2)

plot(c,z)
lines(lowess(c,z))
identify(c,z)
identify(c,z,labels=name)

pairs(rich1~rich5+c+z,panel=panel.smooth,lwd=2)
pairs(rich1~rich2+rich3+rich4+rich5,panel=panel.smooth,lwd=2)

mod0.lm<-lm(z~c)
summary(mod0.lm)
mod1.lm<-aov(z~Error(name/cor))
summary(mod1.lm)
mod2.lm<-aov(z~c+Error(name/cor))
summary(mod2.lm)

zdata<-groupedData(z~c|name/cor,data=endatlv2)
mod1.ls<-lmList(z~1|name/cor,data=zdata)
plot(mod1.ls)
plot(intervals(mod1.ls))

mod1.lme<-lme(z~1,random=~1|name/cor,method="ML")
summary(mod1.lme)
Random effects:
 Formula: ~1 | name
        (Intercept)
StdDev:  0.01451821

 Formula: ~1 | cor %in% name
        (Intercept)   Residual
StdDev:  0.01105779 0.03569639
mod2.lme<-lme(z~c,random=~1|name/cor,method="ML")
anova(mod1.lme,mod2.lme)
summary(mod2.lme)
mod3.lme<-update(mod2.lme,random=~1|name)
anova(mod2.lme,mod3.lme)
##so you can drop the corners as a random effect
plot(mod3.lme)
##examining patters
plot(c,z,col=1,pch=19)
lines(c(min(c),max(c)),c(min(c)*fixef(mod3.lme)[2]+fixef(mod3.lme)[1],max(c)*fixef(mod3.lme)[2]+fixef(mod3.lme)[1]),type='l',col='red',lwd=2)
##just for plots
uni.name<-rownames(ranef(mod3.lme))
for(i in 1:20){
 lines(c(min(c[name==uni.name[i]]),max(c[name==uni.name[i]])),
       c(min(c[name==uni.name[i]])*fixef(mod3.lme)[2]+fixef(mod3.lme)[1]+ranef(mod3.lme,lev=1)[i,1],
         max(c[name==uni.name[i]])*fixef(mod3.lme)[2]+fixef(mod3.lme)[1]+ranef(mod3.lme,lev=1)[i,1]),
        type='l',col='grey',lwd=2)
 lines(lowess(c(c[name==uni.name[i]],c[name==uni.name[i]]),
       c(z[name==uni.name[i]],z[name==uni.name[i]])),
        type='l',col='skyblue',lwd=2)
}
legend("bottomleft",c('fixed effect','plot ranef','lowess'),col=c("red","grey",'skyblue'),lty=1,lwd=2)

##drop three outliers, #336,340,731
z[c(336,340,731)]
c[c(336,340,731)]

mod1.lme<-lme(z~1,random=~1|name/cor,method="ML",subset=c(-336,-340,-731))
mod2.lme<-lme(z~c,random=~1|name,method="ML",subset=c(-336,-340,-731))
mod3.lme<-lme(z~c,random=~1|name/cor,method="ML",subset=c(-336,-340,-731))
anova(mod1.lme,mod2.lme,mod3.lme)
#no need for corner random effect
##graphics
c.new<-c[c(-336,-340,-731)]
z.new<-z[c(-336,-340,-731)]
z[1663:1665]
z.new[1663:1665]
c[1663:1665]
c.new[1663:1665]

##subseting seems to have worked
plot(c.new,z.new,col=1,pch=19)
lines(c(min(c.new),max(c.new)),c(min(c.new)*fixef(mod2.lme)[2]+fixef(mod2.lme)[1],max(c)*fixef(mod2.lme)[2]+fixef(mod2.lme)[1]),type='l',col='red',lwd=3)
#lines(lowess(c.new,z.new),type='l',col='skyblue',lwd=3)
##just for plots
for(i in 1:20){
 lines(c(min(c.new[name==uni.name[i]]),max(c.new[name==uni.name[i]])),
       c(min(c.new[name==uni.name[i]])*fixef(mod2.lme)[2]+fixef(mod2.lme)[1]+ranef(mod2.lme,lev=1)[i,1],
         max(c.new[name==uni.name[i]])*fixef(mod2.lme)[2]+fixef(mod2.lme)[1]+ranef(mod2.lme,lev=1)[i,1]),
        type='l',col='grey',lwd=2)
}
for(i in 1:20){
 lines(lowess(c(c.new[name==uni.name[i]],c.new[name==uni.name[i]]),
       c(z.new[name==uni.name[i]],z.new[name==uni.name[i]])),
        type='l',col='skyblue',lwd=2)
}
legend("bottomleft",c('fixed effect','plot ranef','lowess'),col=c("red","grey",'skyblue'),lty=1,lwd=2)

