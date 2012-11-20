
dat<-read.table('C:/Users/hurlbertlab/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/rich98-09.csv',sep=',',header=T)
attach(dat)
names(dat)
[1] "plot" "year" "L1"

plot(year,L1)
lines(lowess(year,L1),col='red',lwd=2)
points(unique(year),tapply(L1,year,mean),col='blue',pch=19,type='o')

mod<-lm(L1~as.factor(plot))
plot(year,residuals(mod))
lines(lowess(year,residuals(mod)),col='red',lwd=2)

mod<-lm(L1~as.factor(plot)+year+I(year^2))
summary(mod)
anova(mod)

diag.func<-function(x,grp){
 plot(unique(grp),tapply(x,grp,mean))
 arrows(unique(grp),tapply(x,grp,mean),unique(grp),tapply(x,grp,mean)+tapply(x,grp,function(x){sd(x)/sqrt(12)}),angle = 90)
 arrows(unique(grp),tapply(x,grp,mean),unique(grp),tapply(x,grp,mean)-tapply(x,grp,function(x){sd(x)/sqrt(12)}),angle = 90)
}

diag.func(L1,year)

boxplot(L1~plot)
points(1:20,L1[year==2009],col='red')
