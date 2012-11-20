spdat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
attach(endat)
##begin working only at the 100 m2 scale so no corners
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)

library(vegan)
library(nlme)

##prepare variables
spdat.n<-as.matrix(sqrt(spdat[,-1]))
colnames(spdat.n)

yr.f<-factor(yr)
site.m<-endatlv1[,8:26]
yr.m<-endatlv1[,28:37]
dist.m<-cbind(YrsOB,YrsSLB,BP5Yrs)

##big bluestem
spat.mods(spdat.n[,1],cbind(site.m,yr.m,dist.m),yr,name)
bb.mod<-gls(spdat.n[,1]~name+yr.f+dist.m,cor=corSpher(form=~yr|name,n=T))
summary(bb.mod)$t
##no strong resp
##little blue
spat.mods(spdat.n[,2],cbind(site.m,yr.m,dist.m),yr,name)
lb.mod<-gls(spdat.n[,2]~name+yr.f+dist.m,cor=corRatio(form=~yr|name,n=F))
round(summary(lb.mod)$t,3)
anova(lb.mod,type='m')
##increase with BP5
##dropseed
spat.mods(spdat.n[,3],cbind(site.m,yr.m,dist.m),yr,name)
spr.mod<-gls(spdat.n[,3]~name+yr.f+dist.m,cor=corRatio(form=~yr|name,n=F))
round(summary(spr.mod)$t,3)
anova(spr.mod,type='m')
##decrease with BP5
##indian grass
spat.mods(spdat.n[,5],cbind(site.m,yr.m,dist.m),yr,name)
in.mod<-gls(spdat.n[,5]~name+yr.f+dist.m,cor=corRatio(form=~yr|name,n=F))
round(summary(in.mod)$t,3)
anova(in.mod,type='m')
##no resp to dist
##switch grass
spat.mods(spdat.n[,7],cbind(site.m,yr.m,dist.m),yr,name)
sw.mod<-gls(spdat.n[,7]~name+yr.f+dist.m,cor=corRatio(form=~yr|name,n=F))
round(summary(sw.mod)$t,3)
anova(sw.mod,type='m')
##no resp to dist
##goldenrod
spat.mods(spdat.n[,6],cbind(site.m,yr.m,dist.m),yr,name)
go.mod<-gls(spdat.n[,6]~name+yr.f+dist.m,cor=corExp(form=~yr|name,n=F))
round(summary(go.mod)$t,3)
anova(go.mod,type='m')
##no resp to dist
##cumin ragweed (psilo)
spat.mods(spdat.n[,4],cbind(site.m,yr.m,dist.m),yr,name)
rw.mod<-gls(spdat.n[,4]~name+yr.f+dist.m,cor=corExp(form=~yr|name,n=F))
round(summary(rw.mod)$t,3)
anova(rw.mod,type='m')
##baptisia
spat.mods(spdat.n[,41],cbind(site.m,yr.m,dist.m),yr,name)
ba.mod<-gls(spdat.n[,41]~name+yr.f+dist.m,cor=corGaus(form=~yr|name,n=F))
round(summary(ba.mod)$t,3)
anova(ba.mod,type='m')


##abundance distribution across plots for each year
P<-matrix(NA,nrow=11,ncol=ncol(spdat.n))
uni.yrs<-1998:2008
for(i in 1:11){
 for(j in 1:ncol(spdat.n)){
  P[i,j]<-mean(spdat.n[yr==uni.yrs[i],j])
}}
rowsums<-apply(P,1,sum)
rel<-matrix(NA,nrow=11,ncol=ncol(spdat.n))
for(i in 1:11){
 for(j in 1:ncol(spdat.n)){
  rel[i,j]<-P[i,j]/rowsums[i]
}}
plot(rel[1,],type='n',ylim=c(0,.06))
for(i in 1:11){
 points(sort(rel[i,],decreasing=T),type='l',col=i)
}
##log P
plot(log(rel[1,]),type='n',ylim=log(c(.0001,.06)))
for(i in 1:11){
 points(log(sort(rel[i,],decreasing=T)),type='l',col=i)
}
legend('topright',as.character(c(1998:2008)),col=c(1:11),lty=1)


tapply(spdat.n,list(name,yr),mean)
