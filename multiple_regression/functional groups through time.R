endat<-read.table('C:/Users/dmcglinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
#endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
tapply(endat$rich,endat$lv,mean)
endatlv1<-endat[endat$lv==1&endat$cor==1,]
attach(endatlv1)

se<-function(x){sd(x)/sqrt(length(x))}

rich.avg<-tapply(rich,yr,mean)
rich.se<-tapply(rich,yr,se)
F.avg<-tapply(F.p,yr,mean)
F.se<-tapply(F.p,yr,se)
L.avg<-tapply(Leg.p,yr,mean)
L.se<-tapply(Leg.p,yr,se)
G3.avg<-tapply(G3.p,yr,mean)
G3.se<-tapply(G3.p,yr,se)
G4.avg<-tapply(G4.p,yr,mean)
G4.se<-tapply(G4.p,yr,se)
W.avg<-tapply(Wood.p,yr,mean)
W.se<-tapply(Wood.p,yr,se)

sr<-cbind(rich.avg,F.avg,L.avg,G3.avg,G4.avg,W.avg)
se<-cbind(rich.se,F.se,L.se,G3.se,G4.se,W.se)

##to get a sense of skewedness
boxplot(srs[,2]~yr,ylim=c(0,50),boxwex = 0.5)
for(i in 2:6){
 par(new=T)
 boxplot(srs[,i]~yr,ylim=c(0,50),boxwex = 0.5,axes=F,col=i)
}
##the funct groups with small means are a bit more skewed
##either pois estimates of mean should be derived or we should us bootstraped
##medians


##raw plots
cls<-c('black','black','grey40','grey60','grey80','black')
#plot(1998:2008,sr[,1],ylim=c(0,95),type='n',ylab='',xlab='',frame.plot=F)
plot(1998:2008,sr[,1],ylim=c(0,50),type='n',ylab='',xlab='',frame.plot=F,axes=F)
axis(side=1,cex.axis=1.5,lwd=3)
axis(side=2,cex.axis=1.5,lwd=3)
for(i in 2:6){
 points(1998:2008,sr[,i],type='l',col=cls[i],lty=i-1,lwd=3)
 arrows(1998:2008,sr[,i],1998:2008,sr[,i]+se[,i],angle=90,len=.05,col=cls[i],lwd=3)
 arrows(1998:2008,sr[,i],1998:2008,sr[,i]-se[,i],angle=90,len=.05,col=cls[i],lwd=3)
#95% confis
# arrows(1998:2008,sr[,i],1998:2008,sr[,i]+se[,i]*1.96,angle=90,len=.05,col=cls[i],lwd=3)
# arrows(1998:2008,sr[,i],1998:2008,sr[,i]-se[,i]*1.96,angle=90,len=.05,col=cls[i],lwd=3)
}
#legend('topleft',c('all species','forbs','legumes','C3','C4','W'),col=cls,lty=1:5,lwd=3,bty='n')
plot(1:10,1:10,type='n')
legend('center',c('Forbs','Legumes','C3 Graminoids','C4 Graminoids','Shrubs'),cex=1.5,col=cls[-1],lty=1:5,lwd=3,bty='n')

##raw plots in color
cls<-c('black','black','lightblue','lightgreen','orange','pink3')
#plot(1998:2008,sr[,1],ylim=c(0,95),type='n',ylab='',xlab='',frame.plot=F)
plot(1998:2008,sr[,1],ylim=c(0,50),type='n',ylab='',xlab='',frame.plot=F,axes=F)
axis(side=1,cex.axis=1.5,lwd=3)
axis(side=2,cex.axis=1.5,lwd=3)
for(i in 2:6){
 points(1998:2008,sr[,i],type='l',col=cls[i],lty=i-1,lwd=3)
 arrows(1998:2008,sr[,i],1998:2008,sr[,i]+se[,i],angle=90,len=.05,col=cls[i],lwd=3)
 arrows(1998:2008,sr[,i],1998:2008,sr[,i]-se[,i],angle=90,len=.05,col=cls[i],lwd=3)
#95% confis
# arrows(1998:2008,sr[,i],1998:2008,sr[,i]+se[,i]*1.96,angle=90,len=.05,col=cls[i],lwd=3)
# arrows(1998:2008,sr[,i],1998:2008,sr[,i]-se[,i]*1.96,angle=90,len=.05,col=cls[i],lwd=3)
}
#legend('topleft',c('all species','forbs','legumes','C3','C4','W'),col=cls,lty=1:5,lwd=3,bty='n')
plot(1:10,1:10,type='n')
legend('center',c('Forbs','Legumes','C3 Graminoids','C4 Graminoids','Shrubs'),cex=1.5,col=cls[-1],lty=1:5,lwd=3,bty='n')




##diff from mean
rich.avg<-tapply(rich-mean(rich),yr,mean)
rich.se<-tapply(rich-mean(rich),yr,se)
F.avg<-tapply(F.p-mean(F.p),yr,mean)
F.se<-tapply(F.p-mean(F.p),yr,se)
L.avg<-tapply(Leg.p-mean(Leg.p),yr,mean)
L.se<-tapply(Leg.p-mean(Leg.p),yr,se)
G3.avg<-tapply(G3.p-mean(G3.p),yr,mean)
G3.se<-tapply(G3.p-mean(G3.p),yr,se)
G4.avg<-tapply(G4.p-mean(G4.p),yr,mean)
G4.se<-tapply(G4.p-mean(G4.p),yr,se)

sr<-cbind(rich.avg,F.avg,L.avg,G3.avg,G4.avg)
se<-cbind(rich.se,F.se,L.se,G3.se,G4.se)
##r
plot(1998:2008,sr[,1],ylim=c(-20,20),type='n',ylab='',xlab='',frame.plot=F)
for(i in 1:5){
 points(1998:2008,sr[,i],type='l',col=i)
 arrows(1998:2008,sr[,i],1998:2008,sr[,i]+se[,i],angle=90,len=.1,col=i)
 arrows(1998:2008,sr[,i],1998:2008,sr[,i]-se[,i],angle=90,len=.1,col=i)
}

srs<-cbind(rich,F.p,Leg.p,G3.p,G4.p,Wood.p)

grp.rda<-rda(srs[,-1]~dist.m+Condition(name)+Condition(yr.f),scale=F)
plot(grp.rda,display=c("sp","cn"),scaling=0)



Y<-matrix(0,nrow=6,ncol=3)
rownames(Y)<-c('all','F','Leg','G3','G4','W')
colnames(Y)<-c('int','betayr','p-val')
for(i in 1:6){
 mod<-lm(srs[,i]~scale(yr,scale=F))
 Y[i,1:2]<-coef(mod)
 Y[i,3]<-anova(mod)$Pr[1]
}
round(Y,4)
        int betayr  p-val
all 76.2364 0.9332 0.0006
F   37.0364 0.2995 0.1059
Leg 10.2364 0.3109 0.0000
G3  15.8136 0.1377 0.1283
G4  11.7318 0.1600 0.0012
W    1.4227 0.0264 0.1630


library(nlme)

cyr<-scale(yr,scale=F)
best<-matrix(0,ncol=5,nrow=10)
for(i in 1:5){
 sr<-srs[,i]
 mod.exp<-gls(sr~cyr,cor=corExp(form=~yr|name))
 mod.gau<-gls(sr~cyr,cor=corGaus(form=~yr|name))
 mod.lin<-gls(sr~cyr,cor=corLin(form=~yr|name))
 mod.rat<-gls(sr~cyr,cor=corRatio(form=~yr|name))
 mod.sph<-gls(sr~cyr,cor=corSpher(form=~yr|name))
 mod.expn<-gls(sr~cyr,cor=corExp(form=~yr|name,nugget=T))
 mod.gaun<-gls(sr~cyr,cor=corGaus(form=~yr|name,nugget=T))
 mod.linn<-gls(sr~cyr,cor=corLin(form=~yr|name,nugget=F))
 mod.ratn<-gls(sr~cyr,cor=corRatio(form=~yr|name,nugget=T))
 mod.sphn<-gls(sr~cyr,cor=corSpher(form=~yr|name,nugget=T))
 best[,i]<-c(AIC(mod.exp),AIC(mod.gau),AIC(mod.lin),AIC(mod.rat),AIC(mod.sph),
             AIC(mod.expn),AIC(mod.gaun),AIC(mod.linn),AIC(mod.ratn),AIC(mod.sphn))
 best[,i]<-best[,i]-min(best[,i])
 #c(1:10)[best[,i]==min(best[,1])] 
}
rownames(best)<-c('exp','gau','lin','ratio','sph','expN','gauN','lin','ratioN','sphN')
colnames(best)<-c('all','F','Leg','G3','G4')
best
expN    1.847959  0.4189526   0.05820306   0.0000000  0.000000
gauN    1.890518  1.5569160   1.24780573   6.9953014 11.049796
lin    71.860784 59.0251505  47.39474806  34.2003320 32.647499
ratioN  0.000000  0.0000000   0.35023239   5.9142413  4.089250
sphN    2.032776  0.7608165   0.00000000   0.5600589  2.404596

Y<-matrix(0,nrow=5,ncol=4)
rownames(Y)<-c('all','F','Leg','G3','G4')
colnames(Y)<-c('int','betayr','p-val')
for(i in 1:5){
 sr<-srs[,i]
 if(i<=2){
  mod<-gls(sr~cyr,cor=corRatio(form=~yr|name,nugget=T))
 }
 else{
  if(i==3){
   mod<-gls(sr~cyr,cor=corSpher(form=~yr|name,nugget=T))
  } 
  else{
   mod<-gls(sr~cyr,cor=corExp(form=~yr|name,nugget=T))
 }} 
 Y[i,1:2]<-coef(mod)
 Y[i,3]<-mod$model$corStruct[1]
 Y[i,4]<-anova(mod)$p[2]
}
round(Y,4)
       [,1]   [,2]   [,3]   [,4]
all 74.8475 1.2505 1.8650 0.0002
F   36.2747 0.4398 1.8389 0.0483
Leg 10.0808 0.3096 4.2787 0.0000
G3  15.9234 0.1781 3.3985 0.0197
G4  11.5231 0.1782 2.5393 0.0014

mod.all<-gls(rich~1,cor=corExp(form=~yr|name))
plot(ACF(mod.all),alpha=0.05)
mod.F<-gls(F.p~1,cor=corExp(form=~yr|name))
mod.L<-gls(Leg.p~1,cor=corExp(form=~yr|name))
mod.G3<-gls(G3.p~1,cor=corExp(form=~yr|name))
mod.G4<-gls(G4.p~1,cor=corExp(form=~yr|name))
plot(ACF(mod.F),alpha=0.05)
plot(ACF(mod.L),alpha=0.05)
plot(ACF(mod.G3),alpha=0.05)
plot(ACF(mod.G4),alpha=0.05)
anova(gls(Leg.p~1),mod.L)

best<-matrix(0,ncol=5,nrow=10)
for(i in 1:5){
 sr<-srs[,i]
 mod.exp<-gls(sr~1,cor=corExp(form=~yr|name))
 mod.gau<-gls(sr~1,cor=corGaus(form=~yr|name))
 mod.lin<-gls(sr~1,cor=corLin(form=~yr|name))
 mod.rat<-gls(sr~1,cor=corRatio(form=~yr|name))
 mod.sph<-gls(sr~1,cor=corSpher(form=~yr|name))
 mod.expn<-gls(sr~1,cor=corExp(form=~yr|name,nugget=T))
 mod.gaun<-gls(sr~1,cor=corGaus(form=~yr|name,nugget=T))
 mod.linn<-gls(sr~1,cor=corLin(form=~yr|name,nugget=F))
 mod.ratn<-gls(sr~1,cor=corRatio(form=~yr|name,nugget=T))
 mod.sphn<-gls(sr~1,cor=corSpher(form=~yr|name,nugget=T))
 best[,i]<-c(AIC(mod.exp),AIC(mod.gau),AIC(mod.lin),AIC(mod.rat),AIC(mod.sph),
             AIC(mod.expn),AIC(mod.gaun),AIC(mod.linn),AIC(mod.ratn),AIC(mod.sphn))
 best[,i]<-best[,i]-min(best[,i])
 c(1:10)[best[,i]==min(best[,1])] 

}
rownames(best)<-c('exp','gau','lin','ratio','sph','expN','gauN','lin','ratioN','sphN')
colnames(best)<-c('all','F','Leg','G3','G4')
best
              all           F         Leg          G3       G4
expN    2.1703299  0.09600435   5.7105564   0.0000000  0.00000
gauN    0.4817378  1.02345833   0.0000000   5.2406280 10.75636
ratioN  0.0000000  0.00000000   0.6336475   4.6040278  3.54438

mod<-gls(F.p~1,cor=corRatio(form=~yr|name,nugget=T))
plot(Variogram(mod))
plot(Variogram(mod,resType='n'))
mod<-gls(Leg.p~1,cor=corGaus(form=~yr|name,nugget=T))
plot(Variogram(mod))
plot(Variogram(mod,resType='n'))
mod<-gls(G4.p~1,cor=corExp(form=~yr|name,nugget=T))
plot(Variogram(mod))
plot(Variogram(mod,resType='n'))


B<-matrix(0,nrow=5,ncol=3)
rownames(B)<-c('all','F','Leg','G3','G4')
for(i in 1:5){
 mod<-lm(srs[,i]~YrsOB)
 B[i,1:2]<-coef(mod)
 B[i,3]<-anova(mod)$Pr[1]
}
round(B,4)
       [,1]    [,2]   [,3]
all 72.9016  0.9863 0.0000
F   35.1772  0.5499 0.0001
Leg  8.6663  0.4644 0.0000
G3  15.5614  0.0746 0.2711
G4  12.2427 -0.1511 0.0000

