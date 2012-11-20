##this file has just 98 to 08 data 
endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)

##this file has 95 to 08 just the clim data though
clim.var<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/climate(rain&temp).csv',sep=',',header=T)
 [1] "yr"           "rain1"        "rain2"        "rain3"        "rain4"       
 [6] "rain5"        "rain6"        "rain7"        "rain8"        "rain9"       
[11] "rain10"       "rain11"       "rain12"       "temp1"        "temp2"       
[16] "temp3"        "temp4"        "temp5"        "temp6"        "temp7"       
[21] "temp8"        "temp9"        "temp10"       "temp11"       "temp12"      
[26] "sum.rain.tot" "win.rain.tot" "spr.rain.tot" "sum.rain.avg" "win.rain.avg"
[31] "spr.rain.avg" "sum.temp"     "win.temp"     "spr.temp"     "rain.tot"    
[36] "rain.avg"     "temp.avg"    

attach(clim.var)
tot.rain.cen<-(rain.tot-mean(rain.tot))*10 ##center and convert to mm of rain
yrs.all<-1995:2008
detach(clim.var)

attach(endat)

yrs.samp<- 1998:2008
rich1.mat<-as.matrix(tapply(rich1,list(plot,yr),mean))
rich1.plcen<-rich1.mat-apply(rich1.mat,1,mean) ##centered on plot means

G3.mat<-as.matrix(tapply(G3.p,list(plot,yr),mean))
G3.plcen<-G3.mat-apply(G3.mat,1,mean)
G4.mat<-as.matrix(tapply(G4.p,list(plot,yr),mean))
G4.plcen<-G4.mat-apply(G4.mat,1,mean)
Gsum.mat<-as.matrix(tapply(Gsum.p,list(plot,yr),mean))
Gsum.plcen<-Gsum.mat-apply(Gsum.mat,1,mean)
F.mat<-as.matrix(tapply(F.p,list(plot,yr),mean))
F.plcen<-F.mat-apply(F.mat,1,mean)
Leg.mat<-as.matrix(tapply(Leg.p,list(plot,yr),mean))
Leg.plcen<-Leg.mat-apply(Leg.mat,1,mean)
Wood.mat<-as.matrix(tapply(Wood.p,list(plot,yr),mean))
Wood.plcen<-Wood.mat-apply(Wood.mat,1,mean)
##annual, perrenial,tree, and shrub
A.mat<-as.matrix(tapply(A.p,list(plot,yr),mean))
A.plcen<-A.mat-apply(A.mat,1,mean)
P.mat<-as.matrix(tapply(P.p,list(plot,yr),mean))
P.plcen<-P.mat-apply(P.mat,1,mean)
T.mat<-as.matrix(tapply(T.p,list(plot,yr),mean))
T.plcen<-T.mat-apply(T.mat,1,mean)
S.mat<-as.matrix(tapply(S.p,list(plot,yr),mean))
S.plcen<-S.mat-apply(S.mat,1,mean)



###overlay plot ###  ##centered on plot means##
###with groups: c3,c4,forb,ect as subgroups##
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
axis(1,cex.axis=1.5)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,apply(rich1.plcen,2,mean),axes=F,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(-25,10),pch=16,col=1,frame.plot=T,lwd=2,cex.axis=1.5)
axis(2,at=seq(-10,10,5),label=T,cex.axis=1.5)
points(yrs.samp,apply(G3.plcen,2,mean),type='l',col=2,lwd=2)
points(yrs.samp,apply(G4.plcen,2,mean),type='l',col=3,lwd=2)
points(yrs.samp,apply(Gsum.plcen,2,mean),type='l',col=4,lwd=2)
points(yrs.samp,apply(F.plcen,2,mean),type='l',col=5,lwd=2)
points(yrs.samp,apply(Leg.plcen,2,mean),type='l',col=6,lwd=2)
points(yrs.samp,apply(Wood.plcen,2,mean),type='l',col=7,lwd=2)
legend('topleft',c("all sr","C3 grass sr","C4 grass sr","all grass sr","Forb (no Legume) sr","Legume sr","Woody sr"),col=1:7,lty=1,lwd=2,bty='n')
#grid()

###
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
axis(1,cex.axis=1.5)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,apply(rich1.plcen,2,mean),axes=F,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(-25,10),pch=16,col=1,frame.plot=T,lwd=2,cex.axis=1.5)
axis(2,at=seq(-10,10,5),label=T,cex.axis=1.5)
points(yrs.samp,apply(A.plcen,2,mean),type='l',col=2,lwd=2)
points(yrs.samp,apply(P.plcen,2,mean),type='l',col=3,lwd=2)
points(yrs.samp,apply(T.plcen,2,mean),type='l',col=4,lwd=2)
points(yrs.samp,apply(S.plcen,2,mean),type='l',col=5,lwd=2)
legend('topleft',c("all sr","A","P","T","S"),col=1:5,lty=1,lwd=2,bty='n')
#grid()


##
sr1.yr<-by(rich1,yr,mean)
sr2.yr<-by(rich2,yr,mean)
sr3.yr<-by(rich3,yr,mean)
sr4.yr<-by(rich4,yr,mean)
sr5.yr<-by(rich5,yr,mean)

G3.yr<-by(G3.p,yr,mean)
G4.yr<-by(G4.p,yr,mean)
Gsum.yr<-by(Gsum.p,yr,mean)
F.yr<-by(F.p,yr,mean)
Leg.yr<-by(Leg.p,yr,mean)
Wood.yr<-by(Wood.p,yr,mean)

G3.cov.yr<-by(G3.cov,yr,mean)
G4.cov.yr<-by(G4.cov,yr,mean)
Gsum.cov.yr<-by(Gsum.cov,yr,mean)
F.cov.yr<-by(Forb.cov,yr,mean)
Leg.cov.yr<-by(Leg.cov,yr,mean)
Wood.cov.yr<-by(Wood.cov,yr,mean)

A.yr<-by(A.p,yr,mean) ##annual
P.yr<-by(P.p,yr,mean) ##perrenial
S.yr<-by(S.p,yr,mean) ##shrub
T.yr<-by(T.p,yr,mean) ##tree


yrs.samp<-1998:2008

detach(endat)
plot(yrs.all,tot.rain.cen,type='h')


###overlay plot ###
###at each level##
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
#add text to right axis
#text(rep(2008+1,5),seq(-500,500,length.out=3),labels=c('-500','0','500'),srt=90,cex=2)
#allow text function to write to region outside graph
#par(xpd=TRUE)
#text(2008+1,35,'Total Yearly Rainfall (in)',srt=90)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,sr1.yr,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(0,90),pch=16,col=1,frame.plot=T,lwd=2,cex.axis=1.5)
points(yrs.samp,sr2.yr,type='l',col=2,lwd=2)
points(yrs.samp,sr3.yr,type='l',col=3,lwd=2)
points(yrs.samp,sr4.yr,type='l',col=4,lwd=2)
points(yrs.samp,sr5.yr,type='l',col=5,lwd=2)
legend('topleft',legend=c(100,10,1,.1,.01),col=1:5,lty=1,lwd=2,bty='n')

###overlay plot ###
###with groups: c3,c4,forb,ect as subgroups##
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
#add text to right axis
#text(rep(2008+1,5),seq(-500,500,length.out=3),labels=c('-500','0','500'),srt=90,cex=2)
#allow text function to write to region outside graph
#par(xpd=TRUE)
#text(2008+1,35,'Total Yearly Rainfall (in)',srt=90)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,sr1.yr,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(0,90),pch=16,col=1,frame.plot=T,lwd=2,cex.axis=1.5)
points(yrs.samp,G3.yr,type='l',col=2,lwd=2)
points(yrs.samp,G4.yr,type='l',col=3,lwd=2)
points(yrs.samp,Gsum.yr,type='l',col=4,lwd=2)
points(yrs.samp,F.yr,type='l',col=5,lwd=2)
points(yrs.samp,Leg.yr,type='l',col=6,lwd=2)
points(yrs.samp,Wood.yr,type='l',col=7,lwd=2)
legend('topleft',c("all sr","C3 grass sr","C4 grass sr","all grass sr","Forb (no Legume) sr","Legume sr","Woody sr"),col=1:7,lty=1,lwd=2,bty='n')
#grid()

###overlay plot ###
###with groups: c3,c4,forb,ect as subgroups##
###on cover of the groups ###
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
#add text to right axis
#text(rep(2008+1,5),seq(-500,500,length.out=3),labels=c('-500','0','500'),srt=90,cex=2)
#allow text function to write to region outside graph
#par(xpd=TRUE)
#text(2008+1,35,'Total Yearly Rainfall (in)',srt=90)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,sr1.yr,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(0,90),pch=16,col=1,frame.plot=T,lwd=2,cex.axis=1.5)
#points(yrs.samp,apply(cbind(F.yr,Gsum.yr,Leg.yr,Wood.yr),1,sum)+1,col='pink',type='l',lwd=2)
points(yrs.samp,G3.cov.yr,type='l',col=2,lwd=2)
points(yrs.samp,G4.cov.yr,type='l',col=3,lwd=2)
points(yrs.samp,Gsum.cov.yr,type='l',col=4,lwd=2)
points(yrs.samp,F.cov.yr,type='l',col=5,lwd=2)
points(yrs.samp,Leg.cov.yr,type='l',col=6,lwd=2)
points(yrs.samp,Wood.cov.yr,type='l',col=7,lwd=2)
legend('topleft',c("all sr","C3 grass cov","C4 grass cov","all grass cov","Forb (no Legume) cov","Legume cov","Woody cov"),col=1:7,lty=1,lwd=2,bty='n')
#grid()

###overlay plot ###
###with groups: annual, perrenial,ect as subgroups##
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
#add text to right axis
#text(rep(2008+1,5),seq(-500,500,length.out=3),labels=c('-500','0','500'),srt=90,cex=2)
#allow text function to write to region outside graph
#par(xpd=TRUE)
#text(2008+1,35,'Total Yearly Rainfall (in)',srt=90)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,sr1.yr,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(0,90),pch=16,col=1,frame.plot=T,lwd=2,cex.axis=1.5)
#points(yrs.samp,apply(cbind(F.yr,Gsum.yr,Leg.yr,Wood.yr),1,sum)+1,col='pink',type='l',lwd=2)
points(yrs.samp,A.yr,type='l',col=2,lwd=2)
points(yrs.samp,P.yr,type='l',col=3,lwd=2)
points(yrs.samp,S.yr,type='l',col=4,lwd=2)
points(yrs.samp,T.yr,type='l',col=5,lwd=2)
legend('topleft',c("all sr","annual sr","perennial sr","shrub sr","tree sr"),col=1:7,lty=1,lwd=2,bty='n')
#grid()


###overlay plot ###
###without groups###
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
#add text to right axis
#text(rep(2008+1,5),seq(-500,500,length.out=3),labels=c('-500','0','500'),srt=90,cex=2)
#allow text function to write to region outside graph
#par(xpd=TRUE)
#text(2008+1,35,'Total Yearly Rainfall (in)',srt=90)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,sr1.yr,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(50,90),pch=16,col=1,frame.plot=T,lwd=2,cex.axis=1.5)
#grid()

##add bootstrap confid intervals for year means

##create endatlv1 file
endatlv1<-endat[endat$lv==1&endat$cor==1,]
dim(endatlv1)
##get a sense of the skewedness of sr within a year
boxplot(endatlv1$rich~endatlv1$yr)
library(boot)

yr.f<-as.factor(endatlv1$yr)
yr.mod<-lm(endatlv1$rich~-1+yr.f)

fit<-yr.mod$fitted
e<-residuals(yr.mod)
e.c<-e-mean(e)
X<-yr.f
boot.coef<- function(data,indices){
 y <- fit + e.c[indices]
 mod <- lm(y ~-1+ X)
 coef(mod)
}

yr.boot<-boot(endatlv1,boot.coef,4999,strata=endatlv1$yr);alarm()
plot(yr.boot,index=3)
qq.plot(yr.boot$t[,2])

CI<-function(boot.obj){
 types<-c("norm","basic","perc")
 conf<-array(NA,dim=c(11,3,2))##11 = # of yrs, 3=# of cis, 2 =1 high+1 low
 for(i in 1:dim(conf)[1]){
  for(j in 1:dim(conf)[2]){
   if(j == 1){
    ci<-boot.ci(boot.obj,index=i,type="norm")
    conf[i,j,]<-ci$norm[-1]
   }   
   if(j == 2){
    ci<-boot.ci(boot.obj,index=i,type="basic")
    conf[i,j,]<-ci$basic[-(1:3)]
   }
   if(j == 3){
    ci<-boot.ci(boot.obj,index=i,type="perc")
    conf[i,j,]<-ci$perc[-(1:3)]
 }}}
conf}

##all boot output
output<-matrix(NA,ncol=11,nrow=4999)
for(i in 1:11){
 output[,i]<-yr.boot$t[,i]
}
yrmat<-matrix(NA,ncol=11,nrow=4999)
for(i in 1:11){ 
 yrmat[,i]<-yrs.samp[i]
}

##boxplot along with yr avgs
plot(yrs.samp,sr1.yr,pch=19,lwd=2,type='n',ylim=c(55,95),xlab='',ylab='',frame.plot=F,cex.axis=1.5)
boxplot(output~yrmat,ylim=c(55,95),axes=F,xlim=c(2,11),add=T,at=1998:2008,boxwex=.5)

##calcualte confidence interval
yr.ci<-CI(yr.boot)
attach(endatlv1)
yr.se<-tapply(rich,yr,sd)/sqrt(20)
detach(endatlv1)

plot(yrs.samp,sr1.yr,pch=19,ylim=c(55,95))
points(cbind(yrs.samp,yrs.samp,yrs.samp),yr.ci[,,1],pch=1,cex=1,col=2:4)
points(cbind(yrs.samp,yrs.samp,yrs.samp),yr.ci[,,2],pch=1,cex=1,col=2:4)
legend('bottomright',c('norm','basic','perc'),col=2:4,pch=1)

##95% CI for 4999 bootstrap runs
plot(yrs.samp,sr1.yr,pch=19,ylim=c(60,90),ylab='',xlab='',frame.plot=F,type='o',lwd=2)
arrows(yrs.samp,sr1.yr,yrs.samp,yr.ci[,3,1],angle=90,len=.1)
arrows(yrs.samp,sr1.yr,yrs.samp,yr.ci[,3,2],angle=90,len=.1)

#####With frames
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10,frame.plot=F)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,sr1.yr,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(45,90),pch=16,col=1,lwd=2,axes=F,frame.plot=T)
axis(side=1,cex.axis=1.5)
axis(side=2,at=seq(60,90,length.out=4),lablel=T,cex.axis=1.5)
#arrows(yrs.samp,sr1.yr,yrs.samp,yr.ci[,3,1],angle=90,len=.05,lwd=2)
#arrows(yrs.samp,sr1.yr,yrs.samp,yr.ci[,3,2],angle=90,len=.05,lwd=2)
arrows(yrs.samp,sr1.yr,yrs.samp,sr1.yr+yr.se*1.96,angle=90,len=.05,lwd=2,col='black')
arrows(yrs.samp,sr1.yr,yrs.samp,sr1.yr-yr.se*1.96,angle=90,len=.05,lwd=2,col='black')


######### Without frames
plot(yrs.all,tot.rain.cen,type='h',axes=F,ylab='',xlab='',xlim=c(1995,2008),ylim=c(-500,2000),lwd=10,frame.plot=F)
#set right axis
axis(4,at=seq(-500,500,length.out=3),label=T,cex.axis=1.5)
#allow overlaying of next plot
par(new=TRUE)
plot(yrs.samp,sr1.yr,type='l',ylab="",xlab="",xlim=c(1995,2008),ylim=c(45,90),pch=16,col=1,lwd=2,axes=F,frame.plot=F)
axis(side=1,cex.axis=1.5)
axis(side=2,at=seq(60,90,length.out=4),lablel=T,cex.axis=1.5)
arrows(yrs.samp,sr1.yr,yrs.samp,yr.ci[,3,1],angle=90,len=.05,lwd=2)
arrows(yrs.samp,sr1.yr,yrs.samp,yr.ci[,3,2],angle=90,len=.05,lwd=2)

#grid()





