endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
names(endat)
attach(endat)
##work only with level 1
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)
yr.f<-factor(yr)

summary(lm(rich~yr.f+bison+YrsSLB+BP5Yrs))
summary(lm(rich~name+yr.f+bison+YrsSLB+BP5Yrs))
summary(lm(rich~name+yr+YrsOB+YrsSLB+BP5Yrs))
summary(lm(rich~name+yr+YrsSLB+BP5Yrs))


ob.avg <- tapply(rich[OrBis==1],yr[OrBis==1],mean)
oc.avg <- tapply(rich[OrBis==0&ChBis==0],yr[OrBis==0&ChBis==0],mean)

ob.len <- tapply(rich[OrBis==1],yr[OrBis==1],length)
oc.len <- tapply(rich[OrBis==0&ChBis==0],yr[OrBis==0&ChBis==0],length)

bis.avg <- tapply(rich[bison==1],yr[bison==1],mean)
cat.avg <- tapply(rich[bison==0],yr[bison==0],mean)

richc<-residuals(lm(rich~name))
bis.avgc <- tapply(richc[bison==1],yr[bison==1],mean)
cat.avgc <- tapply(richc[bison==0],yr[bison==0],mean)

bis.len <- tapply(rich[bison==1],yr[bison==1],length)
cat.len <- tapply(rich[bison==0],yr[bison==0],length)

##6 orginal bison and 7 orginal cattle
plot(1998:2008,ob.avg,type='o',col='blue',ylim=c(50,100))
arrows(1998:2008,ob.avg,1998:2008,ob.avg+tapply(rich[OrBis==1],yr[OrBis==1],sd)/sqrt(ob.len),angle=90,len=.05,col='blue')
arrows(1998:2008,ob.avg,1998:2008,ob.avg-tapply(rich[OrBis==1],yr[OrBis==1],sd)/sqrt(ob.len),angle=90,len=.05,col='blue')
points(1998:2008,oc.avg,type='o',col='red')
arrows(1998:2008,oc.avg,1998:2008,oc.avg+tapply(rich[OrBis==0&ChBis==0],yr[OrBis==0&ChBis==0],sd)/sqrt(oc.len),angle=90,len=.05,col='red')
arrows(1998:2008,oc.avg,1998:2008,oc.avg-tapply(rich[OrBis==0&ChBis==0],yr[OrBis==0&ChBis==0],sd)/sqrt(oc.len),angle=90,len=.05,col='red')



par(mfrow=c(2,2))
plot(1998:2008,bis.avg,type='o',xlab='',ylab='avg rich',ylim=c(50,100),col='blue')
arrows(1998:2008,bis.avg,1998:2008,bis.avg+tapply(rich[bison==1],yr[bison==1],sd)/sqrt(bis.len),angle=90,len=.05,col='blue')
arrows(1998:2008,bis.avg,1998:2008,bis.avg-tapply(rich[bison==1],yr[bison==1],sd)/sqrt(bis.len),angle=90,len=.05,col='blue')
points(1998:2008,cat.avg,type='o',col='red')
arrows(1998:2008,cat.avg,1998:2008,cat.avg+tapply(rich[bison==0],yr[bison==0],sd)/sqrt(cat.len),angle=90,len=.05,col='red')
arrows(1998:2008,cat.avg,1998:2008,cat.avg-tapply(rich[bison==0],yr[bison==0],sd)/sqrt(cat.len),angle=90,len=.05,col='red')
par(new=T)
plot(1998:2008,bis.len,type='l',lwd=1,lty=2,axes=F,xlab='',ylab='',ylim=c(0,20))
axis(side=4)

plot(1998:2008,bis.avg/bis.avg[1],type='o',xlab='',ylab='avg rich / 1998 avg rich',ylim=c(.8,1.5),col='blue')
arrows(1998:2008,bis.avg/bis.avg[1],1998:2008,bis.avg/bis.avg[1]+tapply(rich[bison==1],yr[bison==1],sd)/sqrt(bis.len)/bis.avg[1],angle=90,len=.05,col='blue')
arrows(1998:2008,bis.avg/bis.avg[1],1998:2008,bis.avg/bis.avg[1]-tapply(rich[bison==1],yr[bison==1],sd)/sqrt(bis.len)/bis.avg[1],angle=90,len=.05,col='blue')
points(1998:2008,cat.avg/cat.avg[1],type='o',col='red')
arrows(1998:2008,cat.avg/cat.avg[1],1998:2008,cat.avg/cat.avg[1]+tapply(rich[bison==0],yr[bison==0],sd)/sqrt(cat.len)/cat.avg[1],angle=90,len=.05,col='red')
arrows(1998:2008,cat.avg/cat.avg[1],1998:2008,cat.avg/cat.avg[1]-tapply(rich[bison==0],yr[bison==0],sd)/sqrt(cat.len)/cat.avg[1],angle=90,len=.05,col='red')
par(new=T)
plot(1998:2008,bis.len,type='l',lwd=1,lty=2,axes=F,xlab='',ylab='',ylim=c(0,20))
axis(side=4)

plot(1998:2008,bis.avgc,type='o',xlab='',ylab='avg rich centered on plot means',ylim=c(-25,15),col='blue')
arrows(1998:2008,bis.avgc,1998:2008,bis.avgc+tapply(richc[bison==1],yr[bison==1],sd)/sqrt(bis.len),angle=90,len=.05,col='blue')
arrows(1998:2008,bis.avgc,1998:2008,bis.avgc-tapply(richc[bison==1],yr[bison==1],sd)/sqrt(bis.len),angle=90,len=.05,col='blue')
points(1998:2008,cat.avgc,type='o',col='red')
arrows(1998:2008,cat.avgc,1998:2008,cat.avgc+tapply(richc[bison==0],yr[bison==0],sd)/sqrt(cat.len),angle=90,len=.05,col='red')
arrows(1998:2008,cat.avgc,1998:2008,cat.avgc-tapply(richc[bison==0],yr[bison==0],sd)/sqrt(cat.len),angle=90,len=.05,col='red')
par(new=T)
plot(1998:2008,bis.len,type='l',lwd=1,lty=2,axes=F,xlab='',ylab='',ylim=c(0,20))
axis(side=4)

plot(1998:2008,bis.avgc/bis.avgc[1],type='o',ylim=c(-1,1),xlab='',ylab='avg rich centered on plot means / 1998 avg centered rich',col='blue')
arrows(1998:2008,bis.avgc/bis.avgc[1],1998:2008,bis.avgc/bis.avgc[1]+tapply(richc[bison==1],yr[bison==1],sd)/sqrt(bis.len)/bis.avgc[1],angle=90,len=.05,col='blue')
arrows(1998:2008,bis.avgc/bis.avgc[1],1998:2008,bis.avgc/bis.avgc[1]-tapply(richc[bison==1],yr[bison==1],sd)/sqrt(bis.len)/bis.avgc[1],angle=90,len=.05,col='blue')
points(1998:2008,cat.avgc/cat.avgc[1],type='o',col='red')
arrows(1998:2008,cat.avgc/cat.avgc[1],1998:2008,cat.avgc/cat.avgc[1]+tapply(richc[bison==0],yr[bison==0],sd)/sqrt(cat.len)/cat.avgc[1],angle=90,len=.05,col='red')
arrows(1998:2008,cat.avgc/cat.avgc[1],1998:2008,cat.avgc/cat.avgc[1]-tapply(richc[bison==0],yr[bison==0],sd)/sqrt(cat.len)/cat.avgc[1],angle=90,len=.05,col='red')
par(new=T)
plot(1998:2008,bis.len,type='l',lwd=1,lty=2,axes=F,xlab='',ylab='',ylim=c(0,20))
axis(side=4)

rich1.mat<-as.matrix(tapply(rich1,list(plot,yr),mean))
rich1.plcen<-rich1.mat-apply(rich1.mat,1,mean) ##centered on plot means

graz.mat<-as.matrix(tapply(bison,list(plot,yr),sum))

