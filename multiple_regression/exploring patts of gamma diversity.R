spdat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
attach(endat)
##begin working only at the 100 m2 scale so no corners
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)

##totrich defined within a year
totrich<-rep(0,11)
icount<-0
for(i in 1998:2008){
 icount<-icount+1
 totrich[icount]<-sum(ifelse(apply(spdat[yr==i,-1],2,sum)>0,1,0))
}

plot(1998:2008,totrich,type='o')
summary(lm(totrich~c(-5:5)))
##not sig
par(mfrow=c(1,2))
plot(1998:2008,totrich/tapply(rich,yr,mean))
lines(lowess(1998:2008,totrich/tapply(rich,yr,mean)),col='red')
plot(1998:2008,totrich-tapply(rich,yr,mean))
lines(lowess(1998:2008,totrich-tapply(rich,yr,mean)),col='red')
##totrich defined within a plot
uni.name<-unique(name)
totrich<-rep(0,20)
for(i in 1:20){
 totrich[i]<-sum(ifelse(apply(spdat[name==uni.name[i],-1],2,sum)>0,1,0))
}

plot(yr,rep(totrich,each=11)-rich)
lines(lowess(yr,rep(totrich,each=11)-rich),lwd=2,col='red')

##for calc beta div
gamma<-rep(totrich,20)

##three ways to calculate beta diversity from gramma and alpha
par(mfrow=c(2,2))
plot(yr,log(gamma)-log(rich))
lines(lowess(yr,log(gamma)-log(rich),f=.1),lwd=2,col='red')
plot(yr,gamma-rich)
lines(lowess(yr,gamma-rich,f=.1),lwd=2,col='red')

##for only plots orginally in bison
#par(mfrow=c(2,2))
plot(yr[OrBis==1],log(gamma[OrBis==1])-log(rich[OrBis==1]),col='blue')
lines(lowess(yr[OrBis==1],log(gamma[OrBis==1])-log(rich[OrBis==1]),f=.1),lwd=2,col='blue')
points(yr[OrBis==0&ChBis==0],log(gamma[OrBis==0&ChBis==0])-log(rich[OrBis==0&ChBis==0]),col='red')
lines(lowess(yr[OrBis==0&ChBis==0],log(gamma[OrBis==0&ChBis==0])-log(rich[OrBis==0&ChBis==0]),f=.1),lwd=2,col='red')
plot(yr[OrBis==1],gamma[OrBis==1]-rich[OrBis==1],col='blue')
lines(lowess(yr[OrBis==1],gamma[OrBis==1]-rich[OrBis==1],f=.1),lwd=2,col='blue')
points(yr[OrBis==0&ChBis==0],gamma[OrBis==0&ChBis==0]-rich[OrBis==0&ChBis==0],col='red')
lines(lowess(yr[OrBis==0&ChBis==0],gamma[OrBis==0&ChBis==0]-rich[OrBis==0&ChBis==0],f=.1),lwd=2,col='red')

plot(yr[OrBis==1],(gamma[OrBis==1])/(rich[OrBis==1]),col='blue')
lines(lowess(yr[OrBis==1],(gamma[OrBis==1])/(rich[OrBis==1]),f=.1),lwd=2,col='blue')
points(yr[OrBis==0&ChBis==0],(gamma[OrBis==0&ChBis==0])/(rich[OrBis==0&ChBis==0]),col='red')
lines(lowess(yr[OrBis==0&ChBis==0],(gamma[OrBis==0&ChBis==0])/(rich[OrBis==0&ChBis==0]),f=.1),lwd=2,col='red')
plot(yr[OrBis==1],gamma[OrBis==1]-rich[OrBis==1],col='blue')
lines(lowess(yr[OrBis==1],gamma[OrBis==1]-rich[OrBis==1],f=.1),lwd=2,col='blue')
points(yr[OrBis==0&ChBis==0],gamma[OrBis==0&ChBis==0]-rich[OrBis==0&ChBis==0],col='red')
lines(lowess(yr[OrBis==0&ChBis==0],gamma[OrBis==0&ChBis==0]-rich[OrBis==0&ChBis==0],f=.1),lwd=2,col='red')



