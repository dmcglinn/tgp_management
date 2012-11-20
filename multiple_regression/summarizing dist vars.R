spdat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
attach(endat)
##begin working only at the 100 m2 scale so no corners
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)

YrsSLB[name=="Muddy Run"]
YrsSLB[name=="Stadia"]


hist(YrsSLB,breaks=0:12,main='',xlab='',ylab='',right=F)
par(new=T)
hist(YrsSLB[bison==1],breaks=0:12,ylim=c(0,80),col='red',axes=F,,main='',xlab='',ylab='',right=F)
par(new=T)
hist(YrsSLB[bison==0],breaks=0:12,ylim=c(0,80),density=12,angle=45,border='blue',axes=F,main='',xlab='',ylab='',right=F)

hist(BP5Yrs,breaks=0:6,main='',xlab='',ylab='',right=F)
par(new=T)
hist(BP5Yrs[bison==1],breaks=0:6,ylim=c(0,80),col='red',axes=F,,main='',xlab='',ylab='',right=F)
par(new=T)
hist(BP5Yrs[bison==0],breaks=0:6,ylim=c(0,80),density=12,angle=45,border='blue',axes=F,main='',xlab='',ylab='',right=F)

panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="grey", ...)
}
panel.cor <- function(x, y, digits=2, prefix="", cex.cor)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- (cor(x, y))
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}

pairs(~YrsOB+YrsSLB+BP5Yrs,lower.panel=panel.smooth,diag.panel=panel.hist,upper.panel=panel.cor,labels=c("Years of Bison","Years since fire","# of burns in past 5 years"))

##table of # of burns each year of the study
bcount<-tapply(burn,list(bison,yr),sum)
  1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008
0    9    6    6    7    5    5    3    6    3    0    5
1    1    0    3    3    2    1    3    5    2    2    3

yrsavg<-tapply(YrsSLB,list(bison,yr),mean)
  1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008
0  1.2  1.3  1.5  1.4  1.7  2.0  2.6  1.3  1.8  2.9  2.1
1  2.5  3.5  3.2  2.0  1.8  2.5  2.7  2.3  2.3  2.8  2.3

se.func<-function(x){ 
 hold<- glm(x~1,family='poisson') 
 exp(summary(hold)$coef[2])
} 
##for a right skewed dist this prob is not very meaningful
round(yrsse<-tapply(YrsSLB,list(bison,yr),se.func),2)
  1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008
0 1.28 1.27 1.27 1.28 1.25 1.27 1.23 1.39 1.32 1.25  1.3
1 1.29 1.24 1.22 1.29 1.30 1.21 1.20 1.20 1.20 1.18  1.2

test<-glm(YrsSLB[bison==1]~1,family='poisson')
exp(coef(test))
exp(summary(test)$coef[2])
test<-glm(YrsSLB[bison==0]~1,family='poisson')
exp(coef(test))
exp(summary(test)$coef[2])

table(bison,yr)
     yr
bison 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008
    0   14   14   12   12   12    9    9    7    7    7    7
    1    6    6    8    8    8   11   11   13   13   13   13


test<-glm(round(YrsSLB)~1,family='poisson')

hist(YrsSLB,breaks=0:15,right=F,ylim=c(0,100))
par(new=T)

X<-matrix(NA,ncol=1000,nrow=220)
for(i in 1:dim(X)[2]){
 X[,i]<-sort(rpois(220,lambda=exp(coef(test))))
}
x<-apply(X,1,mean)
hist(x,breaks=0:15,right=F,border='blue',ylim=c(0,100))

test<-glm(BP5Yrs~1,family='poisson')
hist(BP5Yrs,breaks=0:8,right=F,ylim=c(0,100))
par(new=T)

X<-matrix(NA,ncol=500,nrow=220)
for(i in 1:dim(X)[2]){
 X[,i]<-sort(rpois(220,lambda=exp(coef(test))))
}
x<-apply(X,1,mean)
hist(x,breaks=0:8,right=F,border='blue',ylim=c(0,100))

