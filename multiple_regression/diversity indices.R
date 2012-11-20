endat<-read.table('C:/Users/dmcglinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
spdat<-read.table('C:/Users/dmcglinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
#endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
#spdat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
sp<-spdat[,-1]
library(nlme)
library(car)
names(endat)
attach(endat)

endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)


library(vegan)
H<-diversity(spdat[,-1])
H.avg<-tapply(diversity(spdat[,-1]),yr,mean)

cor(H,rich)
.65
mod<-lm(H~name+dist.m+rain.tot*temp.avg)
summary(mod)$r.s
summary(lm(rich~name+yr.f+dist.m))$r.s

##see Wilsey et al. 2005
rarity<-function(x){ ##as a proportion of total diversity
 x <- x[x>0]
 rich<-length(x)
 rel<-x/sum(x)
 sum(ifelse(rel<(1/rich),1,0))/rich
}
test<-matrix(c(10,0,1,3,0,10,3,1),ncol=2)
apply(test,2,rarity)
##Berg Parker, maximum proportional biomass
BP<-apply(spdat[,-1],1,function(x)max(x)/sum(x))
D<-diversity(spdat[,-1],"simp")
E<-D/rich
R<-apply(spdat[,-1],1,rarity)
plot(H,D)
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
pairs(~rich+R+H+D+E+BP,upper.panel=panel.cor,lower.panel=panel.smooth,diag.panel=NULL,labels=c('Richness','Rarity','Shannon','Simpson','Evenness','Berg Dom'))

par(mfrow=c(1,5))
plot(rich,H,xlab='Richness',ylab="Shannon's diversity (H')")
lines(lowess(rich,H),col='red',lwd=2)
plot(rich,D,xlab='Richness',ylab="Simpson's diversity (D)")
lines(lowess(rich,D),col='red',lwd=2)
plot(rich,R,xlab='Richness',ylab='proportion of rare species (Rarity)')
lines(lowess(rich,R),col='red',lwd=2)
plot(rich,E,xlab='Richness',ylab='Evenness (E)')
lines(lowess(rich,E),col='red',lwd=2)
plot(rich,BP,xlab='Richness',ylab='Berger-Parker dominance (BP)')
lines(lowess(rich,BP),col='red',lwd=2)



par(mfrow=c(1,5))
[1] 5.1 4.1 4.1 2.1
c(bottom, left, top, right) 
par(mar=c(5,2,3,0))
plot(R,rich,xlab='Rarity')
lines(lowess(R,rich),col='red')
par(mar=c(5,0,3,0))
plot(H,rich,xlab='Shannon',axes=F,frame.plot=T)
axis(side=1)
lines(lowess(H,rich),col='red')
plot(D,rich,xlab='Simpson',axes=F,frame.plot=T)
axis(side=1)
lines(lowess(D,rich),col='red')
plot(E,rich,xlab='Evenness',axes=F,frame.plot=T)
axis(side=1)
lines(lowess(E,rich),col='red')
par(mar=c(5,0,3,1))
plot(BP,rich,xlab='Berg Dom',axes=F,frame.plot=T)
axis(side=1)
lines(lowess(BP,rich),col='red')
##.19,.66,.39,-.96,-.28


##Wilsely analysis
pc<-princomp(cbind(rich,R,D,E,BP))
rownames(pc$scores)<-rep(as.character(""),220)
plot(pc)
biplot(pc,scale=T,choices=c(1,2))
biplot(pc,scale=T,choices=c(1,2),ylim=c(-0.001,0.001),xlim=c(-0.001,0.001))
round(pc$loadings,3)
pc2<-prcomp(cbind(rich,R,D,E,BP))
plot(pc2)
biplot(pc2,scale=T)
round(pc2$rotation,3)
        PC1    PC2    PC3    PC4    PC5
rich -1.000 -0.002  0.001  0.000  0.000
R    -0.001 -0.061 -0.998 -0.022  0.000
D    -0.001  0.377 -0.043  0.925 -0.013
E     0.000  0.004  0.000  0.012  1.000
BP    0.002 -0.924  0.048  0.378 -0.001
summary(pc2)
Importance of components:
                          PC1    PC2     PC3    PC4      PC5
Standard deviation     17.621 2.1809 0.08095 0.0123 0.000407
Proportion of Variance  0.985 0.0151 0.00002 0.0000 0.000000
Cumulative Proportion   0.985 1.0000 1.00000 1.0000 1.000000



mod<-lm(D~name+yr.f+dist.m)
summary(mod)$r.s

sr.avg<-tapply(rich,yr,mean)
D.avg<-tapply(D,yr,mean)
plot(1998:2008,H.avg,type='o',ylim=c(3,4.5))
points(1998:2008,(sr.avg/max(sr.avg))+3,pch=19,type='o')
points(1998:2008,D.avg+3,pch=19,col='red',type='o')

site.m<-endatlv1[,8:27]
yr.m<-endatlv1[,28:38]
dist.m<-cbind(YrsOB,bison,YrsSLB,BP5Yrs)
plot(varpart(E,yr.m,site.m,dist.m))

