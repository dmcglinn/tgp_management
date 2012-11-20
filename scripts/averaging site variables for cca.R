
site.m<-as.matrix(endatlv1[,50:66])
site.a<-site.m
for(i in 1:ncol(site.m))
 site.a[,i]<-rep(tapply(site.m[,i],name,mean),each=11)

site.pc<-prcomp(site.a,scale=T)
plot(site.pc)
biplot(site.pc,cex=.5)
summary(site.pc)
site.pc1<-site.pc$x[,1]

soil.pc<-prcomp(site.a[,-(1:3)],scale=T)
biplot(soil.pc)
topo.pc<-prcomp(site.a[,1:2])
biplot(topo.pc)

soil.pc1<-soil.pc$x[,1]
topo.pc1<-topo.pc$x[,1:2]

cca.site<-cca(spdat.n~soil.pc1+topo.pc1+Condition(YrsOB)+
             +Condition(YrsSLB)+Condition(BP5Yrs)+Condition(yr.f))
p.int(cca.site)
plot(cca.site)

test<-cca(spdat.n,site.m,cbind(dist.m,yr.m))
test2<-cca(spdat.n,site.a,cbind(dist.m,yr.m))
p.int(test)
p.int(test2)


ca.a<-rep(tapply(logCa,name,mean),each=11)
fe.a<-rep(tapply(logFe,name,mean),each=11)
ph.a<-rep(tapply(pH,name,mean),each=11)
zn.a<-rep(tapply(logZn,name,mean),each=11)


summary(logCa)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  2.886   3.296   3.423   3.416   3.561   3.699 
plot(NA,ylim=c(0,5),xlim=c(2.8,3.8),ylab='cover',xlab='logCa')
for(i in 1:10){
 lines(lowess(logCa,spdat.n[,i],f=.66))
}
summary(logFe)

plot(NA,ylim=c(0,4),xlim=c(1.8,2.8),ylab='cover',xlab='logFe')
for(i in 1:5){
 lines(lowess(logFe,spdat.n[,i],f=.66))
}

plot(NA,ylim=c(0,5),xlim=c(0,5),ylab='cover',xlab='YrsSLB')
for(i in 1:5){
 lines(lowess(BP5Yrs,spdat.n[,i],f=.1))
}

plot(NA,ylim=c(0,5),xlim=c(0,15),ylab='cover',xlab='YrsSLB')
for(i in 1:5){
 lines(lowess(YrsOB,spdat.n[,i],f=.66))
}

 [1] "andrgera" "schiscop" "sporcomp" "ambrpsil" "sorgnuta" "solicana"
 [7] "panivirg" "ambrbide" "bromjapo" "elymvirg"

