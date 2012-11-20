mod.covars<-lm(rich1~dist+site+grazer)
resids<-residuals(mod.covars)

clim<-as.matrix(cbind(endat[,69:77],endat[,80:82],endat[,84:92]))
clim<-as.matrix(endat[,69:80])
clim<-as.matrix(endat[,sample(69:80)])
colnames(clim)

full.mod<-lm(resids~clim)
min.mod<-lm(resids~1)

scope=list(lower=~1,upper=~clim)

add1(min.mod,test="Chisq",~temp1+temp2+temp3+temp4+temp5+temp6+temp7+temp8+temp9+temp10+temp11+temp12)
step.temp.lm<-step(min.mod,direction='forward',~temp1+temp2+temp3+temp4+temp5+temp6+temp7+temp8+temp9+temp10+temp11+temp12)

step.rain.lm<-step(min.mod,direction='forward',~rain1+rain2+rain3+rain4+rain5+rain6+rain7+rain8+rain9+rain10+rain11+rain12)

summary(full.mod)
step(full.mod)
full.mod<-lm(resids~rain6+rain7+rain8+rain9+rain10+rain11+rain12+rain1+rain2+rain3+rain4+rain5)
full.mod<-lm(resids~temp6+temp7+temp8+temp9+temp10+temp11+temp12+temp1+temp2+temp3+temp4+temp5)
full.mod<-lm(resids~sum.rain*win.rain* spr.rain*sum.temp* win.temp*spr.temp)

clim2<-as.matrix(cbind(sum.rain,win.rain, spr.rain,sum.temp, win.temp,spr.temp,rain.avg , temp.avg))
full.mod<-lm(resids~ clim2)


step.lm<-step(full.mod)
summary(step.lm)

clim2.std<-clim2
for(j in 1:dim(clim2)[2]){
for(i in 1:dim(clim2)[1]){
 clim2.std[i,j]<-(clim2[i,j]-apply(clim2,2,mean)[j])/apply(clim2,2,sd)[j]
}}
apply(clim2.std,2,mean)
apply(clim2.std,2,sd)

resids.std<-(resids-mean(resids))/sd(resids)


full.mod<-lm(resids.std~-1+ clim2.std)
full.mod<-lm(resids.std~-1+clim2)
summary(full.mod)

step.lm<-step(full.mod)
summary(step.lm)


##climate PCA##
rownames(clim2)<-endat[,3]
clim.pca<-prcomp(clim2,scale=T)
plot(clim.pca)
biplot(clim.pca,cex=.5)

eigen<-clim.pca$sdev^2/sum(clim.pca$sdev^2)
sum(eigen[1:2])##axes 1:2 expl 83% of var
sum(eigen[1:3])##axes 1:3 expl 96% of var




scope = list(upper = ~Eth*Sex*Age*Lrn, lower = ~1),


resids<-runif(dim(clim2)[1])
clim2<-matrix(rnorm(dim(clim2)[1]*dim(clim2)[2],mean=10)*runif(dim(clim2)[1]*dim(clim2)[2]),nrow=dim(clim2)[1],ncol=dim(clim2)[2])




