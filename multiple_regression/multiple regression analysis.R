endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env.csv',sep=',',header=T)
spdat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp sp.csv',sep=',',header=T)
names(endat)
  [1] "plot"        "yr"          "plotnum"     "yr.plot"     "X205"       
  [6] "X206"        "X208"        "X220"        "X222"        "X226"       
 [11] "X238"        "X244"        "X254"        "X259"        "X303"       
 [16] "X307"        "X308"        "X309"        "X317"        "X319"       
 [21] "X331"        "X343"        "X346"        "X350"        "X1998"      
 [26] "X1999"       "X2000"       "X2001"       "X2002"       "X2003"      
 [31] "X2004"       "X2005"       "X2006"       "X2007"       "X2008"      
 [36] "OrBis"       "ChBis"       "bison"       "cattle"      "YrsOB"      
 [41] "BP5Yrs"      "YrsSLB"      "burn"        "spring"      "summer"     
 [46] "fall"        "slope"       "northness"   "eastness"    "OM"         
 [51] "pH"          "SolS"        "logP"        "logCa"       "logFe"      
 [56] "logMg"       "logK"        "logNa"       "logB"        "logMn"      
 [61] "logCu"       "logZn"       "logAl"       "PDSIavg"     "SodTemp"    
 [66] "SPI1"        "SPI12"       "SPI24"       "rain1"       "rain2"      
 [71] "rain3"       "rain4"       "rain5"       "rain6"       "rain7"      
 [76] "rain8"       "rain9"       "rain10"      "rain11"      "rain12"     
 [81] "temp1"       "temp2"       "temp3"       "temp4"       "temp5"      
 [86] "temp6"       "temp7"       "temp8"       "temp9"       "temp10"     
 [91] "temp11"      "temp12"      "sum.rain"    "win.rain"    "spr.rain"   
 [96] "sum.temp"    "win.temp"    "spr.temp"    "rain.avg"    "temp.avg"   
[101] "rich1"       "rich2"       "rich3"       "rich4"       "rich5"      
[106] "A.cov"       "P.cov"       "S.cov"       "T.cov"       "A.p"        
[111] "P.p"         "S.p"         "T.p"         "ForbLeg.cov" "Forb.cov"   
[116] "Gsum.cov"    "G3.cov"      "G4.cov"      "Leg.cov"     "Wood.cov"   
[121] "F.p"         "Gsum.p"      "G3.p"        "G4.p"        "L.p"        
[126] "W.p"  

##first examine multicollinarity in expl vars
season<-as.factor(ifelse(endat$spring==1,'spring','other'))
grazer<-as.factor(ifelse(endat$bison==1,'bison','cattle'))
mang.cat<-rep(NA,dim(endat)[1])
for(i in 1:length(mang.cat)){
 if(endat$OrBis[i]==1)
  mang.cat[i]<-'OBis'
 if(endat$ChBis[i]==1)
  mang.cat[i]<-'CBis'
 if(endat$OrBis[i]+endat$ChBis[i]==0)
  mang.cat[i]<-'OCat'
}
mang.cat<-as.factor(mang.cat)

topo<-as.matrix(endat[,47:49])
clim<-as.matrix(endat[,69:92])
dist<-as.matrix(endat[,40:42])
site<-as.matrix(endat[,5:24])

attach(endat)

##climate PCA##
rownames(clim)<-endat[,3]
symnum(cor(clim))
clim.pca<-prcomp(clim,scale=T)
plot(clim.pca)
biplot(clim.pca,cex=.5)

eigen<-clim.pca$sdev^2/sum(clim.pca$sdev^2)
sum(eigen[1:2])##axes 1:2 expl 83% of var
sum(eigen[1:3])##axes 1:3 expl 96% of var

##drop temp wt, drop SPI1&4
##############################
##############################

mod.site<-lm(rich~site)
mod.full<-lm(residuals(mod.site)~season+grazer+dist+OM+logCa+logFe+topo+rain.avg*temp.avg+PDSIavg+I(as.matrix(PDSIavg)^2))
Coefficients:
                          Estimate Std. Error t value Pr(>|t|)    
(Intercept)              732.33842  185.30516   3.952 0.000107 ***
seasonspring              -0.34969    1.47385  -0.237 0.812694    
grazercattle               5.60217    1.63231   3.432 0.000726 ***
distYrsOB                  0.66546    0.19691   3.380 0.000870 ***
distBP5Yrs                -1.56630    0.71816  -2.181 0.030333 *  
distYrsSLB                -1.03538    0.35697  -2.900 0.004136 ** 
OM                         0.52995    0.82447   0.643 0.521090    
logCa                     -1.57204    5.81686  -0.270 0.787238    
logFe                      2.18715    4.08022   0.536 0.592521    
toposlope                  0.02136    0.35602   0.060 0.952228    
toponorthness              0.72030    0.71351   1.010 0.313930    
topoeastness              -0.42651    1.09062  -0.391 0.696152    
rain.avg                -233.26399   51.46322  -4.533 9.94e-06 ***
temp.avg                 -12.48839    3.16259  -3.949 0.000108 ***
PDSIavg                    0.65566    1.06362   0.616 0.538291    
I(as.matrix(PDSIavg)^2)   -1.53772    0.42117  -3.651 0.000332 ***
rain.avg:temp.avg          3.98902    0.88114   4.527 1.02e-05 ***

mod.PDSI<-lm(rich~PDSIavg)
mod.PDSI2<-update(mod.PDSI,.~.+I(as.matrix(PDSIavg)^2))
anova(mod.PDSI,mod.PDSI2)

plot(rich~PDSIavg)
points(by(rich,yr,mean)~by(PDSIavg,yr,mean),pch=17,col='blue')
lines(lowess(PDSIavg,rich),col='green',lwd=2)
pdsi<-seq(min(PDSIavg),max(PDSIavg),.1)
lines(pdsi,coef(mod.PDSI2)[1]+coef(mod.PDSI2)[2]*pdsi+coef(mod.PDSI2)[3]*pdsi^2,col='blue',lwd=2)
lines(pdsi,coef(mod.PDSI)[1]+coef(mod.PDSI)[2]*pdsi,col='red',lwd=2)


mod.covars<-lm(rich1~dist+site+grazer)
mod.PDSI<-lm(residuals(mod.covars)~PDSIavg)
mod.PDSI2<-update(mod.PDSI,.~.+I(as.matrix(PDSIavg)^2))
anova(mod.PDSI,mod.PDSI2)

plot(residuals(mod.covars)~PDSIavg)
points(by(residuals(mod.covars),yr,mean)~by(PDSIavg,yr,mean),pch=17,col='blue')
lines(lowess(PDSIavg,residuals(mod.covars)),col='green',lwd=2)
pdsi<-seq(min(PDSIavg),max(PDSIavg),.1)
lines(pdsi,coef(mod.PDSI2)[1]+coef(mod.PDSI2)[2]*pdsi+coef(mod.PDSI2)[3]*pdsi^2,col='blue',lwd=2)
lines(pdsi,coef(mod.PDSI)[1]+coef(mod.PDSI)[2]*pdsi,col='red',lwd=2)
summary(mod.PDSI2)

mod.SodTemp<-lm(rich~SodTemp)
mod.SodTemp2<-update(mod.SodTemp,.~.+I(as.matrix(SodTemp)^2))
anova(mod.SodTemp,mod.SodTemp2)
plot(rich~SodTemp)
abline(mod.SodTemp)

mod.SodTemp<-lm(residuals(mod.covars)~SodTemp)
mod.SodTemp2<-update(mod.SodTemp,.~.+I(as.matrix(SodTemp)^2))
anova(mod.SodTemp,mod.SodTemp2)
plot(residuals(mod.covars)~SodTemp)
abline(mod.SodTemp)

mod.temp.avg<-lm(residuals(mod.covars)~temp.avg)
mod.temp.avg2<-update(mod.temp.avg,.~.+I(as.matrix(temp.avg)^2))
anova(mod.temp.avg,mod.temp.avg2)
plot(residuals(mod.covars)~temp.avg)
abline(mod.temp.avg)
summary(mod.temp.avg)

mod.temp.rain<-lm(residuals(mod.covars)~temp.avg+rain.avg)
mod.temp.rain.int<-update(mod.temp.rain,.~.+temp.avg:rain.avg)
anova(mod.temp.rain,mod.temp.rain.int)

library(scatterplot3d)
test<-scatterplot3d(temp.avg,rain.avg,residuals(mod.covars))
test$plane3d(mod.temp.rain)

library(lattice)

mod.covars<-lm(rich3~dist+site+grazer)
mod.temp.rain.int<-lm(residuals(mod.covars)~temp.avg*rain.avg)
coefs<-coef(mod.temp.rain.int)
mod.funct<-function(x,y){ coefs[1]+coefs[2]*x+coefs[3]*y+coefs[4]*x*y }

temp<-seq(min(temp.avg),max(temp.avg),len=20)
rain<-seq(min(rain.avg),max(rain.avg),len=20)
g<-expand.grid(temp = temp, rain = rain)

g$z <- mod.funct(g$temp,g$rain)
wireframe(z~temp*rain,data= g, drape = TRUE,  colorkey = TRUE)
#panel.cloud(temp.avg,rain.avg,1,residuals(mod.covars))
#cloud(by(residuals(mod.covars),yr,mean)~by(temp.avg,yr,mean)*by(rain.avg,yr,mean))





z.mod<-outer(temp,rain,mod.funct)
p3d<-persp(temp, rain, z.mod, theta = 115, phi = 5,ticktype = "detailed")


library(rgl)
colorlut <- terrain.colors(round(max(z.mod))) # height color lookup table
#library(dichromat)
colorlut <- colorRampPalette(c("blue","skyblue","lightblue"),space="rgb")(length(rain*2)) #set up colors template
col <- colorlut[ round(-min(z.mod)+z.mod) ] # assign colors to heights for each point

plot3d(temp.avg,rain.avg,residuals(mod.covars),size=3)
rgl.material(alpha=.5)
surface3d(temp,rain,z.mod,col=col)
points3d(by(temp.avg,yr,mean),by(rain.avg,yr,mean),by(residuals(mod.covars),yr,mean),size=5,col='red')
play3d(spin3d(),duration=10)
movie3d( spin3d(), duration=10,dir="C:/R/R-2.7.2/bin/movie",clean=FALSE)


p3d<-persp3d(temp, rain, z.mod, theta = 115, phi = 5,ticktype = "detailed")


contour(temp,rain,z)
image(temp,rain,z,drawlabels = FALSE)

mod.clim<-lm(rich~clim+I(as.matrix(PDSIavg^2)))
summary(mod.clim)


