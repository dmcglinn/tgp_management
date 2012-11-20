endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env.csv',sep=',',header=T)

names(endat)
 [1] "plot"      "yr"        "plotnum"   "yr.plot"   "X205"      "X206"     
 [7] "X208"      "X220"      "X222"      "X226"      "X238"      "X244"     
[13] "X254"      "X259"      "X303"      "X307"      "X308"      "X309"     
[19] "X317"      "X319"      "X331"      "X343"      "X346"      "X350"     
[25] "X1998"     "X1999"     "X2000"     "X2001"     "X2002"     "X2003"    
[31] "X2004"     "X2005"     "X2006"     "X2007"     "X2008"     "OrBis"    
[37] "ChBis"     "bison"     "cattle"    "YrsOB"     "BP5Yrs"    "YrsSLB"   
[43] "burn"      "spring"    "summer"    "fall"      "slope"     "northness"
[49] "eastness"  "OM"        "pH"        "SolS"      "logP"      "logCa"    
[55] "logFe"     "logMg"     "logK"      "logNa"     "logB"      "logMn"    
[61] "logCu"     "logZn"     "logAl"     "PDSIavg"   "SodTemp"   "SPI1"     
[67] "SPI12"     "SPI24"    

###################################################
##Disturbance PCA##
dist<-as.matrix(cbind(endat[,38:42],endat[,44]))
colnames(dist)<-c(colnames(endat)[38:42],colnames(endat)[44])
rownames(dist)<-endat[,3]
dist.pca<-prcomp(dist,scale=T)
plot(dist.pca)
biplot(dist.pca,cex=.5)

##consider combining fall and summer burns - together make up 20% of all burns

###################################################
##soil PCA##
soil<-endat[,50:63]
rownames(soil)<-endat[,3]
soil.pca<-prcomp(soil,scale=T)
plot(soil.pca)
biplot(soil.pca,cex=.5)

eigen<-soil.pca$sdev^2/sum(soil.pca$sdev^2)
sum(eigen[1:2])##axes 1:2 expl 52% of var
sum(eigen[1:3])##axes 1:3 expl 63% of var

##331 or South plot is a clear outlier, higher Zn and Fe
##208 or Seep is also somewhat of an outlier
##appears most variability is with respect to Ca and Fe
##Ca cor Cu
##Fe cor Al
##Zn cor P
##B cor Mg
##So consider dropping Cu, Al, P, Mg
#####################################################
#####################################################
##climate vars PCA##
clim<-as.matrix(endat[,64:68])
rownames(clim)<-endat[,3]
clim.pca<-prcomp(clim,scale=T)
plot(clim.pca)
biplot(clim.pca,cex=.5)

eigen<-clim.pca$sdev^2/sum(clim.pca$sdev^2)
sum(eigen[1:2])##axes 1:2 expl 83% of var
sum(eigen[1:3])##axes 1:3 expl 96% of var

#####################################################
#####################################################
##site vars PCA##
site<-endat[,47:63]
rownames(site)<-endat[,3]
site.pca<-prcomp(site,scale=T)
plot(site.pca)
biplot(site.pca,cex=.5)

eigen<-site.pca$sdev^2/sum(site.pca$sdev^2)
sum(eigen[1:2])##axes 1:2 expl 47% of var
sum(eigen[1:3])##axes 1:3 expl 57% of var

##including slope and aspect with the soils did not change much
##South plot still a clear outlier, as well as 309 or west buzz
##Ca cor OM and pH
##Mg cor B and eastness
##Zn cor P
#####################################################
#####################################################
##change in soil vars through time##
##first see if a yr dummie var shows diffs between years
p.vals<-rep(NA,length(49:62))
icnt<-0
for(i in 1:length(p.vals)){
 icnt<-icnt+1
 fvals<-summary(lm(endat[,48+i]~as.matrix(endat[,24:34])))$fstat
 p.vals[icnt]<-1-pf(fvals[1],fvals[2],fvals[3]) 
}
as.data.frame(rbind(colnames(endat[,49:62]),round(p.vals,3),p.vals<0.05))
1    OM    pH SolS logP logCa logFe logMg logK logNa  logB logMn logCu logZn logAl
2 0.143 0.959    0    0 0.721 0.018 0.242    0 0.145 0.088  0.46 0.002     0     0
3 FALSE FALSE TRUE TRUE FALSE  TRUE FALSE TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE
##check for linear temp trends
yr<-endat[,2]

p.vals<-rep(NA,length(49:62))
icnt<-0
for(i in 1:length(p.vals)){
 icnt<-icnt+1
 p.vals[icnt]<-summary(lm(endat[,48+i]~yr))$coef[2,4] 
}
as.data.frame(rbind(colnames(endat[,49:62]),round(p.vals,3),p.vals<0.05))
1    OM    pH SolS logP logCa logFe logMg logK logNa  logB logMn logCu logZn logAl
2 0.044 0.127    0    0 0.268 0.022 0.343    0 0.509  0.97 0.102 0.026     0 0.004
3  TRUE FALSE TRUE TRUE FALSE  TRUE FALSE TRUE FALSE FALSE FALSE  TRUE  TRUE  TRUE
sum(p.vals<.05) ##8
##graphics for temporal trend
par(mfrow=c(3,5))
par(mar=c(2, 4, 2, 2))
for(i in 1:length(p.vals)){
 if(p.vals[i]<0.05){
  plot(yr,endat[,48+i],ylab=colnames(endat)[48+i],xlab='',main='**')
  abline(lm(endat[,48+i]~yr),col='blue')
  lines(lowess(yr,endat[,48+i]),col='red')
 }
 else{
 plot(yr,endat[,48+i],ylab=colnames(endat)[48+i],xlab='')
 abline(lm(endat[,48+i]~yr),col='blue')
 lines(lowess(yr,endat[,48+i]),col='red')
}}
##Reperform analysis and center by sites##
plots<-unique(endat[,1])
soil<-matrix(NA,nrow=dim(endat)[1],ncol=length(49:62))##centered matrix
icnt<-0
for(i in 1:11){ 
 for(ii in 1:20){
  icnt<-icnt+1
  for(j in 1:length(49:62)){
  soil[icnt,j]<-endat[endat[,1]==plots[ii],48+j][i]-mean(endat[endat[,1]==plots[ii],48+j])
}}}

p.vals<-rep(NA,length(49:62))
icnt<-0
for(i in 1:length(p.vals)){
 icnt<-icnt+1
 p.vals[icnt]<-summary(lm(soil[,i]~yr))$coef[2,4] 
}
as.data.frame(rbind(colnames(endat[,49:62]),round(p.vals,3),p.vals<0.05))
##graphics for temporal trend
par(mfrow=c(3,5))
par(mar=c(2, 4, 2, 2))
for(i in 1:length(p.vals)){
 if(p.vals[i]<0.05){
  plot(yr,soil[,i],ylab=colnames(endat)[48+i],xlab='',main='**')
  abline(lm(soil[,i]~yr),col='blue')
  lines(lowess(yr,soil[,i]),col='red')
 }
 else{
 plot(yr,soil[,i],ylab=colnames(endat)[48+i],xlab='')
 abline(lm(soil[,i]~yr),col='blue')
 lines(lowess(yr,soil[,i]),col='red')
}}
#####################################################
#####################################################
##time since fire and soil##
tsb<-endat[,41]##time since burn

p.vals<-rep(NA,length(49:62))
icnt<-0
for(i in 1:length(p.vals)){
 icnt<-icnt+1
 p.vals[icnt]<-summary(lm(soil[,i]~tsb))$coef[2,4] 
}
as.data.frame(rbind(colnames(endat[,49:62]),round(p.vals,3),p.vals<0.05))
1    OM    pH  SolS  logP logCa logFe logMg  logK logNa  logB logMn logCu logZn logAl
2 0.801 0.274 0.049 0.036 0.974 0.366 0.841 0.263 0.493 0.427 0.679 0.544 0.655 0.746
3 FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

##graphics for temporal trend
par(mfrow=c(3,5))
par(mar=c(2, 4, 2, 2))
for(i in 1:length(p.vals)){
 if(p.vals[i]<0.05){
  plot(tsb,soil[,i],ylab=colnames(endat)[48+i],xlab='',main='**')
  abline(lm(soil[,i]~tsb),col='blue')
  lines(lowess(tsb,soil[,i]),col='red')
 }
 else{
 plot(tsb,soil[,i],ylab=colnames(endat)[48+i],xlab='')
 abline(lm(soil[,i]~tsb),col='blue')
 lines(lowess(tsb,soil[,i]),col='red')
}}
##no strong trends#####
#####################################################
#####################################################