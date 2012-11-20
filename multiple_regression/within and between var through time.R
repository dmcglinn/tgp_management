endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
library(nlme)
library(car)
names(endat)
attach(endat)

##examining within and between corner/plot variance through time
##begin working only at the 100 m2 scale so no corners
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)

holder<-matrix(NA,nrow=11,ncol=3)
uni.yrs<-1998:2008
for(i in 1:11){
 resp<-rich[yr==uni.yrs[i]]
 group<-name[yr==uni.yrs[i]]
 holder[i,1]<-fixef(lme(resp~1,random=~1|group))
 for(j in 1:2)
  holder[i,j+1]<-as.numeric(VarCorr(lme(resp~1,random=~1|group))[j,2])
} 
colnames(holder)<-c('avg sr','between sd','error sd')
rownames(holder)<-1998:2008
holder
plot(holder[,-1],type='o')
plot(holder[,-3],type='p')
plot(holder[,-2],type='p')



rich.98<-rich[yr==1998]
name.98<-name[yr==1998]
test<-lme(rich.98~1,random=~1|name.98)
