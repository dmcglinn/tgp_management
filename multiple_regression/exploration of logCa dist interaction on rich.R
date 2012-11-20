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

summary(lm(residuals(lm(rich~factor(yr)))~logCa*BP5Yrs))

logCa.s<-scale(logCa)
bp5yrs<-scale(BP5Yrs)
summary(lm(rich~logCa+BP5Yrs))
summary(lm(rich~logCa*BP5Yrs))

Ca<-cut(logCa,quantile(logCa),include.lowest=T)
freqB<-cut(BP5Yrs,quantile(BP5Yrs),include.lowest=T)
class(freqB)
interaction.plot(Ca,freqB,rich)
interaction.plot(freqB,Ca,rich)

summary(lm(rich~logCa.s+bp5yrs))
summary(lm(rich~logCa.s*bp5yrs))
Ca<-cut(logCa.s,quantile(logCa.s),include.lowest=T)
freqB<-cut(bp5yrs,quantile(bp5yrs),include.lowest=T)
class(freqB)
interaction.plot(Ca,freqB,rich)
interaction.plot(freqB,Ca,rich)
