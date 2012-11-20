
mang<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/LTREB proposal/vegmon management variables.csv',sep=',',header=T)

names(mang)
 [1] "plotnum" "OB"      "CB"      "bbal"    "plot"    "plotna"  "year"   
 [8] "sr"      "sr.old"  "juld"    "grazer"  "bis"     "cat"     "yob"    
[15] "b5yr"    "yslb"    "season"  "spr"     "gro"     "win"     "pH"     
[22] "log.ca." "log.fe."

attach(mang)

library(nlme)

pl<-factor(plotna)
contrasts(pl)<-contr.helmert ##default is contr.treatment

contrasts(grazer)<-c(1,0) #one for bison, zero for cow

cyear<-year-mean(year)
mod1.lm<-lm(sr~cyear+b5yr+yslb+log.ca.+grazer)
mod2.lm<-update(mod1.lm,.~.+pl)

summary(mod1.lm)
Coefficients:
Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept) 221.3669    17.4823  12.662  < 2e-16 ***
cyear         1.6363     0.4070   4.020 9.66e-05 ***
b5yr         -2.9827     0.9262  -3.220  0.00161 ** 
yslb         -2.4844     0.5996  -4.144 6.02e-05 ***
log.ca.     -40.1982     5.0432  -7.971 6.10e-13 ***
grazer1       5.5747     1.9914   2.799  0.00588 **  ##Bison effect
Residual standard error: 9.322 on 134 degrees of freedom
Multiple R-squared: 0.6835,     Adjusted R-squared: 0.6175 
F-statistic: 10.35 on 24 and 115 DF,  p-value: < 2.2e-16 

anova(mod1.lm)

##the plots are the random variables we want to control for
gmang<-groupedData(sr~cyear+b5yr+yslb+log.ca.|pl)


mod1.lis<-lmList(gmang)
plot(intervals(mod1.lis))
plot(mod1.lis)
pairs(mod1.lis,id=0.01)
summary(mod1.lis)$RSE
[1] 6.614478
summary(mod1.lis)$df.residual
[1] 40

##mixed effects models
mod1.lme<-lme(sr~cyear+b5yr+yslb+log.ca.+grazer,random=~1|pl,data=gmang)
mod2.lme<-update(mod1.lme,random=~cyear|pl)
mod3.lme<-update(mod2.lme,random=pdDiag(~cyear))

anova(mod1.lme,mod2.lme,mod3.lme,mod1.lm)
         Model df       AIC      BIC    logLik   Test  L.Ratio p-value
mod1.lme     1  8  998.5859 1021.769 -491.2930                        
mod2.lme     2 10  992.2501 1021.229 -486.1250 1 vs 2 10.33588  0.0057
mod3.lme     3  9  991.6140 1017.694 -486.8070 2 vs 3  1.36392  0.2429
mod1.lm      4  7 1019.4313 1039.716 -502.7156 3 vs 4 31.81732  <.0001

##mod3.lme is the superior model
#non common slope and intercept, pos diagonal corr matrix

###noncommon slope is superior
summary(mod3.lme)
 Data: gmang 
      AIC      BIC   logLik
  991.614 1017.695 -486.807

Random effects:
 Formula: ~cyear | pl
 Structure: Diagonal
        (Intercept)    cyear Residual
StdDev:    6.819195 1.678588 6.770327

Fixed effects: sr ~ cyear + b5yr + yslb + log.ca. + grazer 
                Value Std.Error  DF   t-value p-value
(Intercept) 197.76317 28.086127 115  7.041312  0.0000
cyear         1.88623  0.491993 115  3.833845  0.0002
b5yr         -2.89965  1.046101 115 -2.771868  0.0065
yslb         -1.44957  0.569693 115 -2.544482  0.0123
log.ca.     -32.98744  8.245874 115 -4.000478  0.0001
grazer1      -1.58164  2.589583 115 -0.610770  0.5426

##diagnostic plots
plot(mod3.lme)
qqnorm(mod3.lme)
##both dianostic plots look reasonably good
plot(mod3.lme,resid(.)~fitted(.)|pl) ##plots do appear to have diff variances
plot(mod3.lme,resid(.,type='n')~fitted(.)|pl) ##plots do appear to have diff variances


mod4.lme<-update(mod3.lme,weights=varIdent(form=~1|pl))
anova(mod3.lme,mod4.lme)

plot(ACF(mod3.lme,maxLag=5),alpha=0.01)
varo<-Variogram(mod3.lme,maxDist=5,form=~cyear)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5),type='o',xlab='temporal lag',ylab='semivariance')


##there does not appear to be a large degree of autocorrelation

mod4.lme<-update(mod3.lme,corr=corARMA(p=1))
mod5.lme<-update(mod3.lme,corr=corARMA(p=2))
mod6.lme<-update(mod3.lme,corr=corARMA(p=3))
anova(mod3.lme,mod4.lme,mod5.lme,mod6.lme)
anova(mod3.lme,mod5.lme)
##autocorrelation functions do not appear necessary

plot(rep(1998:2004,20),predict(mod2.lme),type='n')
for(i in 1:20){
 points(c(1998:2004),predict(mod2.lme)[(7*(i-1)+1):(i*7)],col=i,pch=i,type='o')
# abline(lm(predict(mod2.lme)[(7*(i-1)+1):(i*7)]~c(1998:2004)),col=i)
}

plot(augPred(mod1.lme,primary=~cyear)) ##common slope
plot(augPred(mod2.lme,primary=~cyear)) ##non common slope
plot(augPred(mod2.lme,primary=~cyear,level = 0:1,length.out=2 ),xlab='year centered',ylab='species richness')

mod5.lme<-update(mod5.lme,.~.-grazer)
plot(augPred(mod5.lme,primary=~cyear,level = 0:1,length.out=2 ),xlab='year centered',ylab='species richness')

mod3.lme<-update(mod3.lme,.~.-grazer)
plot(augPred(mod3.lme,primary=~cyear)) ##non common slope


plot(augPred(mod2.lme,primary=~yslb,level = 0:1,length.out=2),xlab='years since last burn',ylab='species richness')
plot(augPred(mod2.lme,primary=~b5yr,level = 0:1,length.out=2),xlab='# burns in past 5 years',ylab='species richness')
plot(augPred(mod2.lme,primary=~log.ca.,level = 0:1,length.out=2),xlab='log calcium',ylab='species richness')

) ##non common slope

#examine lack of fit
plot(sr,predict(mod2.lme),type='n')
for(i in 1:20){
 points(sr[(7*(i-1)+1):(i*7)],predict(mod2.lme)[(7*(i-1)+1):(i*7)],col=i,pch=i)
}

###########################################
#############################################



plot(groupedData(sr~yr|plot + b5yr|plot,data=data))



interaction.plot(b5yr,season,sr)

summary(lm(sr~b5yr))

summary(m.mixed<-lme(sr~b5yr,random=~b5yr|plot))
plot(m.mixed)
xyplot(sr~b5yr|plot)
plot(groupedData(sr~b5yr|plot,data=data))

plot((sr~year+pl))
plot(mod<-lm(sr~year+pl))

summary(lm(residuals(mod)~yslb))
plot(residuals(mod)~yslb)
abline(lm(residuals(mod)~yslb))


test<-groupedData(sr~year|pl)
plot(test,outer=~grazer)

##################################################
mod1.lm<-lm(sr~1)
Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)   75.157      1.052   71.43   <2e-16 ***
Residual standard error: 12.45 on 139 degrees of freedom

mod2.lm<-lm(sr~pl-1)
Residual standard error: 9.155 on 120 degrees of freedom

mod1.lme<-lme(sr~1,random=~1|pl)
Linear mixed-effects model fit by REML
 Data: NULL 
       AIC      BIC    logLik
  1058.534 1067.337 -526.2668

Random effects:
 Formula: ~1 | pl
        (Intercept) Residual
StdDev:    8.624743 9.155404

Fixed effects: sr ~ 1 
               Value Std.Error  DF  t-value p-value
(Intercept) 75.15714  2.077988 120 36.16823       0

Standardized Within-Group Residuals:
         Min           Q1          Med           Q3          Max 
-2.124796027 -0.630523231  0.001783000  0.617962814  2.114921982 

Number of Observations: 140
Number of Groups: 20 

###Variation between plots is slightly less than variation within plots

mod2.lme<- lme(sr~yr,random=~1|pl)

mod3.lme<- lme(sr~pl,random=~1|yr)

mod4.lme<- lme(sr~pl,random=~1|yr/pl)

######################################################

interaction.plot(yr,pl,sr,las=1)
##the lines are not all parallel so there is likely an interaction between
##plot and year
mod1.lme<-lme(sr~yr,random=~1|pl)
##make new model that considers var of plots within years
mod2.lme<-update(mod1.lme,random=~1|yr/pl)
anova(mod1.lme,mod2.lme)

mod3.lme<-update(mod1.lme,random=~pl-1|yr)

#####################################################
gobject<-groupedData(sr~year|pl)
plot(gobject)

mod1.lm<-lm(sr~yr-1)
summary(mod1.lm)
Residual standard error: 11.5 on 133 degrees of freedom

gobject<-groupedData(sr~year|pl)
mod1.lis<-lmList(sr~year,data=gobject)
coef(mod1.lis)
plot(intervals(mod1.lis))
mod2.lis<-update(mod1.lis,sr~I(year-mean(year)))
plot(intervals(mod2.lis))

mod1.lme<-lme(sr~year,random=~1|pl,data=gobject)
summary(mod1.lme)
Linear mixed-effects model fit by REML
 Data: gobject 
       AIC      BIC    logLik
  1034.670 1046.379 -513.3352
Random effects:
 Formula: ~1 | pl
        (Intercept) Residual
StdDev:    8.754942 8.245101
##centering decreases the std error estimate for the incercept, and changes
#interpretation but nothing else

#new model will consider noncommon slopes
mod2.lme<-update(mod1.lme,random=~year|pl)

anova(mod1.lme,mod2.lme)
##a common slope seems better than noncommon slopes

random.effects(mod1.lme)
##or
ranef(mod1.lme)
plot(gobject,aspect='xy')
plot(augPred(mod1.lme),aspect='xy',grid=T)
###########################################
###Grazing investigation##################
yr<-year-1998
pl<-factor(plot)
contrasts(pl)<-contr.helmert

contrasts(grazer)<-c(1,0)
colnames(contrasts(grazer))<-c('bis')


mod1.lm<-lm(sr~year+pl+grazer)
summary(mod1.lm)
anova(mod1.lm)

mod2.lm<-lm(sr~year+pl+grazer+log.ca.)
summary(mod2.lm)
anova(mod2.lm)

srb<-sr[bbal==1&pl!=343]
yrb<-yr[bbal==1&pl!=343]
plb<-pl[bbal==1&pl!=343]
grb<-grazer[bbal==1&pl!=343]
interaction.plot(yrb,grb,srb)
interaction.plot(yr,grazer,sr)

summary(lm(srb~yrb*grb))

gobject<-groupedData(srb~grb|plb)



mod1.lis<-lmList(sr~grazer|pl,data=gobject)

mod1.lme<-lme(sr~grazer*yr,random=~1|pl)
summary(mod1.lme)
interaction.plot(yr, grazer, sr)

plot(gobject,outer=~grazer)
mod1.lme<-lme(sr~grazer,random=~1|yr,data=gobject)
plot(augPred(mod1.lme),aspect='xy',grid=T)


##########################################
#######autocorrelation functions##########
#from previous analysis (above) common slope was best idea
pl<-factor(plot)
contrasts(pl)<-contr.helmert
cyear<-year-mean(year)
gobject<-groupedData(sr~cyear|pl)
plot(gobject,inner=~b5yr)
mod1.lme<-lme(sr~cyear,random=~cyear|pl)
ACF( mod1.lme)
plot(ACF(mod1.lme,maxLag=5),alpha=0.01)

mod2.lme<-update(mod1.lme,corr=corAR1())
anova(mod1.lme,mod2.lme)
##incorperating 1st lv autoregression model sig. improved fit

mod3.lme<-update(mod1.lme,corr=corARMA(q=2))

anova(mod2.lme,mod3.lme,test=F)

mod4.lme<-update(mod1.lme,corr=corARMA(p=1,q=1))
anova(mod2.lme,mod3.lme,mod4.lme,test=F)

plot(ACF(mod3.lme,maxLag=5,resType="normalized"),alpha=.01)
##model 3 with a 2nd order moving average was best

####instead of correlation functions try out variograms
varo<-Variogram(mod1.lme,form=~cyear)
plot(varo[,2],varo[,1],xlim=c(0,5),ylim=c(0,1.5))
lines(lowess(varo[,2],varo[,1]))

mod5.lme<-update(mod1.lme,corr=corExp(form=~cyear,nugget=T))
plot(Variogram(mod5.lme,form=~cyear,maxDist=5))
plot(Variogram(mod5.lme,form=~cyear,maxDist=5,resType='n'))

anova(mod1.lme,mod5.lme,mod3.lme)


mod6.lme<-update(mod3.lme,weights=varPower())
anova(mod3.lme,mod6.lme)
##adjusting the variance function is not necessary
##analysis may not be proper b/c corr structure and var comp compete with each other

mod1.gls<-gls(sr~cyear,corr=corARMA(q=2))

#######################

##########################################
detach(data)
rich<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/LTREB proposal/rich(98-07).csv',sep=',',header=T)
names(rich)
[1] "plot"    "yr"      "sr"      "plot.n"  "yr.n"    "plot.yr" "names"
attach(rich)

pl<-factor(names)
contrasts(pl)<-contr.helmert

cyear<-yr.n-mean(yr.n)

mod1.aov<-aov(sr~cyear+Error(pl))
summary(mod1.aov)
mod1.lm<-lm(sr~cyear)
summary(mod1.lm)
Coefficients:
            Estimate Std. Error t value Pr(>|t|)    
(Intercept)  76.1100     0.8895  85.564  < 2e-16 ***
cyear         1.0400     0.3097   3.358  0.00094 ***
Residual standard error: 12.58 on 198 degrees of freedom

#library(lattice) #needed for bwplot function
gobject<-groupedData(sr~cyear|pl)

bwplot(getGroups(gobject)~resid(mod1.lm))
anova(mod1.lm)
plot(yr,sr)
abline(lm(sr~yr))

plot(gobject)
mod1.lis<-lmList(gobject)
summary(mod1.lis)
Residual standard error: 7.786794 on 160 degrees of freedom
plot(intervals(mod1.lis))
bwplot(getGroups(gobject)~resid(mod1.lis))
pairs(mod1.lis,id=0.01)

mod1.lme<-lme(sr~cyear,random=~1|pl,data=gobject)
mod2.lme<-update(mod1.lme,random=~cyear|pl)
mod3.lme<-update(mod2.lme,random=pdDiag(~cyear))
anova(mod1.lme,mod2.lme,mod3.lme,mod1.lm)
#mod2 is best fit, but mod 3 with assumption of independence between plots is as good (p=0.07)
summary(mod3.lme)
plot(Variogram(mod3.lme,form=~cyear))
plot(ACF(mod3.lme,maxLag=5),alpha=0.01)

mod4.lme<-update(mod3.lme,corr=corAR1())
anova(mod3.lme,mod4.lme)


plot(augPred(mod3.lme))
plot(mod3.lme)
plot(cyear,sr)
abline(fixef(mod3.lme))
qqnorm(mod3.lme,abline=c(0,1))
# normal plot of standardized residuals by plot
qqnorm(mod3.lme, ~ resid(., type = "p") | pl, abline = c(0, 1))
# normal plots of random effects
qqnorm(mod3.lme, ~ranef(.))

plot(compareFits(coef(mod1.lis),coef(mod3.lme)),mark=fixef(mod3.lme))
plot(comparePred(mod1.lis,mod3.lme,length.out=2),mark=fixef(mod3.lme))
plot(comparePred(mod2.lme,mod3.lme,length.out=2),mark=fixef(mod3.lme))
bwplot(getGroups(gobject)~resid(mod3.lme))
bwplot(cyear~resid(mod3.lme),xlim=c(-32,32))
bwplot(cyear~resid(mod1.lm))

