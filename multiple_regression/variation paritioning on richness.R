endat<-read.table('C:/Users/Daniel McGlinn/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env corners.csv',sep=',',header=T)
names(endat)
attach(endat)
##work only with level 1
endatlv1<-endat[lv==1&cor==1,] ##to only work at the 100 m2 scale
detach(endat)
attach(endatlv1)
yr.f<-factor(yr)

##for var part with package:vegan##
site.m<-endatlv1[,8:27]
yr.m<-endatlv1[,28:38]
dist.m<-cbind(YrsOB,bison,YrsSLB,BP5Yrs)

library(vegan)
varpar<-varpart(rich,site.m,yr.m,dist.m)
Partition of variation in RDA
Call:
varpart(Y = rich, X = site.m, yr.m, dist.m) 
Explanatory tables:
X1:  site.m
X2:  yr.m
X3:  dist.m 
No. of explanatory tables: 3 
Total variation (SS): 35894 
            Variance: 163.90 
No. of observations: 220 

Partition table:
                      Df R.square Adj.R.square Testable
[a+d+f+g] = X1        19  0.55536      0.51312     TRUE
[b+d+e+g] = X2        10  0.17846      0.13916     TRUE
[c+e+f+g] = X3         4  0.15942      0.14378     TRUE
[a+b+d+e+f+g] = X1+X2 29  0.73383      0.69320     TRUE
[a+c+d+e+f+g] = X1+X3 23  0.64351      0.60168     TRUE
[b+c+d+e+f+g] = X2+X3 14  0.30158      0.25389     TRUE
[a+b+c+d+e+f+g] = All 33  0.77387      0.73375     TRUE
Individual fractions                                   
[a] = X1 | X2+X3      19               0.47986     TRUE
[b] = X2 | X1+X3      10               0.13207     TRUE
[c] = X3 | X1+X2       4               0.04055     TRUE
[d]                    0              -0.02197    FALSE
[e]                    0               0.04801    FALSE
[f]                    0               0.07419    FALSE
[g]                    0              -0.01896    FALSE
[h] = Residuals                        0.26625    FALSE


plot(varpart(rich,site.m,yr.m,dist.m))
##this method does not completely agree with the manual method below
##this may due to subtle diffs in MR and RDA##

###Variation paritioning
##Legendre p205... to get the unadjusted components..
mod.full<-lm(rich~name+yr.f+dist.m)
mod.pl<-lm(rich~name)
mod.yr<-lm(rich~yr.f)
mod.ds<-lm(rich~dist.m)
mod.plyr<-lm(rich~name+yr.f)
mod.plds<-lm(rich~name+dist.m)
mod.yrds<-lm(rich~yr.f+dist.m)

##plot only
r2.diff<-function(modfull,modsub){
  summary(modfull)$r.sq-summary(modsub)$r.sq
}

r2.full<-summary(mod.full)$r.sq
pl.r2<-r2.diff(mod.full,mod.yrds)
yr.r2<-r2.diff(mod.full,mod.plds)
ds.r2<-r2.diff(mod.full,mod.plyr)

##shared components
plyr.r2<-r2.diff(mod.full,mod.ds)-pl.r2-yr.r2
plds.r2<-r2.diff(mod.full,mod.yr)-pl.r2-ds.r2
yrds.r2<-r2.diff(mod.full,mod.pl)-yr.r2-ds.r2
plyrds.r2<-r2.full-pl.r2-yr.r2-ds.r2-plyr.r2-plds.r2-yrds.r2

r2.diff.adj<-function(modfull,modsub){
 summary(modfull)$adj-summary(modsub)$adj
}

r2.full.adj<-summary(mod.full)$adj
pl.r2.adj<-r2.diff.adj(mod.full,mod.yrds)
yr.r2.adj<-r2.diff.adj(mod.full,mod.plds)
ds.r2.adj<-r2.diff.adj(mod.full,mod.plyr)

##shared components
plyr.r2.adj<-r2.diff.adj(mod.full,mod.ds)-pl.r2.adj-yr.r2.adj
plds.r2.adj<-r2.diff.adj(mod.full,mod.yr)-pl.r2.adj-ds.r2.adj
yrds.r2.adj<-r2.diff.adj(mod.full,mod.pl)-yr.r2.adj-ds.r2.adj
plyrds.r2.adj<-r2.full.adj-pl.r2.adj-yr.r2.adj-ds.r2.adj-plyr.r2.adj-plds.r2.adj-yrds.r2.adj

parts<-c(r2.full,pl.r2,yr.r2,ds.r2,plyr.r2,yrds.r2,plds.r2,plyrds.r2)
parts.adj<-c(r2.full.adj,pl.r2.adj,yr.r2.adj,ds.r2.adj,plyr.r2.adj,yrds.r2.adj,plds.r2.adj,plyrds.r2.adj)
names<-c('full','pl','yr','ds','plyr','plds','yrds','plyrds')
data.frame(names=names,parts=parts,adj=parts.adj)
   names       parts         adj
1   full  0.77386843  0.73374832
2     pl  0.47228369  0.47986023
3     yr  0.13035831  0.13207119
4     ds  0.04004118  0.04054743
5   plyr  0.01180421 -0.02196667
6   plds  0.04810497  0.04800614
7   yrds  0.08308028  0.07418535
8 plyrds -0.01180421 -0.01895537

###Drop binary bison variable
dist.m<-cbind(YrsOB,BP5Yrs,YrsSLB)
varpar<-varpart(rich,site.m,yr.m,dist.m)
No. of explanatory tables: 3 
Total variation (SS): 35894 
            Variance: 163.90 
No. of observations: 220 

Partition table:
                      Df R.square Adj.R.square Testable
[a+d+f+g] = X1        19  0.55536      0.51312     TRUE
[b+d+e+g] = X2        10  0.17846      0.13916     TRUE
[c+e+f+g] = X3         3  0.15470      0.14296     TRUE
[a+b+d+e+f+g] = X1+X2 29  0.73383      0.69320     TRUE
[a+c+d+e+f+g] = X1+X3 22  0.64350      0.60369     TRUE
[b+c+d+e+f+g] = X2+X3 13  0.29555      0.25109     TRUE
[a+b+c+d+e+f+g] = All 32  0.77385      0.73515     TRUE
Individual fractions                                   
[a] = X1 | X2+X3      19               0.48405     TRUE
[b] = X2 | X1+X3      10               0.13145     TRUE
[c] = X3 | X1+X2       3               0.04195     TRUE
[d]                    0              -0.02332    FALSE
[e]                    0               0.04862    FALSE
[f]                    0               0.06999    FALSE
[g]                    0              -0.01760    FALSE
[h] = Residuals                        0.26485    FALSE

plot(varpar)

###Variation paritioning w/o Bison binary
###Variation paritioning
##Legendre p205... to get the unadjusted components..
mod.full<-lm(rich~name+yr.f+dist.m)
mod.pl<-lm(rich~name)
mod.yr<-lm(rich~yr.f)
mod.ds<-lm(rich~dist.m) ##note that 3 auto corr terms are appr here not just 1
mod.plyr<-lm(rich~name+yr.f)
mod.plds<-lm(rich~name+dist.m)
mod.yrds<-lm(rich~yr.f+dist.m)

##plot only
r2.diff<-function(modfull,modsub){
  summary(modfull)$r.sq-summary(modsub)$r.sq
}

r2.full<-summary(mod.full)$r.sq
pl.r2<-r2.diff(mod.full,mod.yrds)
yr.r2<-r2.diff(mod.full,mod.plds)
ds.r2<-r2.diff(mod.full,mod.plyr)

##shared components
plyr.r2<-r2.diff(mod.full,mod.ds)-pl.r2-yr.r2
plds.r2<-r2.diff(mod.full,mod.yr)-pl.r2-ds.r2
yrds.r2<-r2.diff(mod.full,mod.pl)-yr.r2-ds.r2
plyrds.r2<-r2.full-pl.r2-yr.r2-ds.r2-plyr.r2-plds.r2-yrds.r2

r2.diff.adj<-function(modfull,modsub){
 summary(modfull)$adj-summary(modsub)$adj
}
summary(mod.full)$r.s
summary(mod.full)$adj
summary(mod.yrds)$r.s
summary(mod.yrds)$adj

r2.full.adj<-summary(mod.full)$adj
pl.r2.adj<-r2.diff.adj(mod.full,mod.yrds)
yr.r2.adj<-r2.diff.adj(mod.full,mod.plds)
ds.r2.adj<-r2.diff.adj(mod.full,mod.plyr)

##shared components
plyr.r2.adj<-r2.diff.adj(mod.full,mod.ds)-pl.r2.adj-yr.r2.adj
plds.r2.adj<-r2.diff.adj(mod.full,mod.yr)-pl.r2.adj-ds.r2.adj
yrds.r2.adj<-r2.diff.adj(mod.full,mod.pl)-yr.r2.adj-ds.r2.adj
plyrds.r2.adj<-r2.full.adj-pl.r2.adj-yr.r2.adj-ds.r2.adj-plyr.r2.adj-plds.r2.adj-yrds.r2.adj

parts<-c(r2.full,pl.r2,yr.r2,ds.r2,plyr.r2,yrds.r2,plds.r2,plyrds.r2)
parts.adj<-c(r2.full.adj,pl.r2.adj,yr.r2.adj,ds.r2.adj,plyr.r2.adj,yrds.r2.adj,plds.r2.adj,plyrds.r2.adj)
names<-c('full','pl','yr','ds','plyr','plds','yrds','plyrds')
data.frame(names=names,parts=parts,adj=parts.adj)
   names       parts         adj
1   full  0.77384715  0.73514720
2     pl  0.47829906  0.48405481
3     yr  0.13034269  0.13145442
4     ds  0.04001990  0.04194632
5   plyr  0.01050320 -0.02332399
6   plds  0.04812059  0.04862292
7   yrds  0.07706490  0.06999076
8 plyrds -0.01050320 -0.01759804


##test a few components
pl.rda<-rda(rich~name+Condition(yr.f)+Condition(dist.m))
yr.rda<-rda(rich~yr.f+Condition(name)+Condition(dist.m))
ds.rda<-rda(rich~YrsOB+YrsSLB+BP5Yrs+Condition(name)+Condition(yr.f))

anova(pl.rda,strata=yr,step=400)
Model: rda(formula = rich ~ name + Condition(yr.f) + Condition(dist.m))
          Df    Var      F N.Perm Pr(>F)   
Model     19  78.39 20.815 399.00 0.0025 **
Residual 187  37.07  

##F value can be hand calculated as so
(78.392/19)/(37.07/(220-19-10-3-1)) ##MSE/RSE

anova(yr.rda,strata=name,step=400)
Model: rda(formula = rich ~ yr.f + Condition(name) + Condition(dist.m))
          Df    Var      F N.Perm Pr(>F)   
Model     10  21.36 10.778 399.00 0.0025 **
Residual 187  37.07

anova(ds.rda,strata=name,step=400)
Model: rda(formula = rich ~ dist.m + Condition(name) + Condition(yr.f))
          Df    Var      F N.Perm Pr(>F)   
Model      3   6.56 11.030 399.00 0.0025 **
Residual 187  37.07
anova(ds.rda,strata=yr,step=400)
##same result as above
anova(ds.rda,strata=name,step=2000,by='margin')
Model: rda(formula = rich ~ YrsOB + YrsSLB + BP5Yrs + Condition(name) + Condition(yr.f))
YrsOB      1    4.08 20.585 1999.00  5e-04 ***
YrsSLB     1    1.79  9.038 1999.00  4e-03 ** 
BP5Yrs     1    1.48  7.480 1999.00  7e-03 **
Residual 187   37.07 
anova(ds.rda,strata=yr,step=2000,by='margin')
          Df     Var      F  N.Perm Pr(>F)    
YrsOB      1    4.08 20.585 1999.00 0.0010 ***
YrsSLB     1    1.79  9.038 1999.00 0.0035 ** 
BP5Yrs     1    1.48  7.480 1999.00 0.0080 **

anova(rda(rich~name+yr.f+YrsOB+YrsSLB+BP5Yrs),by='margin',step=2000)

r2adj<-function(r2,p){1-((220-1)/(220-p-1)*(1-r2))}
##Code below is obsolete, varpart does correspond with the above approach
##################################################################
###why does varpart with rda not equal my approach?
##because varpart calculates diffs in R2 values to get 
##indvidual components rather than regressions on residuals

set.seed(12)
x1<-rnorm(100)
x2<-rnorm(100)+.2*x1
y<-.5*x1+.9*x2+runif(100)
y<-y-mean(y)
plot(varpar<-varpart(y,x1,x2))

ss.func<-function(mod){
 sum((mod$fit-mean(mod$fit))^2)
}
tot.ss<-sum((y-mean(y))^2)
modx1<-lm(y~x1)
modx2<-lm(y~x2)
mod.full<-lm(y~x1+x2)
x1.ss<-ss.func(modx1)
x2.ss<-ss.func(modx2)
full.ss<-ss.func(mod.full)
x1.ss/tot.ss
x2.ss/tot.ss
full.ss/tot.ss
r2adj<-function(reg.ss,p){1-((100-1)/(100-p-1)*(1-reg.ss/tot.ss))}
r2adj(x1.ss,1)
r2adj(x2.ss,1)
r2adj(full.ss,2)
r2adj(full.ss,2)-r2adj(x2.ss,1)##individual comp of x1
r2adj(full.ss,2)-r2adj(x1.ss,1)##indiv comp of x2 
r2adj(full.ss,2)-(r2adj(full.ss,2)-r2adj(x2.ss,1))-(r2adj(full.ss,2)-r2adj(x1.ss,1))##shared compon
varpar$par
##all above agress with rda varpart
x1.ss.ind<-ss.func(lm(resid(modx2)~x1))
x2.ss.ind<-ss.func(lm(resid(modx1)~x2))
x1.ss.ind/tot.ss
x2.ss.ind/tot.ss
r2adj(x1.ss.ind,1)
r2adj(x2.ss.ind,1)
#r2adj((full.ss/tot.ss)-(x1.ss.ind/tot.ss)-(x2.ss.ind/tot.ss),1)

y.c<-y-mean(y)
(y.c-mod.full$fit)%*%(y.c-mod.full$fit)
sqrt(sum((y-mean(mod.full$fit))^2))

vegandocs(doc = c("NEWS", "ChangeLog", "FAQ-vegan.pdf",
    "intro-vegan.pdf", "diversity-vegan.pdf", "decision-vegan.pdf",
    "partitioning.pdf"))

vegandocs(doc = "decision-vegan.pdf")



