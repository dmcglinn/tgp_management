

varML<-function(x){
 sum( (x- 75.11689)^2)/length(x)
}

Neff<-function(n,rho){
 n*(1-rho)/(1+rho)
}

sqrt(varML(rich)/Neff(220,.640706))
sd(rich)/sqrt(Neff(222,0.6464833))

summary(test3)
summary(test2)

(sd(rich)/1.628694)^2
(Neff(222,0.6464833))

site.m<-endatlv1[,8:26]
armod<-arima(x=rich,order=c(1,0,0),xreg=site.m)
tsdiag(armod)
Box.test(rich)

##simulate time series
set.seed(1)
x=1:100
y=x+arima.sim(n = 100, list(ar = c(0.1)))
plot(y,type='l')
arima(y,order=c(1,0,0))
arima(y,order=c(1,0,0),xreg=x)
tsdiag(arima(y,order=c(1,0,0),xreg=x))

rich.c<-rich-mean(rich)
((t(rich.c)%*%rich.c)*(220-1)^-1*(220)^-1)^.5
sd(rich)/sqrt(220)
summary(test)

sum((rich-mean(rich))^2)
var(rich)
t(rich.c)%*%rich.c


dist.m<-cbind(YrsOB,YrsSLB,BP5Yrs)
mod<-gls(rich~name+yr.f+dist.m)
modar<-gls(rich~name+yr.f+dist.m,cor=corAR1(form=~yr|name))
anova(mod,type='m')
Denom. DF: 187 
            numDF  F-value p-value
(Intercept)     1 551.0345  <.0001
name           19  20.8154  <.0001
yr.f           10  10.7777  <.0001
dist.m          3  11.0305  <.0001
anova(modar,type='m')
Denom. DF: 187 
            numDF  F-value p-value
(Intercept)     1 365.3284  <.0001
name           19  11.5282  <.0001
yr.f           10  13.1652  <.0001
dist.m          3   6.4018   4e-04
summary(mod)

library(vegan)
yr.f
srs
test<-rda(srs[,-1]~YrsOB+YrsSLB+BP5Yrs+Condition(name)+Condition(yr.f))
anova(test,by='m')
Permutation test for rda under reduced model
Marginal effects of terms
Model: rda(formula = srs[, -1] ~ YrsOB + YrsSLB + BP5Yrs + Condition(name) + Condition(yr.f))
          Df     Var       F  N.Perm Pr(>F)   
YrsOB      1   1.325  8.8495 199.000  0.005 **
YrsSLB     1   1.831 12.2286 199.000  0.005 **
BP5Yrs     1   0.965  6.4454 199.000  0.010 **
Residual 187  28.004                          
anova(test,by='m',strata=yr)
coef(test)

plot(test,display = c("sp","cn"),scaling=1)

test<-rda(srs[,1]~YrsOB+YrsSLB+BP5Yrs+Condition(name)+Condition(yr.f))
plot(test)
anova(test,by='m')
Model: rda(formula = srs[, 1] ~ YrsOB + YrsSLB + BP5Yrs + Condition(name) + Condition(yr.f))
          Df     Var      F  N.Perm Pr(>F)   
YrsOB      1   4.080 20.585 199.000  0.005 **
YrsSLB     1   1.791  9.038 199.000  0.010 **
BP5Yrs     1   1.483  7.480 199.000  0.015 * 
Residual 187  37.066                              
anova(test,by='m',strata='yr')
Model: rda(formula = srs[, 1] ~ YrsOB + YrsSLB + BP5Yrs + Condition(name) + Condition(yr.f))
          Df     Var      F  N.Perm Pr(>F)   
YrsOB      1   4.080 20.585 199.000  0.005 **
YrsSLB     1   1.791  9.038 199.000  0.010 **
BP5Yrs     1   1.483  7.480 199.000  0.005 **
Residual 187  37.066                         
anova(test,by='m',strata='name')
Model: rda(formula = srs[, 1] ~ YrsOB + YrsSLB + BP5Yrs + Condition(name) + Condition(yr.f))
          Df     Var      F  N.Perm Pr(>F)   
YrsOB      1   4.080 20.585 199.000  0.005 **
YrsSLB     1   1.791  9.038 199.000  0.015 * 
BP5Yrs     1   1.483  7.480 199.000  0.010 **
Residual 187  37.066


