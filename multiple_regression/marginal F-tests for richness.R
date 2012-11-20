library(nlme)

summary(lm(rich~name+yr.f+dist.m))
mod<-gls(rich~name+yr.f+dist.m,meth="ML")
anova(mod,type='m')##marginal F-tests
Denom. DF: 187 
            numDF  F-value p-value
(Intercept)     1 468.3793  <.0001
name           19  17.6931  <.0001
yr.f           10   9.1610  <.0001
dist.m          3   9.3759  <.0001
summary(lm(rich~name+yr.f+dist.m))
F-statistic:    20 on 32 and 187 DF,  p-value: < 2.2e-16

mod.ar1<-update(mod,.~.,cor=corAR1(form=~yr|name))
anova(mod.ar1,type='m')
Denom. DF: 187 
            numDF  F-value p-value
(Intercept)     1 389.9952  <.0001
name           19  13.4791  <.0001
yr.f           10  10.2736  <.0001
dist.m          3   7.2492   1e-04

plot(ACF(mod.ar1,maxLag=5,resType='r'),alpha=0.05)
plot(ACF(mod.ar1,maxLag=5,resType='n'),alpha=0.05)

##understanding F-value
qf(1-1e-4, 3, 220-3-1)
1-pf(7.2492, 3, 220-3-1)
##I would like to calculate a F-value for the whole model




mod.int<-gls(rich~1,meth="ML")
anova(mod.int,mod,mod.ar1)
        Model df      AIC      BIC    logLik   Test  L.Ratio p-value
mod.int     1  2 1749.165 1755.952 -872.5824                        
mod         2 34 1486.125 1601.508 -709.0625 1 vs 2 327.0397  <.0001
mod.ar1     3 35 1483.096 1601.873 -706.5481 2 vs 3   5.0287  0.0249


