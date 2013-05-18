library(vegan)

setwd('~/tgp_management/')

source('./scripts/tgp_functions.R')
source('./scripts/tgp_grid_data_import.R')

dir.create('./results/')

## test significance of model terms
print('OLS model of richness------------------------------------')
## ols model of richness
full = rda(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
soil_ind = rda(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
mang_ind = rda(env$sr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

summary(lm(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB))

anova(full, by='margin')
anova(soil_ind, by='margin')
anova(mang_ind, by='margin')

print('RDA model of composition-------------------------------')
### rda model of composition
full = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
soil_ind = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
mang_ind = rda(comm_sqr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

anova(full, by='margin')
anova(soil_ind, by='margin')
anova(mang_ind, by='margin')

print('CCA model of composition-------------------------------')
### cca model of composition
full = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
soil_ind = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
rain_ind = cca(comm_sqr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
mang_ind = cca(comm_sqr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

anova(full, by='margin')
anova(soil_ind, by='margin')
anova(mang_ind, by='margin')


## examine for residual spatial dependence
tgp_xy = env[ , c('easting', 'northing')]

ols_lm = rda(env$sr, cbind(soil_mat, mang_mat))
rda_lm = rda(comm_sqr, cbind(soil_mat, mang_mat))
cca_lm = cca(comm_sqr, cbind(soil_mat, mang_mat))

ols_mso = mso(ols_lm, tgp_xy, grain=1000, permutations = 999)
rda_mso = mso(rda_lm, tgp_xy, grain=1000, permutations = 999)
cca_mso = mso(cca_lm, tgp_xy, grain=1000, permutations = 999)

pdf('./figs/grid_model_mso.pdf')
  msoplot(ols_mso, main='Grid, OLS')
  msoplot(rda_mso, main='Grid, RDA')
  msoplot(cca_mso, main='Grid, CCA')
dev.off()
