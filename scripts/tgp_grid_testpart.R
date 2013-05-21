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
summary(lm(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB))

anova(full, by='margin')

print('RDA model of composition-------------------------------')
### rda model of composition
full = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

print('CCA model of composition-------------------------------')
### cca model of composition
full = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

anova(full, by='margin')

## examine for residual spatial dependence
tgp_xy = env[ , c('easting', 'northing')]

ols_lm = rda(env$sr, cbind(soil_mat, mang_mat))
rda_lm = rda(comm_sqr, cbind(soil_mat, mang_mat))
cca_lm = cca(comm_sqr, cbind(soil_mat, mang_mat))

ols_mso = mso(ols_lm, tgp_xy, grain=1000, permutations = 999)
rda_mso = mso(rda_lm, tgp_xy, grain=1000, permutations = 999)
cca_mso = mso(cca_lm, tgp_xy, grain=1000, permutations = 999)

pdf('./figs/grid_model_mso.pdf')
  msoplot(ols_mso, main='Grid, OLS', ylim=c(50, 250))
  msoplot(rda_mso, main='Grid, RDA', ylim=c(15, 30))
  msoplot(cca_mso, main='Grid, CCA', ylim=c(3, 6))
dev.off()
