library(vegan)

setwd('~/tgp_management/')

source('./scripts/tgp_functions.R')
source('./scripts/tgp_repeat_data_import.R')

dir.create('./results/')

rain_mat = as.data.frame(rain_mat)
## test significance of model terms
print('OLS model of richness-----------------------------------')
## ols model of richness
full = rda(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
                    rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
                    mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
soil_ind = rda(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
               Condition(rain_mat$sum_rain) + Condition(rain_mat$win_rain) + Condition(rain_mat$spr_rain) +
               Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
rain_ind = rda(env$sr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
mang_ind = rda(env$sr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 Condition(rain_mat$sum_rain) + Condition(rain_mat$win_rain) + Condition(rain_mat$spr_rain) +
                 mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

summary(lm(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB))

anova(full, by='margin')
anova(soil_ind, by='margin')
anova(rain_ind, by='margin')
anova(mang_ind, by='margin')

print('RDA model of composition-------------------------------')
### rda model of composition
full = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
soil_ind = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
                 Condition(rain_mat$sum_rain) + Condition(rain_mat$win_rain) + Condition(rain_mat$spr_rain) +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
rain_ind = rda(comm_sqr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
mang_ind = rda(comm_sqr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 Condition(rain_mat$sum_rain) + Condition(rain_mat$win_rain) + Condition(rain_mat$spr_rain) +
                 mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

anova(full, by='margin')
anova(soil_ind, by='margin')
anova(rain_ind, by='margin')
anova(mang_ind, by='margin')

print('CCA model of composition-------------------------------')
### cca model of composition
full = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
soil_ind = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
                 Condition(rain_mat$sum_rain) + Condition(rain_mat$win_rain) + Condition(rain_mat$spr_rain) +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
rain_ind = cca(comm_sqr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
                 Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) + Condition(mang_mat$YrsSLB))
mang_ind = cca(comm_sqr ~ 
                 Condition(soil_mat[,1]) + Condition(soil_mat[,2]) + Condition(soil_mat[,3]) +
                 Condition(rain_mat$sum_rain) + Condition(rain_mat$win_rain) + Condition(rain_mat$spr_rain) +
                 mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

anova(full, by='margin')
anova(soil_ind, by='margin')
anova(rain_ind, by='margin')
anova(mang_ind, by='margin')


## examine for residual spatial dependence
tgp_xy = env[ , c('easting', 'northing')]

ols_lm = rda(env$sr, cbind(soil_mat, rain_mat, mang_mat))
rda_lm = rda(comm_sqr, cbind(soil_mat, rain_mat, mang_mat))
cca_lm = cca(comm_sqr, cbind(soil_mat, rain_mat, mang_mat))

ols_mso_spat = mso(ols_lm, tgp_xy, grain=1000, permutations = 999)
rda_mso_spat = mso(rda_lm, tgp_xy, grain=1000, permutations = 999)
cca_mso_spat = mso(cca_lm, tgp_xy, grain=1000, permutations = 999)

ols_mso_temp = mso(ols_lm, env$yr, permutations = 999)
rda_mso_temp = mso(rda_lm, env$yr, permutations = 999)
cca_mso_temp = mso(cca_lm, env$yr, permutations = 999)

pdf('./figs/repeat_model_mso.pdf', width=7, height=7*2)
  par(mfrow=c(2,1))
  msoplot(ols_mso_spat, main='Repeat, OLS, spatial')
  msoplot(ols_mso_temp, main='Repeat, OLS, temporal')
  msoplot(rda_mso_spat, main='Repeat, RDA, spatial')
  msoplot(rda_mso_temp, main='Repeat, RDA, temporal')
  msoplot(cca_mso_spat, main='Repeat, CCA, spatial')
  msoplot(cca_mso_temp, main='Repeat, CCA, temporal')
dev.off()
