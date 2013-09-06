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

summary(lm(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB))

anova(full, by='margin')

soil = rda(env$sr, soil_mat, cbind(rain_mat, mang_mat))
permutest(soil, permutations=999)

rain = rda(env$sr, rain_mat, cbind(soil_mat, mang_mat))
permutest(rain, permutations=999)

mang = rda(env$sr, mang_mat, cbind(rain_mat, soil_mat))
permutest(mang, permutations=999)

print('RDA model of composition-------------------------------')
### rda model of composition
full = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

soil = rda(comm_sqr, soil_mat, cbind(rain_mat, mang_mat))
permutest(soil, permutations=999)

rain = rda(comm_sqr, rain_mat, cbind(soil_mat, mang_mat))
permutest(rain, permutations=999)

mang = rda(comm_sqr, mang_mat, cbind(rain_mat, soil_mat))
permutest(mang, permutations=999)


print('CCA model of composition-------------------------------')
### cca model of composition
full = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             rain_mat$sum_rain + rain_mat$win_rain + rain_mat$spr_rain +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

soil = cca(comm_sqr, soil_mat, cbind(rain_mat, mang_mat))
permutest(soil, permutations=999)

rain = cca(comm_sqr, rain_mat, cbind(soil_mat, mang_mat))
permutest(rain, permutations=999)

mang = cca(comm_sqr, mang_mat, cbind(rain_mat, soil_mat))
permutest(mang, permutations=999)

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

pdf('./figs/repeat_model_mso.pdf', width=7*2, height=7)
  par(mfrow=c(1,2))
  msoplot(ols_mso_spat, main='Repeat, OLS, spatial', ylim=c(50, 250))
  msoplot(ols_mso_temp, main='Repeat, OLS, temporal', ylim=c(50, 225))
  msoplot(rda_mso_spat, main='Repeat, RDA, spatial', ylim=c(10, 40))
  msoplot(rda_mso_temp, main='Repeat, RDA, temporal', ylim=c(20, 35))
  msoplot(cca_mso_spat, main='Repeat, CCA, spatial', ylim=c(1, 4))
  msoplot(cca_mso_temp, main='Repeat, CCA, temporal', ylim=c(2, 3.5))
dev.off()
