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

resid_list = list(residuals(ols_lm), residuals(rda_lm))

## Mantel tests
spat_mantel = sapply(1:2, function(x) mantel(dist(tgp_xy), dist(resid_list[[x]])))

temp_mantel = sapply(1:2, function(x) mantel(dist(env$yr), dist(resid_list[[x]])))

mod_names = c('OLS', 'RDA')
xlabs = c('Spatial Distance (m)', 'Temporal Distance (yr)')

pdf('./figs/repeat_model_mantel.pdf', width=7*2, height=7)
  par(mfrow=c(1,2))
  for(i in 1:2) {
    plot(dist(tgp_xy), dist(resid_list[[i]]), ylab='Residual Distance', xlab=xlabs[1],
         main=paste('Repeat, ', mod_names[i], ', Spatial Corr', sep=''),
         type='n')
    abline(lm(dist(resid_list[[i]]) ~ dist(tgp_xy)), col='red', lwd=2)
    lines(lowess(dist(tgp_xy), dist(resid_list[[i]])), col='blue', lwd=2)
    mtext(side=3,paste('Mantel, p=', spat_mantel[4, i], sep=''))
    ##
    plot(dist(env$yr), dist(resid_list[[i]]), ylab='Residual Distance', xlab=xlabs[2],
         main=paste('Repeat, ', mod_names[i], ', Temporal Corr', sep=''),
         type='n')
    abline(lm(dist(resid_list[[i]]) ~ dist(env$yr)), col='red', lwd=2)
    lines(lowess(dist(env$yr), dist(resid_list[[i]])), col='blue', lwd=2)
    mtext(side=3,paste('Mantel, p=', temp_mantel[4, i], sep=''))
    legend('topright', c('Linear', 'Lowess'), col=c('red', 'blue'), lty=1, 
           lwd=4, bty='n')    
  }  
dev.off()

