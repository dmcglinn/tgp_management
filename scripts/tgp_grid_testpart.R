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

soil = rda(env$sr, soil_mat, mang_mat)
permutest(soil, permutations=999)

mang = rda(env$sr, mang_mat, soil_mat)
permutest(mang, permutations=999)

print('RDA model of composition-------------------------------')
### rda model of composition
full = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

soil = rda(comm_sqr, soil_mat, mang_mat)
permutest(soil, permutations=999)

mang = rda(comm_sqr, mang_mat, soil_mat)
permutest(mang, permutations=999)

print('CCA model of composition-------------------------------')
### cca model of composition
full = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

soil = cca(comm_sqr, soil_mat, mang_mat)
permutest(soil, permutations=999)

mang = cca(comm_sqr, mang_mat, soil_mat)
permutest(mang, permutations=999)

## examine for residual spatial dependence
tgp_xy = env[ , c('easting', 'northing')]

ols_lm = rda(env$sr, cbind(soil_mat, mang_mat))
rda_lm = rda(comm_sqr, cbind(soil_mat, mang_mat))

resid_list = list(residuals(ols_lm), residuals(rda_lm))

geod = dist(tgp_xy)
maxd = max(geod) / 2
geod[geod > maxd] = NA

## Mantel tests
spat_mantel = sapply(1:2, function(x) mantel(geod, dist(resid_list[[x]]),
                                             na.rm=T))

mod_names = c('OLS', 'RDA')

pdf('./figs/grid_model_mantel.pdf')
  true = !is.na(geod)
  for(i in 1:2) {
    plot(geod[true], dist(resid_list[[i]])[true], ylab='Residual Distance', xlab='Spatial Distance (m)',
         main=paste('Grid, ', mod_names[i], ', Spatial Corr', sep=''),
         type='n')
    abline(lm(dist(resid_list[[i]]) ~ geod), col='red', lwd=2)
    lines(lowess(geod[true], dist(resid_list[[i]])[true]),
          col='blue', lwd=2)
    mtext(side=3,paste('Mantel, p=', spat_mantel[4, i], sep=''))
    legend('topright', c('Linear', 'Lowess'), col=c('red', 'blue'), lty=1, 
           lwd=4, bty='n') 
  }  
dev.off()