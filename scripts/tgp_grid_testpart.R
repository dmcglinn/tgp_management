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

