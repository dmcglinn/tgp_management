library(vegan)

source('./scripts/tgp_functions.R')
source('./scripts/tgp_repeat_data_import.R')

dir.create('./results/')

## test significance of model terms
print('OLS model of richness-----------------------------------')
## ols model of richness
full = rda(env$sr ~ plot_mat + year_mat + 
                    mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

summary(lm(env$sr ~ plot_mat + year_mat + 
                    mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB))

anova(full, by='margin')

site = rda(env$sr, plot_mat, cbind(year_mat, mang_mat))
permutest(site, permutations=999)

year = rda(env$sr, year_mat, cbind(plot_mat, mang_mat))
permutest(year, permutations=999)

mang = rda(env$sr, mang_mat, cbind(year_mat, plot_mat))
permutest(mang, permutations=999)

print('RDA model of composition-------------------------------')
### rda model of composition
full = rda(comm_sqr ~ plot_mat + year_mat +
                      mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

site = rda(comm_sqr, plot_mat, cbind(year_mat, mang_mat))
permutest(site, permutations=999)

year = rda(comm_sqr, year_mat, cbind(plot_mat, mang_mat))
permutest(year, permutations=999)

mang = rda(comm_sqr, mang_mat, cbind(year_mat, plot_mat))
permutest(mang, permutations=999)


print('CCA model of composition-------------------------------')
### cca model of composition
full = cca(comm_sqr ~ plot_mat + year_mat +
                      mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

site = cca(comm_sqr, plot_mat, cbind(year_mat, mang_mat))
permutest(site, permutations=999)

year = cca(comm_sqr, year_mat, cbind(plot_mat, mang_mat))
permutest(year, permutations=999)

mang = cca(comm_sqr, mang_mat, cbind(year_mat, plot_mat))
permutest(mang, permutations=999)



