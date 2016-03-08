library(vegan)

source('./scripts/tgp_functions.R')
source('./scripts/tgp_repeat_data_import.R')
source('./scripts/termplot.R')

dir.create('./results/')

## test significance of model terms
print('OLS model of richness-----------------------------------')
## ols model of richness
print('non-parametric tests')
full = rda(env$sr ~ plot_mat + year_mat + 
                    mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

site = rda(env$sr, plot_mat, cbind(year_mat, mang_mat))
permutest(site, permutations=999)

year = rda(env$sr, year_mat, cbind(plot_mat, mang_mat))
permutest(year, permutations=999)

mang = rda(env$sr, mang_mat, cbind(year_mat, plot_mat))
permutest(mang, permutations=999)

print('parametric tests')
full = lm(env$sr ~ as.factor(env$plot) + as.factor(env$yr) + 
          mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
Anova(full)
summary(full)
avPlots(full)

dat= data.frame(sr=env$sr, plot_mat, year_mat, mang_mat)
full = lm(sr ~ ., data=dat)

png('./figs/repeat_partial_regression.png', width=480*3, height=480*1,
    res=100)
par(mfrow=c(1,3))
for(i in c("YrsOB","BP5Yrs","YrsSLB")) {
    termplot(full, partial.resid = T, terms=i,
             se=T, smooth=panel.smooth,
             col.term = 'blue', lwd.term=2,
             col.se = 'blue', lwd.se = 2, lty.se=3,
             col.smth = 'red', lty.smth = 2,
             xlab='', ylab='', axes=F, bty='n')
    axis(side=1, cex.axis=1.75, padj=0.5)
    axis(side=2, cex.axis=1.75)
    mtext(side=1, i, padj=2.5, cex=1.5)
    mtext(side=2, "Partial Residuals", padj=-1.75, cex=1.5)
    if (i == 'YrsSLB')
        legend('topright', c('model fit', '95% CI', 'lowess line'), 
               col=c('blue','blue', 'red'), lty=c(1,3,2), lwd=2,
               cex=2, bty='n')
}
dev.off()

grid_mod = lm(scale(env$sr) ~ scale(plot_mat) + scale(year_mat) + scale(mang_mat))
grid_coef = coef(grid_mod)[32:34]
grid_err = 1.96 * summary(grid_mod)$coef[32:34 , 2]

png('./figs/repeat_standardized_betas.png', res=100)
plot(1:3, grid_coef, xlab='', ylab='', pch=19,
     xlim=c(0,3), ylim=c(-.4, .8), cex=1.25, frame.plot=F,
     axes=F)
axis(side=2, cex.axis=1, lwd=2)
axis(side=1, at=1:3, lab=names(dat)[32:34], tick=F, cex.axis=1)
mtext(side=2, expression('Partial Standarized '*beta), padj=-1.75, cex=1.5)
x = 1:3
y1 = grid_coef + grid_err
y2 = grid_coef - grid_err
arrows(x, y1, x, y2, angle=90, length=.1, code=3, lwd=2) 
abline(h=0, col='grey', lty=2, lwd=2)
dev.off()

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



