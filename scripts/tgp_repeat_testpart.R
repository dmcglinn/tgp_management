library(vegan)
library(car)

source('./scripts/tgp_functions.R')
source('./scripts/tgp_repeat_data_import.R')
source('./scripts/termplot.R')

env$plot = as.factor(env$plot)
env$yr = as.factor(env$yr)

dir.create('./results/')

## test significance of model terms
print('OLS model of richness-----------------------------------')
## ols model of richness-------------
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
plotID = factor(env$plot)
year = as.factor(env$yr)
full = lm(env$sr ~ plotID + year + 
          mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)

plotID_res = residuals(full, 'partial')[ , 1]
year_res = residuals(full, 'partial')[ , 2]
ylim = range(plotID_res, year_res)
#ID = factor(ID, levels(ID)[order(tapply(ID_res, ID, mean))])
pdf('./figs/repeat_partial_regression_factors.pdf', width=7*1.33)
xlabs = c('Plot ID', 'Year')
par(mfrow=c(1,2))
#for(i in seq_along(xlabs)) {
    boxplot(plotID_res ~ plotID, frame.plot=F, ylim=ylim)
    boxplot(year_res ~ year, frame.plot=F, ylim=ylim)
#    axis(side=1, labels= cex.axis=1.75, padj=0.5)
#    axis(side=2, cex.axis=1.75)
#    mtext(side=1, xlabs[i], padj=2.5, cex=1.5)
#    mtext(side=2, "Partial Residuals", padj=-1.75, cex=1.5)
#}
dev.off()

Anova(full)
summary(full)
avPlots(full)

dat = data.frame(sr=env$sr, plot_mat, year_mat, mang_mat)
full = lm(sr ~ ., data=dat)
full_std = lm(sr ~ . , data = data.frame(scale(dat)))
full_stats = data.frame(summary(full_std)$coefficients)
names(full_stats) = c('beta', 'se', 't', 'p')

# fix up the p-values 
full_stats$p[32] = " < 0.001"
full_stats$p[33] = " = 0.01"
full_stats$p[34] = " = 0.001"


pdf('./figs/repeat_partial_regression.pdf', width = 7*2, height = 7*1)
par(mfrow = c(1, 3))
xlabs = c('Years of Bison', '# of Burns in Past 5 Years', 'Years Since Burn')
for(i in seq_along(xlabs)) {
    termplot(full, partial.resid = T, terms=names(dat)[i + 31],
             se=T, smooth=panel.smooth, col.res = 'black',
             col.term = 'dodgerblue', lwd.term=2,
             col.se = 'dodgerblue', lwd.se = 2, lty.se=3,
             col.smth = 'red', lty.smth = 2,
             xlabs='', ylabs='', axes=F, bty='n', ylim = c(-25, 25))
    stats = substitute(paste(beta, ' = ', beta_std, ", ", italic(p), pv), 
                  list(beta_std = round(full_stats$beta[i + 31], 2),
                       pv = full_stats$p[i + 31]))
    legend('bottomleft', legend=stats, bty='n', cex=2)
    axis(side=1, cex.axis=1.75, padj=0.5)
    axis(side=2, cex.axis=1.75)
    mtext(side=1, xlabs[i], padj=2.5, cex=1.5)
    mtext(side=2, "Partial Residuals", padj=-1.75, cex=1.5)
    mtext(side=3, paste(letters[i], ")", sep=''), cex=1.5, 
          at = par("usr")[1] + 0.025*diff(par("usr")[1:2]))
    if (i == 3)
        legend('topright', c('model fit', '95% CI', 'lowess line'), 
               col=c('dodgerblue','grey', 'red'), lty=c(1,1,2), lwd=c(2, 8, 2),
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
### rda model of composition------------
full = rda(comm_sqr ~ plot_mat + year_mat +
                      mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

site = rda(comm_sqr ~ plot + Condition(yr) + Condition(YrsOB) + Condition(BP5Yrs) +
           Condition(YrsSLB), data=env)
permutest(site, permutations=999)

year = rda(comm_sqr ~ Condition(plot) + yr + Condition(YrsOB) + Condition(BP5Yrs) +
           Condition(YrsSLB), data=env)
permutest(year, permutations=999)

mang = rda(comm_sqr ~ Condition(plot) + Condition(yr) + YrsOB + BP5Yrs + YrsSLB,
           data=env)
permutest(mang, permutations=999)

### plot RDA partial plots--------

pdf('./figs/repeat_rda.pdf', width=7*1.5)
par(mfrow=c(1,2))
xlims = c(-3.75, 5.75)
ylims = c(-6, 4)
results = list(site, year, mang)
yr_labels = substr(as.character(1998:2009), 3,4)
yr_cols = colorRampPalette(c('black', 'red'))
for (i in seq_along(results)) {
    plot(results[[i]], type='n', scaling=1, axes=F, xlab='', ylab='', 
         xlim=xlims, ylim=ylims)
    if (i == 1) {
        orditorp(results[[i]], 'sp', cex=0.75, scaling=1, col=1,
               pcol='dodgerblue', pch=19)  
        points(results[[i]], 'cn', pch=2, col='red', cex=1.5)
    } 
    if (i == 2) { 
        orditorp(results[[i]], 'sp', cex=0.75, scaling=1, col=1,
               pcol='dodgerblue', pch=19)        
        points(year, 'cn', pch=2, col='red', cex=1.5)
    } 
    if (i == 3) {
        orditorp(results[[i]], 'sp', cex=0.75, scaling=1, col=1,
                 pcol='dodgerblue', pch=19) 
        text(mang, 'bp', lwd=2, 
             labels = c('Years of Bison', '# Burns Past 5 years',
                        'Years Since Burn'),
             col='red', cex=1, axis.bp = FALSE)
 
    }
    axis(side=1, cex.axis=1.5, padj=0.5)
    axis(side=2, cex.axis=1.5)
    mtext(side=1, 'partial RDA axis 1', cex=1.5, padj=3.5)
    mtext(side=2, 'partial RDA axis 2', cex=1.5, padj=-3)
    mtext(side=3, paste(letters[i], ")", sep=''), cex=1.5, 
          at = par("usr")[1]+0.025*diff(par("usr")[1:2]))
}
dev.off()


print('CCA model of composition-------------------------------')
### cca model of composition------------
full = cca(comm_sqr ~ plot_mat + year_mat +
                      mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

site = cca(comm_sqr ~ plot + Condition(yr) + Condition(YrsOB) + Condition(BP5Yrs) +
           Condition(YrsSLB), data=env)
permutest(site, permutations=999)

year = cca(comm_sqr ~ Condition(plot) + yr + Condition(YrsOB) + Condition(BP5Yrs) +
           Condition(YrsSLB), data=env)
permutest(year, permutations=999)

mang = cca(comm_sqr ~ Condition(plot) + Condition(yr) + YrsOB + BP5Yrs + YrsSLB,
           data=env)
permutest(mang, permutations=999)

pdf('./figs/repeat_cca.pdf', width=7*1.5)
par(mfrow=c(1,2))
xlims = c(-15, 15)
ylims = c(-15, 15)
results = list(site, year, mang)
yr_labels = substr(as.character(1998:2009), 3,4)
yr_cols = colorRampPalette(c('black', 'red'))
for (i in seq_along(results)) {
    plot(results[[i]], type='n', scaling=1, axes=F, xlab='', ylab='', 
         xlim=xlims, ylim=ylims)
    orditorp(results[[i]], 'sp', cex=0.75, scaling=1, col=1,
             pcol='dodgerblue', pch=19)
    if (i == 1) 
        points(results[[i]], 'cn', pch=2, col='red', cex=1.5)
    if (i == 2)
        points(year, 'cn', pch=2, col='red', cex=1.5)
    if (i == 3)
        text(mang, 'bp',
             labels = c('Years of Bison', '# Burns Past 5 years',
                        'Years Since Burn'),
             col='red', cex=1, axis.bp = FALSE)
    axis(side=1, cex.axis=1.5, padj=0.5)
    axis(side=2, cex.axis=1.5)
    mtext(side=1, 'partial CCA axis 1', cex=1.5, padj=3.5)
    mtext(side=2, 'partial CCA axis 2', cex=1.5, padj=-3)
}
dev.off()


