library(vegan)
library(car)

source('./scripts/tgp_functions.R')
source('./scripts/tgp_grid_data_import.R')
source('./scripts/termplot.R')

dir.create('./results/')

## test significance of model terms
print('OLS model of richness------------------------------------')
## ols model of richness
print('non-parametric tests')
full = rda(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

soil = rda(env$sr, soil_mat, mang_mat)
permutest(soil, permutations=999)

mang = rda(env$sr, mang_mat, soil_mat)
permutest(mang, permutations=999)

print('parametric tests')
full = lm(env$sr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
            mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
Anova(full)
summary(full)
avPlots(full)

dat = data.frame(sr=env$sr, soil_mat, mang_mat)
full = lm(sr ~ ., data=dat)
full_std = lm(sr ~ . , data = data.frame(scale(dat)))
full_stats = data.frame(summary(full_std)$coefficients)
names(full_stats) = c('beta', 'se', 't', 'p')

# fix up the p-values 
full_stats$p = round(full_stats$p, 3)
full_stats$p[2] = " < 0.001"
full_stats$p[-2] = paste(" =", full_stats$p[-2])


pdf('./figs/grid_partial_regression.pdf', width=7*2, height=7*1.5)
par(mfrow=c(2,3))
ylabs = c('Soil PC1', 'Soil PC2' , 'Soil PC3', 
          'Years of Bison', '# of Burns in Past 5 Years', 'Years Since Burn')
for(i in seq_along(ylabs)) {
    termplot(full, partial.resid = T, terms=names(dat)[i + 1],
             se=T, smooth=panel.smooth, pch=19, cex=0.75,
             col.res = 'black', col.term = 'dodgerblue', lwd.term=2,
             col.se = 'dodgerblue', lwd.se = 2, lty.se=3,
             col.smth = 'red', lty.smth = 2,
             xlabs='', ylabs='', axes=F, bty='n', ylim=c(-32, 40))
    stats = substitute(paste(beta, ' = ', beta_std, ", ", italic(p), pv), 
                     list(beta_std = round(full_stats$beta[i + 1], 2),
                          pv = full_stats$p[i + 1]))
    legend('bottomleft', legend=stats, bty='n', cex=2)
    axis(side=1, cex.axis=1.75, padj=0.5)
    axis(side=2, cex.axis=1.75)
    mtext(side=1, ylabs[i], padj=2.5, cex=1.5)
    mtext(side=2, "Partial Residuals", padj=-1.75, cex=1.5)
    mtext(side=3, paste(letters[i], ")", sep=''), cex=1.5, 
          at = par("usr")[1]+0.025*diff(par("usr")[1:2]))
    if (i == 3)
        legend('topright', c('model fit', '95% CI', 'lowess line'), 
               col=c('dodgerblue','grey', 'red'), lty=c(1,1,2), lwd=c(2, 8, 2),
               cex=2, bty='n')
}
dev.off()

grid_mod = lm(scale(env$sr) ~ scale(soil_mat) + scale(mang_mat))
summary(grid_mod)
grid_coef = coef(grid_mod)[-1]
grid_err = 1.96 * summary(grid_mod)$coef[-1 , 2]


png('./figs/grid_standardized_betas.png', res=100)
plot(1:6, grid_coef, xlab='', ylab='', pch=19,
     xlim=c(0,6), ylim=c(-.6, .4), cex=1.25, frame.plot=F,
     axes=F)
axis(side=2, cex.axis=1, lwd=2)
axis(side=1, at=1:6, lab=names(dat)[-1], tick=F, cex.axis=.75)
mtext(side=2, expression('Partial Standarized '*beta), padj=-1.75, cex=1.5)
x = 1:6
y1 = grid_coef + grid_err
y2 = grid_coef - grid_err
arrows(x, y1, x, y2, angle=90, length=.1, code=3, lwd=2) 
abline(h=0, col='grey', lty=2, lwd=2)
dev.off()


print('RDA model of composition-------------------------------')
### rda model of composition --------------------
full = rda(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

soil = rda(comm_sqr ~ soil_mat$Comp.1 + soil_mat$Comp.2 + soil_mat$Comp.3 +
           Condition(mang_mat$YrsOB) + Condition(mang_mat$BP5Yrs) +
           Condition(mang_mat$YrsSLB))
permutest(soil, permutations=999)
u
mang = rda(comm_sqr, mang_mat, soil_mat)
permutest(mang, permutations=999)

sp_scr = soil$CCA$v[ , 1:2]
sp_core = c('andrgera', 'sorgnuta', 'panivirg')
indices = which(row.names(sp_scr) %in% sp_core)


rda_stat = list(list(R2 = RsquareAdj(soil)), list())

pdf('./figs/grid_rda_partial.pdf', width=7*1.5, height=7)
par(mfrow=c(1, 2))

xlim = c(-4, 3)
ylim = c(-3, 3)
plot(soil, type='n', scaling=1, axes=F, xlim=xlim, ylim=ylim,
     xlab='', ylab='')
text(soil, display='bp', labels=rep('', 3), col='red', lwd=2)
orditorp(soil, display='sp', cex=0.75, scaling=1,
         pcol='dodgerblue', pch=19)
points(sp_scr[indices, 1], sp_scr[indices, 2], pch=2, col='blue')
text(sp_scr[indices, 1], sp_scr[indices, 2], labels = sp_core, col='purple') 
x = c(2.85, 0.69, 0.55)
y = c(0.25, -2.5, 1.1)
text(x, y, labels = paste0('PC', 1:3), col='red', lwd=2)
axis(side=1, cex.axis=1.5, padj=0.5)
axis(side=2, cex.axis=1.5)
mtext(side=1, 'partial RDA axis 1', cex=1.5, padj=3.5)
mtext(side=2, 'partial RDA axis 2', cex=1.5, padj=-3)
mtext(side=3, "a)", cex=1.5, 
      at = par("usr")[1]+0.025*diff(par("usr")[1:2]))
stats = substitute(paste(italic(R^2), ' = ', r2, ", ", italic(p), ' = ', pv), 
                   list(r2 = round(lm_stat[[i]]$r.squared, 2),
                        pv = 1 - round(pf(lm_stat[[i]]$f[1], 
                                          lm_stat[[i]]$f[2],
                                          lm_stat[[i]]$f[3]), 3))) 
legend('bottomleft', legend=stats, bty='n', cex=2)


##
plot(mang, type='n', scaling=1, axes=F, xlim=xlim, ylim=ylim,
     xlab='', ylab='')
text(mang, display='bp', labels = rep('', 3), col='red', lwd=2)
orditorp(mang, display='sp', cex=0.75, scaling=1,
         pcol='dodgerblue', pch=19)
x = c(-1.6, 2.33, -3.05)
y = c(3.5, -2.45, -1.35)
text(x, y, 
     labels = c('Years\nof Bison', '# Burns Past\n5 years', 
                'Years Since Burn'),
     col='red', lwd=2)
axis(side=1, cex.axis=1.5, padj=0.5)
axis(side=2, cex.axis=1.5)
mtext(side=1, 'partial RDA axis 1', cex=1.5, padj=3.5)
mtext(side=2, 'partial RDA axis 2', cex=1.5, padj=-3)
mtext(side=3, "b)", cex=1.5, 
      at = par("usr")[1]+0.025*diff(par("usr")[1:2]))
dev.off()

print('CCA model of composition-------------------------------')
### cca model of composition
full = cca(comm_sqr ~ soil_mat[,1] + soil_mat[,2] + soil_mat[,3] +
             mang_mat$YrsOB + mang_mat$BP5Yrs + mang_mat$YrsSLB)
anova(full, by='margin')

soil = cca(comm_sqr, soil_mat, mang_mat)
permutest(soil, permutations=999)

mang = cca(comm_sqr, mang_mat, soil_mat)
permutest(mang, permutations=999)

