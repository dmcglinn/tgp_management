
setwd('~/tgp_management/')

## estimate scaled parital effect size

source('./scripts/tgp_grid_data_import.R')

grid_mod = lm(scale(env$sr) ~ scale(mang_mat) + scale(soil_mat))
grid_coef = coef(grid_mod)[2:4]
grid_err = 1.96 * summary(grid_mod)$coef[2:4 , 2]

source('./scripts/tgp_repeat_data_import.R')
rep_mod = lm(scale(env$sr) ~ scale(mang_mat) + scale(plot_mat) + scale(year_mat))
rep_coef = coef(rep_mod)[2:4]
rep_err = 1.96 * summary(rep_mod)$coef[2:4 , 2]

pdf('./figs/bison_effect_comparison.pdf')
plot(1:2, c(grid_coef[1], rep_coef[1]), xlab='', ylab='',
     pch=19, xlim= c(0, 3), ylim=c(-.1, .6), cex = 1.25,
     frame.plot=F, axes=F)
axis(side=2, cex.axis=1, lwd=2, at = seq(-.1, .6, .1))
axis(side=1, at=c(1, 2), lab=c('Grid','Repeat'), tick=F, 
     cex.axis=1.5)
mtext(side=2, 'Partial Standarized Effect Size', padj=-3, cex=1.5)
x = 1:2
y1 = c(grid_coef[1] + grid_err[1], rep_coef[1] + rep_err[1])
y2 = c(grid_coef[1] - grid_err[1], rep_coef[1] - rep_err[1])
arrows(x, y1, x, y2, angle=90, length=.1, code=3, lwd=2) 
abline(h=0, col='grey', lty=2, lwd=2)
dev.off()


pdf('./figs/manag_effect_comparison.pdf')
incr = .2/2
plot(1:3 - incr, grid_coef, xlab='', ylab='',
     pch=19, xlim= c(.5, 3.5), ylim=c(-.4, .8), cex = 1.25,
     frame.plot=F, axes=F, col='grey')
points(1:3 + incr, rep_coef, pch=19, cex=1.25)
axis(side=2, cex.axis=1, at=seq(-0.4, 0.8, .2), lwd=2)
axis(side=1, at=1:3,
     lab=c('Years of\n bison', '# of burns\n past 5 years', 'Years since\n last burn'),
     tick=F,   cex.axis=1.25)
mtext(side=2, 
      expression('Partial Standarized Coefficient, '* hat(beta)),
      padj=-1.5, cex=1.25)
x = 1:3
y1 = grid_coef + grid_err
y1b = rep_coef + rep_err
y2 = grid_coef - grid_err
y2b = rep_coef - rep_err
arrows(x-incr, y1, x-incr, y2, angle=90, length=.1, code=3, lwd=2, col='grey') 
arrows(x+incr, y1b, x+incr, y2b,angle=90, length=.1, code=3, lwd=2)
abline(h=0, col='grey', lty=2, lwd=2)
legend('topright',c('Grid Analysis', 'Repeat Analysis'), col=c('grey','black'),
       bty='n', pch=19, cex=1.25)
dev.off()


##
barplot(c(grid_coef, rep_coef), width=.25, space=.6,
        xlim=c(0, 1), ylim=c(-.1, .7), names.arg=rep('',2) )
axis(side=1)
x = c(.275, .675)
y1 = c(grid_coef + grid_err, rep_coef + rep_err)
y2 = c(grid_coef - grid_err, rep_coef - rep_err)
arrows(x, y1, x, y2, angle=90, length=.1, code=3)  

## temporal extent analysis
source('./scripts/tgp_repeat_data_import.R')
yrs = unique(env$yr)

output = matrix(NA, ncol=10, nrow=length(yrs) - 1)
irow = 1
for (i in yrs[-1]) {
  true = env$yr <= i
  mod = lm(scale(env$sr[true]) ~ scale(soil_mat[true,]) + scale(mang_mat[true,]) +
             scale(rain_mat[true,]))
  rep_coef = coef(mod)[5:7]
  rep_err = 1.96 * summary(mod)$coef[5:7 , 2]
  output[irow , 1] = i
  tmp = as.vector(t(matrix(c(rep_coef - rep_err, rep_coef, rep_coef + rep_err),
                           ncol=3)))
  output[irow , 2:10] = tmp
  irow=irow+1
}

par(mfrow=c(1,3))
for(i in 1:3) {
  index = 3 * (i - 1)
  plot(output[, 1 ], output[, index + 3], ylim=c(-1, 1), type='o',
       ylab=names(mang_mat)[i], xlab='year')
  arrows(output[,1], output[, index + 2],
         output[,1], output[, index + 4],
         code=3, angle=90, length=.1)
  points(1999.5, grid_coef, pch=19)
  arrows(1999.5, grid_coef - grid_err, 
         1999.5, grid_coef + grid_err,
         code=3, angle=90, length=.1)
}

## with moving window
output = matrix(NA, ncol=10, nrow=sum(11:1))
irow = 1
lyrs = length(yrs)
for (w in 2:lyrs) {
  for (i in 1:(lyrs - w + 1)) {
    indices = i:(i + w - 1)
    true = env$yr %in% yrs[indices]
    mod = lm(scale(env$sr[true]) ~ scale(soil_mat[true,]) + scale(mang_mat[true,]) +
             scale(rain_mat[true,]))
    rep_coef = coef(mod)[5:7]
    rep_err = 1.96 * summary(mod)$coef[5:7 , 2]
    output[irow , 1] = w
    tmp = as.vector(t(matrix(c(rep_coef - rep_err, rep_coef, rep_coef + rep_err),
                             ncol=3)))
    output[irow , 2:10] = tmp
    irow=irow+1
  }
}

pdf('./figs/repeat_mang_effects_time_scale.pdf', width=7 * 2, height=7*1.5)
par(mfrow=c(2,3))
for(i in 1:3) {
  index = 3 * (i - 1)
  plot(output[, 1 ], output[, index + 3], ylim=c(-1, 1), type='p',
       ylab=names(mang_mat)[i], xlab='time span (yr)')
  lines(lowess(output[, 1], output[, index + 3]), col='red', lwd=2)
}


output = aggregate(output, by = list(output[,1]), mean)[ , -1]

for(i in 1:3) {
  index = 3 * (i - 1)
  plot(output[, 1 ], output[, index + 3], ylim=c(-1, 1), type='p',
       ylab=names(mang_mat)[i], xlab='time span (yr)')
  arrows(output[,1], output[, index + 2],
         output[,1], output[, index + 4],
         code=3, angle=90, length=.1)
}
dev.off()
