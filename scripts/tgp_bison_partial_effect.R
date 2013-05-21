
setwd('~/tgp_management/')

## estimate scaled parital effect size

source('./scripts/tgp_grid_data_import.R')

grid_mod = lm(scale(env$sr) ~ scale(soil_mat) + scale(mang_mat))
grid_coef = coef(grid_mod)[5]
grid_err = 1.96 * summary(grid_mod)$coef[5 , 2]

source('./scripts/tgp_repeat_data_import.R')
rep_mod = lm(scale(env$sr) ~ scale(soil_mat) + scale(mang_mat) + scale(rain_mat))
rep_coef = coef(rep_mod)[5]
rep_err = 1.96 * summary(rep_mod)$coef[5 , 2]

pdf('./figs/Bison_effect_comparison.pdf')
plot(1:2, c(grid_coef, rep_coef), xlab='', ylab='', pch=19, 
     xlim= c(0, 3), ylim=c(-.1, .6), cex = 1.25,
     frame.plot=F, axes=F)
axis(side=2, cex.axis=1, lwd=2, at = seq(-.1, .6, .1))
axis(side=1, at=c(1, 2), lab=c('Grid','Repeat'), tick=F, 
     cex.axis=1.5)
x = 1:2
y1 = c(grid_coef + grid_err, rep_coef + rep_err)
y2 = c(grid_coef - grid_err, rep_coef - rep_err)
arrows(x, y1, x, y2, angle=90, length=.1, code=3, lwd=2) 
abline(h=0, col='grey', lty=2, lwd=2)
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
