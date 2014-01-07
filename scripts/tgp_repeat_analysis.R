library(vegan)
library(sp)
library(nlme)

source('./scripts/tgp_functions.R')

env = read.csv('./data/tgp_utm_env_complete.csv')
comm = read.csv('./data/tgp_comm_mat_all.csv')

env = env[env$repeat_plot == 1, ]

comm = comm[match(env$plot_yr, comm$plot.yr), ]
## drop species that don't occur
comm = comm[ , apply(comm, 2, sum) > 0 ] 

nrow(env)
nrow(comm)

all.equal(env$plot_yr, comm$plot.yr)

row.names(comm) = comm$plot.yr
comm = comm[ , -1]
comm_sqr = sqrt(comm)

env$sr = rowSums(comm > 0)

plot_id = sort(unique(env$plot))
year_id = sort(unique(env$yr))
plot_mat = matrix(0, ncol=length(plot_id), nrow=nrow(env))
year_mat = matrix(0, ncol=length(year_id), nrow=nrow(env))

for(i in 1:nrow(env)) {
  plot_mat[i, match(env$plot[i], plot_id)] = 1
  year_mat[i, match(env$yr[i], year_id)] = 1 
}

## drop first columns so no singular variables in models
plot_mat = plot_mat[ , -1]
year_mat = year_mat[ , -1]

## define management variables
mang_vars = c('YrsOB', 'BP5Yrs', 'YrsSLB')
mang_mat = env[ , mang_vars]

## define temporal precip variables
## summer : june-sept : 6, 7, 8 ,9
## winter : oct-jan : 10, 11, 12, 1
## spring : feb-may : 2, 3, 4, 5
env$rain2 = ifelse(is.na(env$rain2), 0, env$rain2)
sum_rain = apply(env[ , paste('rain', 6:9, sep='')], 1, sum)
win_rain = apply(env[ , paste('rain', c(10:12, 1), sep='')], 1, sum)
spr_rain = apply(env[ , paste('rain', 2:5, sep='')], 1, sum)

rain_mat = cbind(sum_rain, win_rain, spr_rain)

## alternatively define using pca - does not change result much
rain_vars = paste('rain', 1:12, sep='')
rain_pca = princomp(scale(env[ , rain_vars]))
par(mfrow=c(1,2))
plot(rain_pca)
biplot(rain_pca)

rain_mat = as.data.frame(rain_pca$scores[ , 1:3])

## define spatial soil variables
soil_vars = c("P","CA","MG","K","NA","B","FE","MN","CU","ZN","AL")
soil_vars = paste('log', soil_vars, sep='')
soil_pca = princomp(scale(env[ , soil_vars]))
summary(soil_pca)
sum(soil_pca$sdev)
par(mfrow=c(1,2))
plot(soil_pca)
biplot(soil_pca)
par(mfrow=c(1,1))
biplot(soil_pca, col=c('white', 'black'))
biplot(soil_pca, choices = c(3,4), col=c('white', 'black'))

soil_mat = as.data.frame(soil_pca$scores[ , 1:3])

## define site variables
site_mat = env[ , c('logCA', 'slope', 'northness')]

## composition analysis ------------------------------------------------

ordi_part(comm_sqr, site_mat, rain_mat, mang_mat, method='rda')$part
                     R2  R2adj  % expl
all = X1+X2+X3    0.224  0.193 100.000
[a] = X1 | X2+X3  0.124  0.117  60.480
[b] = X2 | X1+X3  0.020  0.010   5.150
[c] = X3 | X1+X2  0.066  0.057  29.524
[d]              -0.001 -0.002  -1.151
[e]               0.001  0.000  -0.039
[f]               0.014  0.012   6.327
[g]               0.000 -0.001  -0.291
[h] = Residuals   0.776  0.807 417.800

ordi_part(comm_sqr, soil_mat, rain_mat, mang_mat, method='rda')$part
                     R2  R2adj  % expl
all = X1+X2+X3    0.247  0.217 100.000
[a] = X1 | X2+X3  0.147  0.141  64.886
[b] = X2 | X1+X3  0.022  0.013   5.935
[c] = X3 | X1+X2  0.063  0.054  24.888
[d]              -0.003 -0.005  -2.383
[e]               0.001  0.000   0.079
[f]               0.017  0.015   6.967
[g]              -0.001 -0.001  -0.372
[h] = Residuals   0.753  0.783 360.075

rda_part = ordi_part(comm_sqr, plot_mat, year_mat, mang_mat, method='rda')
$part
                     R2  R2adj  % expl
all = X1+X2+X3    0.610  0.548 100.000
[a] = X1 | X2+X3  0.460  0.450  82.206
[b] = X2 | X1+X3  0.060  0.043   7.803
[c] = X3 | X1+X2  0.013  0.008   1.462
[d]               0.010 -0.014  -2.529
[e]               0.012  0.012   2.178
[f]               0.066  0.063  11.455
[g]              -0.010 -0.014  -2.575
[h] = Residuals   0.390  0.452  82.507

cca_part = ordi_part(comm_sqr, plot_mat, year_mat, mang_mat,
                     method='cca', nperm=100)
$part
R2  R2adj  % expl
all = X1+X2+X3    0.511  0.432 100.000
[a] = X1 | X2+X3  0.405  0.382  88.351
[b] = X2 | X1+X3  0.043  0.018   4.254
[c] = X3 | X1+X2  0.014  0.008   1.930
[d]               0.005 -0.014  -3.237
[e]               0.001  0.001   0.242
[f]               0.047  0.043  10.012
[g]              -0.005 -0.007  -1.552
[h] = Residuals   0.489  0.568 131.318

## richness analysis ----------------------------------------------------

ols_part = ordi_part(env$sr, plot_mat, year_mat, mang_mat, method='rda')
$part
R2  R2adj  % expl
all = X1+X2+X3    0.755  0.716 100.000
[a] = X1 | X2+X3  0.481  0.487  68.027
[b] = X2 | X1+X3  0.127  0.126  17.594
[c] = X3 | X1+X2  0.038  0.040   5.548
[d]               0.002 -0.031  -4.294
[e]               0.036  0.035   4.915
[f]               0.073  0.066   9.262
[g]              -0.002 -0.008  -1.053
[h] = Residuals   0.245  0.284  39.719

mod = lm(sr ~ plot_mat + year_mat + 
           YrsOB + YrsSLB + BP5Yrs, data=env)
summary(mod)
new = data.frame(YrsOB = seq(0, 15, 1))
tst = predict(mod, interval='confidence',
              type='terms', terms='YrsOB')
tst = predict(mod, new, interval='confidence')


termplot(mod, terms='YrsOB', partial=T, se=T,
         lwd.term=3,lwd.se=3,
         col.term='blue', col.se='blue',
         col.res = 'grey', col.smth = "red",
         frame.plot=F, axes=F, xlab='', ylab='',
         ylim=c(-30, 30))
axis(side=1, cex.axis=2)
axis(side=2, cex.axis=2)
res = residuals(mod, 'partial')
res = res[ , 'YrsOB', drop=FALSE]
lines(lowess(env$YrsOB, res), col='red', lty=2, lwd=3)
