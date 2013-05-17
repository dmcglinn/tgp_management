library(vegan)
library(sp)
library(nlme)


setwd('~/Lab data/tgp_management/')

source('./scripts/tgp_functions.R')

load('./data/tgp_shpfiles.Rdata')

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

mang_vars = c('YrsOB', 'BP5Yrs', 'YrsSLB')
mang_mat = env[ , mang_vars]

## drop first columns
plot_mat = plot_mat[ , -1]
year_mat = year_mat[ , -1]


## composition analysis ------------------------------------------------

full = cca(comm_sqr, plot_mat, year_mat, mang_mat)

ordi_part(comm_sqr, plot_mat, year_mat, mang_mat, method='rda')
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
                     method='cca', nperm=999)


## richness analysis ----------------------------------------------------
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
