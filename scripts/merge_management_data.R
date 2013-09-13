setwd('~/tgp_management/')

env = read.csv('./data/tgp_utm_env.csv')
burn = read.csv('./data/plot_burn.csv')
graze = read.csv('./data/plot_graze.csv')

class(env$date_samp)

env$date_samp = as.Date(env$date_samp)
burn$burn_date = as.Date(burn$burn_date)
graze$bison_date= as.Date(graze$bison_date)

yrs_since_burn = rep(NA, nrow(env))
bp5yrs = rep(NA, nrow(env))
yrs_of_bison = rep(NA, nrow(env))

for (i in seq_along(env$plot)) {
  burn_date_tmp = burn$burn_date[burn$plot == env$plot[i]]
  bison_date_tmp = graze$bison_date[graze$plot == env$plot[i]]
  burn_daydiff = env$date_samp[i] - burn_date_tmp
  bison_daydiff = env$date_samp[i] - bison_date_tmp
  burn_daydiffpos = burn_daydiff[burn_daydiff > 0]
  bison_daydiffpos = bison_daydiff[bison_daydiff > 0]
  if (length(burn_daydiffpos) > 0) {
    days_since_burn = as.numeric(min(burn_daydiffpos))
    yrs_since_burn[i] = days_since_burn / 365  
    bp5yrs[i] = sum(burn_daydiffpos/365 <= 5)
  }  
  if (length(bison_daydiffpos) > 0) {
    days_of_bison = as.numeric(min(bison_daydiffpos))
    yrs_of_bison[i] = days_of_bison / 365  
  }  
}

yrs_of_bison = ifelse(is.na(yrs_of_bison), 0, yrs_of_bison)
bison = (yrs_of_bison >= 1) * 1
ungrazed = graze$ungrazed[match(env$plot, graze$plot)]

## four plots have yrs_since_burn of NA
env[is.na(yrs_since_burn), ]
## these plots were all sampled early in the devlopment of the 
## preserve
## we will simply assign these plots the mean yrs_since_burn 
## value
yrs_since_burn[is.na(yrs_since_burn)] = mean(yrs_since_burn, na.rm=T)
bp5yrs[is.na(bp5yrs)] = round(mean(bp5yrs, na.rm=T))
burn = (yrs_since_burn < 1) * 1

## look for differences with old variables
par(mfrow=c(1,3))
plot(env$YrsOB - yrs_of_bison, ylim=c(-.1,.1))
abline(h=0)
plot(env$YrsSLB - yrs_since_burn)
abline(h=0)
plot(env$BP5Yrs - bp5yrs)
abline(h=0)

## replace old data with newly derived variables
env$bison = bison
env$YrsOB = yrs_of_bison
env$BP5Yrs = bp5yrs
env$YrsSLB = yrs_since_burn
env$burn = burn

## export newly merged envio data
write.csv(env, file='./data/tgp_utm_env_complete.csv', row.names=F)

## error checking -------------------------------------------------------------------
load('./data/tgp_shpfiles.Rdata')
library(sp)

head(burns$'1997'@data)
spplot(burns$'1997', 'DATE')

## visually check plot 206
plot(pasture, axes=T, xlim=c(7.25e5, 7.35e5), ylim=c(4.0755e6, 4.0765e6))
#points(env$easting, env$northing, pch=19)
par(new=TRUE)
plot(burns$'2007', col='green', lty=2, xlim=c(7.25e5, 7.35e5), ylim=c(4.0755e6, 4.0765e6))
points(env$easting[2], env$northing[2], pch=19, col='red')

