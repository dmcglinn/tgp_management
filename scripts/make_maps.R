
setwd('~/Lab data/tgp_management/')

library(dichromat)

source('./scripts/tgp_functions.R')

load('./data/tgp_shpfiles.Rdata')

ls()

veg = read.csv('./data/tgp_utm_env_complete.csv')


## make burn maps -------------------------------------------
cls = rep(colorschemes$Categorical.12[seq(2, 12, 2)][-2], 3)
pdf('./figs/fire_93_03.pdf')
  for (i in 9:19) {
    plot(burns[[i]], col=cls[i-8], border=NA)
    plot(tgpBnd, lwd=3, add=T)
    points(veg$easting, veg$northing, pch=19, cex=.5)
  }
dev.off()

## make bison map --------------------------------------------
cls = colorRampPalette(c('green3','khaki1'))(9)
  
cls = cls[order(bison@data$BEG_DATE)]
xlims = bbox(tgpBnd)[1,]
ylims = bbox(tgpBnd)[2,]
pdf('./figs/bison_history.pdf')
  plot(bison, col=cls, xlim=xlims, ylim=ylims)
  plot(tgpBnd, lwd=3, add=T)
  points(veg$easting, veg$northing, pch=19, cex=.5)
dev.off()
