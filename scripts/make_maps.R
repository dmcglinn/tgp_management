
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
sort(bison@data$BEG_DATE)
nlvs = length(seq(1993, 2004, 2))
cls = colorRampPalette(c('green3','khaki1'))(nlvs)
cls = rep(cls, each =2)
yrs = strsplit(as.character(bison@data$BEG_DATE), '/')
yrs = sapply(yrs, function(x) as.numeric(x[1]))

cls_ord = cls[match(yrs, 1993:2004)]
xlims = bbox(tgpBnd)[1,]
ylims = bbox(tgpBnd)[2,]
pdf('./figs/bison_history.pdf')
  plot(bison, col=cls_ord, xlim=xlims, ylim=ylims)
  plot(tgpBnd, lwd=3, add=T)
  par(mar=c(2, 2, 2, 6))
  image(1 , 1:nlvs, matrix(1:nlvs, nrow=1),col=rev(cls), axes=F, 
        xlab='', ylab='')
  axis(side=4, at=1:nlvs, labels=seq(2003, 1993, -2),
       las=2, cex.axis=2)
dev.off()

pdf('./figs/plot_map.pdf')
  plot(tgpBnd, lwd=3)
  points(veg$easting, veg$northing, pch=19, cex=.5)
  plot(tgpBnd, lwd=3)
  points(veg$easting[veg$repeat_plot==1],
         veg$northing[veg$repeat_plot==1], pch=19, cex=1,
         col='dodgerblue')
dev.off()