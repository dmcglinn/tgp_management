library(vegan)
library(dichromat)
library(sp)

source('./scripts/tgp_functions.R')

load('./data/tgp_shpfiles.Rdata')

source('./scripts/tgp_grid_data_import.R')

tgp_xy = cbind(env$easting, env$northing)
pdf('./figs/grid_analysis_variables.pdf', width=7*3, height=7*2)
  par(mfrow=c(2,3))
  ordisurf(tgp_xy, env$BP5Yrs, bubble=5)
  ordisurf(tgp_xy, env$YrsSLB, bubble=5)
  ordisurf(tgp_xy, env$bison, bubble=5)
  ordisurf(tgp_xy, soil_mat[,1], bubble=5)
  ordisurf(tgp_xy, soil_mat[,2], bubble=5)
  ordisurf(tgp_xy, soil_mat[,3], bubble=5)
dev.off()

env = read.csv('./data/tgp_utm_env_complete.csv')
env[is.na(env$waterpct), ]$waterpct = 0
env[is.na(env$rockpct), ]$rockpct = 0

grassland = env$waterpct == 0 &
            env$rockpct <= 20 &
            env$woodypct <= 20

## make burn maps -------------------------------------------
cls = rep(colorschemes$Categorical.12[seq(2, 12, 2)][-2], 3)
pdf('./figs/fire_93_03.pdf')
  for (i in 9:19) {
    plot(burns[[i]], col=cls[i-8], border=NA)
    plot(tgpBnd, lwd=3, add=T)
    points(env$easting, env$northing, pch=19, cex=.5)
  }
dev.off()

## make bison map --------------------------------------------
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
  ## legend (up and down)
  par(mar=c(2, 2, 2, 6))
  image(1 , 1:nlvs, matrix(1:nlvs, nrow=1),col=rev(cls), axes=F, 
        xlab='', ylab='')
  axis(side=4, at=1:nlvs, labels=seq(2003, 1993, -2),
       las=2, cex.axis=2, tick=FALSE)
  ## legend (side to side)
  par(mar=c(2, 2, 2, 2))
  image(1:nlvs , 1, matrix(1:nlvs, ncol=1),col=cls, axes=F, 
        xlab='', ylab='')
  axis(side=1, at=1:nlvs, labels=seq(1993, 2003, 2),
       las=1, cex.axis=2, tick=FALSE)
dev.off()

pdf('./figs/plot_map.pdf')
  plot(tgpBnd, lwd=3)
  points(env$easting, env$northing, pch=19, cex=.5)
  plot(tgpBnd, lwd=3)
  points(env$easting[env$repeat_plot==1],
         env$northing[env$repeat_plot==1], pch=19, cex=1,
         col='dodgerblue')
dev.off()

pdf('./figs/fig1_maps_for_ms.pdf', width=7*3, height=7)
#  x = 724500
#  y = 4086000
  par(mfrow=c(1,3))
  ## management map
  plot(bison, col=cls_ord, xlim=xlims, ylim=ylims)
  plot(tgpBnd, lwd=3, add=T)
#  text(x, y, '(a)', cex=5)
  ## grid plot map
  plot(tgpBnd, lwd=3)
  points(env$easting[grassland], env$northing[grassland],
         pch=19, cex=1.5)
#  text(x, y, '(b)', cex=5)
  ## repeat plot map
  plot(tgpBnd, lwd=3)
  points(env$easting[env$repeat_plot==1],
         env$northing[env$repeat_plot==1], pch=19, cex=1.5)
#  text(x, y, '(c)', cex=5)
dev.off()
