
## Purpose
## to create a burn matrix 
## rows = one fire event in a specific plot
## columns = plot, date, 

library(sp)
library(dichromat)


setwd('~/Lab data/tgp_management/')

load('./data/tgp_shpfiles.Rdata')

plots = read.csv('./data/tgp_plot_utm.csv')
plots_sp = SpatialPointsDataFrame(plots[ , c('easting', 'northing')],
                                  data = plots)
proj4string(plots_sp) = CRS(proj4string(pasture))


plot(pasture)
points(plots$easting, plots$northing, pch=19, col='red')


## extract burn data ----------------------------------------------------------------
## loop through each plot and each burn layer
## evaluate if the plot fell within any of the burn polygons
## record the date when it did

## here is one quick and potentially incorrect method b/c it does
## not consider the possibility that the burn polygons 
## in a single year overlap each other (which is not always 
## true in this case)
tst = sapply(burns, function(x) over(plots_sp, x))

## for a slower and more through method we need this loop
date_column = sapply(burns, function(x) grep('DATE', names(x)))
flag = 0
for (i in seq_along(plots$plot)) { 
  for (j in seq_along(burns)) {
    for (k in seq_along(burns[[j]])) {
      poly_coords = burns[[j]]@polygons[[k]]@Polygons[[1]]@coords
      pt_in = point.in.polygon(plots$easting[i], plots$northing[i], 
                               poly_coords[ , 1], poly_coords[ , 2])
      if (pt_in == 1) {
        burn_date = burns[[j]]@data[k, date_column[j]]
        if (flag == 0) {
          plot_burn = data.frame(plot=plots$plot[i], burn_date)
          flag = 1 
        }  
        else {
          plot_burn = rbind(plot_burn, 
                            data.frame(plot=plots$plot[i], burn_date))
        }  
      }
    }
  }  
}  

## this sort suggests these results are plausible
sort(table(plot_burn[,1]), dec=T)

## convert dates into a consistent format
bd = as.character(plot_burn$burn_date)
bd = gsub('/', '', bd)
bd = as.Date(bd, '%Y%m%d')
plot_burn$burn_date = bd

sum(sapply(tst, function(x) sum(!is.na(x[,1]))))
nrow(plot_burn)
## so it appears that there were 3 burns that 
## overlapped in a single year period thus justifing the 
## more complex and inefficient looping procedure

## export the burn data product
write.csv(plot_burn, file = './data/plot_burn.csv', row.names=F)

## extract grazing data--------------------------------------------------------------
spplot(pasture, 'BISON_2000')
plot(pasture)
points(coordinates(plots_sp), pch=19, col='red')
spplot(pasture, 'USAGE', col.regions=colorschemes$Categorical.12[1:9])
## so some of the plots occur in the ungrazed area

spplot(bison, 'BEG_DATE', col.regions=colorschemes$Categorical.12[1:9])
spplot(bison, 'UNIT')

## first extract the start date of bison grazing for each plot
## Note
## it is unclear how to treat the bison trap area, the bison
## layer indicates that grazing begain and then was stopped by bison
## I'm not sure how this area was used 

plot_bison = over(plots_sp, bison)
bison_date = plot_bison$BEG_DATE

plot_pasture = over(plots_sp, pasture)
table(plot_pasture$USAGE)
plots[which(plot_pasture$USAGE == 'ungrazed'), ]

plot_graze = data.frame(plot = plots$plot, 
                        bison_date = as.Date(bison_date),
                        ungrazed = plot_pasture$USAGE == 'ungrazed')

## export the grazing data product
write.csv(plot_graze, file = './data/plot_graze.csv', row.names=F)


## Overlay examples -----------------------------------------------------------------
## the example below demonstrates that the functions
## overlay and the newer function over only return the 
## first polygon a point intersects, in other words polygons
## are not expected to overlap one another
r1 = cbind(c(0, 2, 2, 0, 0), c(0, 0, 2, 2, 0))
r2 = r1 + 1
sr1=Polygons(list(Polygon(r1)),"r1")
sr2=Polygons(list(Polygon(r2)),"r2")
sr=SpatialPolygons(list(sr1,sr2))
srdf=SpatialPolygonsDataFrame(sr, data.frame(cbind(1:2,3:4), row.names=c("r1","r2")))
proj4string(srdf) = CRS("+proj=longlat +datum=WGS84")

plot(srdf, axes=T)
pts = cbind(c(1.5), c(1.5))
pts = SpatialPoints(coords = pts)
over(pts, srdf)
srdf@data

