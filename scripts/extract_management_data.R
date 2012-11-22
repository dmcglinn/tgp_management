
## Purpose
## to create a burn matrix 
## rows = one fire event in a specific plot
## columns = plot, date, 

library(sp)

setwd('~/Lab data/tgp_management/')

load('./data/tgp_shpfiles.Rdata')

plots = read.csv('./data/tgp_plot_utm.csv')
plots_sp = SpatialPointsDataFrame(plots[ , c('easting', 'northing')],
                                  data = plots)
proj4string(plots_sp) = CRS(proj4string(pasture))


plot(pasture)
points(plots$easting, plots$northing, pch=19, col='red')

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
        date_of_burn = burns[[j]]@data[k, date_column[j]]
        if (flag == 0) {
          plot_burn = data.frame(plot=plots$plot[i], date_of_burn)
          flag = 1 
        }  
        else {
          plot_burn = rbind(plot_burn, 
                            data.frame(plot=plots$plot[i], date_of_burn))
        }  
      }
    }
  }  
}  

## this sort suggests these results are plausible
sort(table(plot_burn[,1]), dec=T)

as.Date('20050829', "%Y%m%d")
as.Date('2005/08/29', "%Y/%m/%d")

sum(sapply(tst, function(x) sum(!is.na(x[,1]))))
nrow(plot_burn)
## so it appears that there were 3 burns that 
## overlapped in a single year period thus justifing the 
## more complex and inefficient looping procedure
  
point.in.polygon(plots$easting[1], plots$northing[1], 
                 xy[,1], xy[,2])


xy = burns[[1]]@polygons[[1]]@Polygons[[1]]@coords
                 
str(burns[[1]])
                 
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

