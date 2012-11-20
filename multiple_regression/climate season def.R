colnames(clim)
 [1] "PDSIavg"  "SodTemp"  "SPI1"     "SPI12"    "SPI24"    "rain1"   
 [7] "rain2"    "rain3"    "rain4"    "rain5"    "rain6"    "rain7"   
[13] "rain8"    "rain9"    "rain10"   "rain11"   "rain12"   "temp1"   
[19] "temp2"    "temp3"    "temp4"    "temp5"    "temp6"    "temp7"   
[25] "temp8"    "temp9"    "temp10"   "temp11"   "temp12"   "sum.rain"
[31] "win.rain" "spr.rain" "sum.temp" "win.temp" "spr.temp" "rain.avg"
[37] "temp.avg"

sum<-6:8
fal<-9:11
win<-c(12,1,2)
spr<-3:5

rain.sum<-0
rain.fal<-0
rain.win<-0
rain.spr<-0
temp.sum<-0
temp.fal<-0
temp.win<-0
temp.spr<-0


for(i in 1:3){
 rain.sum<-rain.sum+clim[,5+sum[i]]/3
 rain.fal<-rain.fal+clim[,5+fal[i]]/3
 rain.win<-rain.win+clim[,5+win[i]]/3
 rain.spr<-rain.spr+clim[,5+spr[i]]/3
 temp.sum<-temp.sum+clim[,17+sum[i]]/3
 temp.fal<-temp.fal+clim[,17+fal[i]]/3
 temp.win<-temp.win+clim[,17+win[i]]/3
 temp.spr<-temp.spr+clim[,17+spr[i]]/3
}


full.mod<-lm(resids~rain.sum*rain.fal*rain.win*rain.spr*
 temp.sum*temp.fal*temp.win*temp.spr)
summary(full.mod)

step.lm<-step(full.mod)
summary(step.lm)

full.mod2<-lm(resids~rain.sum+rain.fal+rain.win+
 temp.sum+temp.win+temp.spr+rain.sum:rain.win)
summary(full.mod2)

