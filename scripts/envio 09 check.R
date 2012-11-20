endat<-read.table('C:/Users/hurlbertlab/Documents/Grad Ecology/TGP - vegmon/management manuscript/data folder/tgp env98-09.csv',sep=',',header=T)
attach(endat)
names(endat)
 [1] "plot"               "plotnum"            "plot.1"            
 [4] "yr"                 "plot_yr"            "date_samp"         
 [7] "jul_samp"           "easting"            "northing"          
[10] "grassht"            "forbht"             "woodyht"           
[13] "woodypct"           "waterpct"           "rockpct"           
[16] "barepct"            "slope"              "aspect"            
[19] "TEC"                "PH"                 "ORG"               
[22] "S"                  "P"                  "CA"                
[25] "MG"                 "K"                  "NA."               
[28] "BCA"                "BMG"                "BK"                
[31] "BNA"                "BH"                 "B"                 
[34] "FE"                 "MN"                 "CU"                
[37] "ZN"                 "AL"                 "rain6"             
[40] "rain7"              "rain8"              "rain9"             
[43] "rain10"             "rain11"             "rain12"            
[46] "rain1"              "rain2"              "rain3"             
[49] "rain4"              "rain5"              "temp6"             
[52] "temp7"              "temp8"              "temp9"             
[55] "temp10"             "temp11"             "temp12"            
[58] "temp1"              "temp2"              "temp3"             
[61] "temp4"              "temp5"              "OrBis"             
[64] "bison"              "YrsOB"              "BP5Yrs"            
[67] "YrsSLB"             "burn"               "date.of.burn"      
[70] "julian.day.of.burn"

diag.func<-function(x,grp){
 plot(unique(grp),tapply(x,grp,mean))
 arrows(unique(grp),tapply(x,grp,mean),unique(grp),tapply(x,grp,mean)+tapply(x,grp,function(x){sd(x)/sqrt(12)}),angle = 90)
 arrows(unique(grp),tapply(x,grp,mean),unique(grp),tapply(x,grp,mean)-tapply(x,grp,function(x){sd(x)/sqrt(12)}),angle = 90)
}


#CA looks fine
diag.func(log(CA),unique(yr))
diag.func(CA,unique(yr))

#
diag.func(log(FE),unique(yr))
diag.func(FE,unique(yr))

#
diag.func(PH,unique(yr))
#

diag.func(temp6,unique(yr))
diag.func(rain6,unique(yr))

boxplot(log(CA)~plot)
points(1:20,log(CA[yr==2009]),col='red')
#all in the lower quantile
boxplot(log(FE)~plot)
points(1:20,log(FE[yr==2009]),col='red')
#all in the upper quantile

boxplot(YrsSLB~plot)
points(1:20,YrsSLB[yr==2009],col='red')



