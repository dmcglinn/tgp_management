library(okmesonet)

print('Begin downloading mesonet data-----------------------------------')

fora.mts = okmts(begintime="1996-06-01 00:00:00",
                 endtime="2009-05-31 23:55", station="fora")

fora.mts.avg = avgokmts(fora.mts, by="month")

write.csv(fora.mts.avg, file='../data/mesonet_mo_avg.csv', row.names=F)

print('Downloading mesonet data Compete!---------------------------------')
