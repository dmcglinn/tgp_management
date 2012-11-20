test<-read.table('C:/Users/Daniel McGlinn/Desktop/test.csv',sep=',',header=T)
test<-as.matrix(test[,2:13])

eigen(cor(test))

y<-rich1[plot==205]

summary(lm(y~test))

y<-runif(1000)
X<-matrix(NA,nrow=length(y),ncol=10)
for(i in 1:dim(X)[2]){
 X[,i]<-runif(length(y))
}
test<-lm(y~X)
summary(test)

summary(apply(X,1,sum))
summary(apply(X,1,mean))

eigen(cor(X))
eigen(cor(clim))



X<-matrix(NA,nrow=length(y),ncol=dim(clim)[2])
for(i in 1:dim(clim)[2]){
 X[,i]<-rnorm(length(y),mean=20,sd=10)
}
y<-X%*%runif(dim(X)[2])+rnorm(length(y),mean=2,sd=1)

test<-lm(y~X)
test<-lm(y~clim,singular.ok=FALSE)

