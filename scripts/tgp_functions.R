
read.shape = function(shpName, path=NULL) {
  require(rgdal)
  if (is.null(path))
    path = getwd()
  fileName = paste(shpName, '.shp', sep='')    
  shp = readOGR(file.path(path, fileName), shpName)
}  

r2_adj<-function(Y,X,Z,reps,method,dummy=0) {
  ##purpose: returns R2, R2adj, and all R2adj replicates that result from permuations
  ##Y,X,Z are spdata, expl mat, covar mat
  ##reps the number of permutations to perform, if no reps then only r2 and r2adj returned 
  ##if reps is not provided then it calculates an analytical R2adj (only ok for RDA)
  ##method specifies "cca" or "rda"
  ##dummy is a number 0, 1 or 2 depending on how many collinear variables are in the explanatory matrix
  ##dummy is only necessary for the analytical R2adj calculation
  Y<-as.matrix(Y)
  X<-as.matrix(X)
  if(missing(Z)){
    cca.emp<-eval(parse(text= paste(method,'(Y,X)')))
    r2<-summary(cca.emp)$constr.chi/cca.emp$tot.chi 
    if(missing(reps)){
      n<-nrow(Y)
      p<-ncol(X)-dummy
      out = c(r2, 1-(((n-1)/(n-p-1))*(1-r2)))
    }
    else{
      rand.r2<-rep(NA,reps)
      for(i in 1:reps){
        Xrand<-X[sample(nrow(X)),]
        rand.r2[i]<-summary( eval(parse(text= paste(method,'(Y,Xrand)'))))$constr.chi
        if (i %% 100 == 0)
          print(i)
      }
      out = c(r2, 
              1-(1/(1-mean(rand.r2/cca.emp$tot.chi)))*(1-r2),
              1-(1/(1-rand.r2/cca.emp$tot.chi))*(1-r2))
    }
  }  
  else{
    Z<-as.matrix(Z)
    cca.emp<-eval(parse(text= paste(method,'(Y,X,Z)')))
    r2<-summary(cca.emp)$constr.chi/cca.emp$tot.chi
    if(missing(reps)){
      n<-nrow(Y)
      p<-ncol(X)-dummy
      out = c(r2,1-(((n-1)/(n-p-1))*(1-r2)))
    }
    else{
      rand.r2<-rep(NA,reps)
      for(i in 1:reps){
        rhold<-sample(nrow(X))
        Xrand<-X[rhold,]
        Zrand<-Z[rhold,]
        rand.r2[i]<-summary( eval(parse(text= paste(method,'(Y,Xrand,Zrand)'))))$constr.chi
        if (i %% 100 == 0)
          print(i)
      }
      out = c(r2,
              1-(1/(1-mean(rand.r2/cca.emp$tot.chi)))*(1-r2),
              1-(1/(1-rand.r2/cca.emp$tot.chi))*(1-r2))
    }  
  }
  return(out)
}