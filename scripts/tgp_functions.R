
read.shape = function(shpName, path=NULL) {
  require(rgdal)
  if (is.null(path))
    path = getwd()
  fileName = paste(shpName, '.shp', sep='')    
  shp = readOGR(file.path(path, fileName), shpName)
}  

r2_adj = function(Y, X, Z, reps, method, dummy=0) {
  ## Returns
  ## a vector of R2, R2adj, and all R2adj replicates that result from permuations
  ## Arguments
  ## Y, X, Z: are species matrix, expl matrix, covar matrix
  ## reps: the number of permutations to perform, 
  ##  if reps not specified the analytical r2 and/or r2adj is returned 
  ##  if method (see next) is rda then r2adj also returned
  ## method: specifies "cca" or "rda"
  ## dummy: a number 0, 1 or 2 depending on how many collinear variables are in the
  ##  explanatory matrix
  ##  dummy is only necessary for the analytical R2adj calculation
  Y = as.matrix(Y)
  X = as.matrix(X)
  if (missing(Z)) {
    cca.emp = eval(parse(text= paste(method, '(Y,X)')))
    r2 = summary(cca.emp)$constr.chi / cca.emp$tot.chi 
    if (missing(reps)) {
      if (method == 'rda') {
        n = nrow(Y)
        p = ncol(X) - dummy
        out = c(r2, 1 - (((n - 1) / (n - p - 1)) * (1 - r2)))
      }
      else 
        out = c(r2, NA)
    }
    else {
      if (reps <= 0)
        stop('reps argument must either be a positive integer or not specified')
      rand.r2 = rep(NA, reps)
      for(i in 1:reps){
        Xrand = X[sample(nrow(X)), ]
        rand.r2[i] = summary(eval(parse(text=paste(method, '(Y,Xrand)'))))$constr.chi
        if (i %% 100 == 0)
          print(i)
      }
      out = c(r2, 
              1 - (1 / (1 - mean(rand.r2 / cca.emp$tot.chi))) * (1 - r2),
              1 - (1 / (1 - rand.r2 / cca.emp$tot.chi)) * (1 - r2))
    }
  }  
  else{
    Z = as.matrix(Z)
    cca.emp = eval(parse(text=paste(method, '(Y,X,Z)')))
    r2 = summary(cca.emp)$constr.chi / cca.emp$tot.chi
    if (missing(reps)) {
      if (method == 'rda') {
        n = nrow(Y)
        p = ncol(X) - dummy
        out = c(r2, 1 - (((n - 1)/(n - p - 1)) * (1 - r2)))
      }
      else
        out = c(r2, NA)
    }
    else{
      if (reps <= 0)
        stop('reps argument must either be a positive integer or not specified')
      rand.r2 = rep(NA, reps)
      for(i in 1:reps){
        rhold = sample(nrow(X))
        Xrand = X[rhold, ]
        Zrand = Z[rhold, ]
        rand.r2[i] = summary( eval(parse(text=paste(method, '(Y,Xrand,Zrand)'))))$constr.chi
        if (i %% 100 == 0)
          print(i)
      }
      out = c(r2,
              1 - (1 / (1 - mean(rand.r2 / cca.emp$tot.chi))) * (1 - r2),
              1 - (1 / (1 - rand.r2 / cca.emp$tot.chi)) * (1 - r2))
    }  
  }
  return(out)
}

partition_r2 = function(full, X, Y, Z, X_Y, X_Z, X_YZ, Y_Z, Y_XZ, Z_XY,
                        adj=TRUE, digit=3) {
  ## Partition R2 values between two (XY) or three (XYZ) classes
  ## Returns:
  ## the independent and shared components of variation
  ## the two class partitioning is based on Legendre and Legendre 1998, p770-775
  ## the three class paritioning is based on Anderson & Gribble 1998
  ## Arguments:
  ## full: r2 for . ~ X + Y or . ~ X + Y + Z
  ## X : r2 for . ~ X
  ## Y : r2 for . ~ Y
  ## Z : r2 for . ~ Z
  ## X_Y : r2 for (. ~ Y) ~ X
  ## X_Z : r2 for (. ~ Z) ~ X
  ## X_YZ : r2 for (. ~ Y + Z) ~ X
  ## Y_Z : r2 for (. ~ Z) ~ Y
  ## Y_XZ : r2 for (. ~ X + Z) ~ Y
  ## Z_XY : r2 for (. ~ X + Y) ~ Z
  ## adj : boolean, if true then it expects adjusted r2 are also in the previous arguments
  ## digit : positive integer where to round the output table at
  ## Examples:
  ## from Legendre and Legendre (1998)
  ## partition_r2(.784, .450, .734, adj=F)
  ## from Anderson & Gribble (1998)
  ## partition_r2(.5050, .3467, .3772, .0794, .1073, .3004, .0889, .3367, .1252, .0205, adj=F, digit=6) * 100
  ## Citations:
  ## Anderson, M. J., and N. A. Gribble. 1998. Partitioning the variation among 
  ##    spatial, temporal and environmental components in a multivariate data set.
  ##    Austral Ecology 23:158–167.
  ## Legendre, P., and L. Legendre. 1998. Numerical ecology. Elsevier, Boston, Mass., USA.
  if (missing(Z)) {
    ## Legendre and Legendre (1998) p770-775
    abc = full
    ab = X
    bc = Y
    a = abc - bc
    c = abc - ab
    b = abc - a - c
    d = 1 - abc
    part = rbind(abc, a, b, c, d)
    rownames(part) = c('all', 'X indep.', 'X & Y shared', 'Y indep.', 'resid.')
  }
  else {
    ## Anderson & Gribble (1998)
    XYZ = X_YZ + (X - X_Y) + (X - X_Z) - X
    part = rbind(full, X_YZ, Y_XZ, Z_XY,
                 X - X_Y - XYZ,
                 X - X_Z - XYZ,
                 Y - Y_Z - XYZ,
                 XYZ, 1 - full)
    rownames(part) = c('all', 'X indep.' , 'Y indep.' ,'Z indep.',
                       'X & Y shared', 'X & Z shared', 'Y & Z shared',
                       'X & Y & Z shared', 'resid.')
  }
  if (adj)
    colnames(part) = c('R2','R2adj')
  else
    colnames(part) = c('R2')
  part = round(part, digit)
  return(part)
}

