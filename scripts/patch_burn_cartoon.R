## the purpose of this script is to consider an idealized
## outcome of a patch-burn management scheme
## The setup:
## there are six patches each sampled across six years
## each patch is burned once
## there are six species each of which prefers a given
## time since burn (i.e., time since burn is the gradient).
## What we learn:
## site and year dummy variables explain 0 variance in composition
## years since burn explains a large portion of the variance 
## however, the degree to which this is true depends on how 
## specialized species are on a particular time science burn
## if they only completely occur in a given time since burn
## then at most this variable explains around 20% 
## this will increase to approx 80% if species show a unimodal
## relationship with their specific year of preference.

library(vegan)

site = rep(1:6, each=6)

yr = rep(1:6, 6)

burn_yr = NULL
for(i in 0:5) { 
  x = NULL
  for(j in i:5) {
    x = c(x, j)
  }
  if(length(x) < 6) {
    x = c(x, 0:(min(x) - 1))
  }
  burn_yr = c(burn_yr, x)
}  

sp_mat = matrix(0, ncol=6, nrow=length(burn_yr))
for(i in 0:5) {
  sp_mat[ , i + 1] = ifelse(burn_yr == i, 1, 0)  
}

gauss.niche<-function(m,z,u,s){
  ##Purpose: to provide an exponential unimodal function 
  ##which is a model for a species response to the enviornment
  ##Called within the functions 'peaks' and 'sim.init.uni'
  ##from Palmer 1992
  ##Arguments:
  ##'m' is max perform
  ##'z' is enviornment at given coordinate
  ##'u' is env optim
  ##'s' is habitat breadth
  m*exp(-.5*(z-u)^2/s^2)
}

sp_mat = matrix(0, ncol=6, nrow=length(burn_yr))
for(i in 0:5) {
  sp_mat[ , i + 1] = gauss.niche(1, burn_yr, i, 3) 
}

plot_id = sort(unique(site))
year_id = sort(unique(yr))
plot_mat = matrix(0, ncol=length(plot_id), nrow=length(burn_yr))
year_mat = matrix(0, ncol=length(year_id), nrow=length(burn_yr))
               
for(i in 1:length(burn_yr)) {
  plot_mat[i, match(site[i], plot_id)] = 1
  year_mat[i, match(yr[i], year_id)] = 1 
}

## drop first columns so no singular variables in models
plot_mat = plot_mat[ , -1]
year_mat = year_mat[ , -1]


ord = rda(sp_mat, burn_yr)
anova(ord)
plot(ord)

ord = rda(sp_mat ~ plot_mat + year_mat + burn_yr)
ord = cca(sp_mat ~ plot_mat + year_mat + burn_yr)
ord
anova(ord, by='terms')
plot(ord)



