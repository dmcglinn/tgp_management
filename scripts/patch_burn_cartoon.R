# the purpose of this script is to consider an idealized outcome of a patch-burn
# management scheme The setup: there are six patches each sampled across six
# years each patch is burned once there are S species each of which prefers a
# given time since burn (i.e., time since burn is the gradient). The species
# optima are evenly spaced along this gradient What we learn: site and year
# dummy variables explain 0 to very little variance in composition years since
# burn explains a large portion of the variance however, the degree to which
# this is true depends on how specialized species are on a particular time since
# burn if they have very narrow niches then time since fire explains about 40%. 
# This increases to about 80% when species have unimodal responses to time since fire

library(vegan)
library(dummies)

gauss.niche <- function(m,z,u,s){
  ##Purpose: to provide an exponential unimodal function 
  ##which is a model for a species response to the enviornment
  ##Called within the functions 'peaks' and 'sim.init.uni'
  ##from Palmer 1992
  ##Arguments:
  ##'m' is max perform
  ##'z' is enviornment at given coordinate
  ##'u' is env optim
  ##'s' is habitat breadth
  m * exp(-0.5 * (z - u)^2 / s^2)
}


setup_env = function(site, yr) {
  env = as.data.frame(expand.grid(yr, site))
  names(env) = c('yr', 'site')
  env$yr = as.factor(env$yr)
  env$site = as.factor(env$site)
 
  icount = 1
  YrsSLB = NULL
  for (i in site) {
      tmp = i - 1
      for (j in yr) {
          if (j > 1) {
              tmp = tmp + 1
              tmp = ifelse(tmp > 5, 0, tmp)
          }
          YrsSLB[icount] = tmp
          icount = icount + 1
      }
  }
  env$YrsSLB = YrsSLB
  return(env)
}

get_BP5Yrs = function(env, n_patches) {
  BP5Yrs = rep(NA, nrow(env))
  icount = 1
  sites = unique(env$site)
  yrs = unique(env$yr)
  for (i in sites) {
     env_sub = subset(env, site == i)
     for (j in yrs) {
         start_index = j - n_patches + 1
         start_index = ifelse(start_index < 1, 1, start_index)
         BP5Yrs[icount] = sum(env_sub$YrsSLB[start_index:j] == 0)
         icount = icount + 1
     }
  }
  return(BP5Yrs)
}


burn_type = c("reg", "ran") 
niche_width = c(0.01, 0.1, 1, 10)
analysis_type = c("rda", "cca")

n_patches = 6
site = 1:n_patches

n_years = 6
yr = 1:n_years

S = 20
max_abu = 10

results = data.frame()
icount = 1
for (n in niche_width) {
  for (b in burn_type) {
    for (a in analysis_type) {
      env = setup_env(site, yr)
      if (b == "ran") {
        env$YrsSLB = sample(env$YrsSLB)
        #env$BP5Yrs = get_BP5Yrs(env, n_patches)
      }
      sp_optima = seq(0, n_patches - 1, length.out = S)
      sp_mat = matrix(0, ncol=S, nrow=nrow(env))
      sp_mat = sapply(1:S, function(i) 
                      gauss.niche(max_abu, env$YrsSLB, sp_optima[i], n))
      richness = rowSums(sp_mat > 0)
      lmod = rda(richness ~ site + yr + YrsSLB, env)
      plmod = rda(richness ~ Condition(site) + Condition(yr) + YrsSLB, env)
      if (a == 'rda') {
        ord = rda(sp_mat ~ site + yr + YrsSLB, env)
        pord = rda(sp_mat ~ Condition(site) + Condition(yr) + YrsSLB, env)
      }
      if (a == 'cca') {
        ord = cca(sp_mat ~ site + yr + YrsSLB, env)
        pord = rda(sp_mat ~ Condition(site) + Condition(yr) + YrsSLB, env)
      }
      results[icount, 'niche_width'] = n
      results[icount, 'burn_type'] = b
      results[icount, 'analysis_type'] = a
      results[icount, 'r2_S_all'] = round(RsquareAdj(lmod)[[2]], 2)
      results[icount, 'r2_S_burn'] = round(RsquareAdj(plmod)[[2]], 2)
      results[icount, 'r2_comp_all'] = round(RsquareAdj(ord)[[2]], 2)
      results[icount, 'r2_comp_burn'] = round(RsquareAdj(pord)[[2]], 2)
      tst = anova(ord, by='terms')
      results[icount, 'Fsite'] = round(tst$F[1])
      results[icount, 'Fyr'] = round(tst$F[2])
      results[icount, 'Fburn'] = round(tst$F[3])
      icount = icount + 1
    }
  }
}
results
   niche_width burn_type analysis_type r2_S_all r2_S_burn r2_comp_all r2_comp_burn Fsite  Fyr  Fburn
1         0.01       reg           rda    -0.46     -0.06        0.17         0.57     0    0     18
2         0.01       reg           cca    -0.46     -0.06       -0.06         0.57     0    0      6
3         0.01       ran           rda     0.03     -0.04        0.41         0.39     2    2     18
4         0.01       ran           cca    -0.03     -0.04        0.02         0.39     1    1      2
5         0.10       reg           rda    -0.46     -0.06       -0.09         0.31     0    0      8
6         0.10       reg           cca    -0.46     -0.06       -0.16         0.31     0    0      6
7         0.10       ran           rda     0.13      0.00        0.22         0.23     1    1      8
8         0.10       ran           cca     0.13     -0.03        0.18         0.28     1    2      7
9         1.00       reg           rda      NaN       NaN        0.33         0.73     0    0     28
10        1.00       reg           cca      NaN       NaN        0.54         0.73     0    0     50
11        1.00       ran           rda      NaN       NaN        0.48         0.42     3    1     21
12        1.00       ran           cca      NaN       NaN        0.65         0.62     2    2     56
13       10.00       reg           rda      NaN       NaN        0.73         1.13     0    0    105
14       10.00       reg           cca      NaN       NaN        1.00         1.13     0    0 244155
15       10.00       ran           rda      NaN       NaN        0.81         0.83     4    6    110
16       10.00       ran           cca      NaN       NaN        1.00         0.72 19793 9174 256578

# one off variation partitioning
site_mat = dummy(env$site) 
year_mat = dummy(env$yr)
varpart(sp_mat, site_mat, year_mat, env$YrsSLB)


