
setwd('~/tgp_management/scripts')

dir.create('../results')
dir.create('./log_files')

## perpare datasets ------------------------------------------------------------

system('Rscript tgp_mesonet_download.R > ./log_files/mesonet_download.log 2>&1',
       wait=FALSE)


## conduct analysis ------------------------------------------------------------
## variation partitioning 
## on grid plots
system('Rscript tgp_grid_varpart.R > ./log_files/grid_varpart.log 2>&1', 
       wait=FALSE)

## on repeat plots
system('Rscript tgp_repeat_varpart.R > ./log_files/repeat_varpart.log 2>&1', 
       wait=FALSE)

## test for significance of model terms and check for residual autocorrelation
## on grid plots
system('Rscript tgp_grid_testpart.R > ./log_files/grid_testpart.log 2>&1', 
       wait=FALSE)

## on repeat plots
system('Rscript tgp_repeat_testpart.R > ./log_files/repeat_testpart.log 2>&1', 
       wait=FALSE)
