
setwd('~/tgp_management/')

dir.create('./results')
dir.create('./log_files')

## variation partitioning of repeat plots

system('Rscript tgp_repeat_varpart.R > ./log_files/repeat_varpart.log 2>&1', 
       wait=FALSE)