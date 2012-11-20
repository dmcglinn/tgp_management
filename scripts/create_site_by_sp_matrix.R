
setwd('~/Lab data/mangement/')

dat = read.csv('./data/Tgp_species_data_to_2009(all)_-_no_unk_sp_-_no_leaning_sp.csv')

cover_classes = read.csv('./data/cover_key.csv')

cov_vals = cover_classes$midcov[match(dat$cover, cover_classes$cov)]

comm = tapply(cov_vals, INDEX = list(dat$plot.yr, dat$name), sum)
comm = ifelse(is.na(comm), 0, comm)

comm_dat = data.frame(plot.yr = rownames(comm), comm)

write.csv(comm_dat, file='./data/tgp_comm_mat_all.csv', row.names=FALSE)