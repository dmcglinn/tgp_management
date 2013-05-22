

traits = read.csv('./data/species_traits.csv')

source('./scripts/tgp_repeat_data_import.R')

sp_codes = names(comm)

codes_nt = sp_codes[!(sp_codes %in% traits$shrt)]

sp_data = read.csv('./data/Tgp_species_data_to_2009(all)_-_no_unk_sp_-_no_leaning_sp.csv')

head(sp_data)

names_nt = as.character(sp_data$spname[match(codes_nt, sp_data$name)])

genera_nt = as.character(sapply(names_nt, function(x)
                                strsplit(x, ' ')[[1]][1]))
                          
genera_wt = sapply(as.character(traits$sp), function(x)
                   strsplit(x, ' ')[[1]][1])

families_nt = as.character(traits$family[match(genera_nt, genera_wt)])

need_trait = cbind(codes_nt, names_nt, families_nt,
                   NA, NA, NA)

need_trait[1 , 4:6] = c('N', 'A', 'F')
need_trait[2 , 3:6] = c('Scrophulariaceae', 'N', 'P', 'F')
need_trait[3 , 3:6] = c('Poaceae', 'N', 'P', 'C3')
need_trait[4 , 4:6] = c('N', 'P', 'F')
need_trait[5 , 4:6] = c('N', 'P', 'C3')
need_trait[6 , 3:6] = c('Asteraceae', 'N', 'A', 'F')
need_trait[7 , 3:6] = c('Apiaceae', 'N', 'A', 'F')
need_trait[8 , 3:6] = c('Asteraceae', 'N', 'A', 'F')
need_trait[9 , 3:6] = c('Asteraceae', 'N', 'P', 'F')
need_trait[10 , 3:6] = c('Fabaceae', 'N', 'T', 'W')
need_trait[11 , 3:6] = c('Asteraceae', 'N', 'A', 'F')
need_trait[12 , 4:6] = c('N', 'P', 'F')
need_trait[13 , 4:6] = c('N', 'P', 'C3')
need_trait[14 , 4:6] = c('N', 'P', 'C3')
need_trait[15 , 3:6] = c('Onagraceae', '', '', 'F')
need_trait[16 , 3:6] = c('Lamiaceae', 'N', 'P', 'F')
need_trait[17 , 3:6] = c('Boraginaceae', 'N', 'P', 'F')
need_trait[18 , 3:6] = c('Urticaceae', 'N', 'A', 'F')
need_trait[19 , 4:6] = c('N', 'T', 'W')
need_trait[20 , 4:6] = c('N', 'P', 'F')
need_trait[21 , 3:6] = c('Ranunculaceae', 'N', 'P', 'F')
need_trait[22 , 4:6] = c('N', 'S', 'W')
need_trait[23 , 3:6] = c('Asteraceae', 'I', 'A', 'F')
need_trait[24 , 4:6] = c('N', 'P', 'C4')
need_trait[25 , 3:6] = c('Lamiaceae', 'N', 'P', 'F')
need_trait[26 , 3:6] = c('Scrophulariaceae', 'I', 'P', 'F')

need_trait = as.data.frame(need_trait)
names(need_trait) = names(traits)

new_traits = rbind(traits, need_trait)

write.csv(new_traits, './data/species_traits.csv', row.names=F)



