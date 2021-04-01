### ================================== ###
### Summarizing + density cleaned data ###
### ================================== ###
# march 2020

library(data.table)
library(lubridate)

cleanfiles <- list.files(path = './20190305_cleaned_data/', pattern = c('cleaned'),
                         recursive = T, full.names = TRUE)

clean <- matrix(1, nrow = length(cleanfiles), ncol = 4)
cleansummary <- as.data.table(clean)
cleansummary <- cleansummary[, V1 := as.integer(V1)]
cleansummary <- cleansummary[, V4 := as.character(V4)]
rm(clean)

metadat <- list.files(path = './processed_info_files/', 
                      pattern = c('_Info'), recursive = T, full.names = TRUE)

for (index in 1L:length(metadat)) {
  print(index)
  cl <-       fread(cleanfiles[index], select = 1)
  fileinf <-  read.table(metadat[index], fill = TRUE, skip = 28)
  set(cleansummary, i = index, j = 1L, dim(cl)[1])
  set(cleansummary, i = index, j = 2L, as.numeric(as.character(fileinf[1,5]))) # total unclean counts
  set(cleansummary, i = index, j = 3L, as.numeric(as.character(fileinf[2,2]))) # density unclean (/uL)
  set(cleansummary, i = index, j = 4L, 
      as.character(substr(cleanfiles[index], 25, 32))) #  sample name
}

setnames(cleansummary, c('count_clean','count_unclean','dens_unclean','sample_name'))
# Change units from per microlitre to per millilitre
cleansummary$dens_unclean <- cleansummary$dens_unclean * 1000

cleansummary = cleansummary[,clean_percent := count_clean / count_unclean]
cleansummary = cleansummary[,dens_clean := clean_percent * dens_unclean]

densdat <- cleansummary
rm(cleansummary)

write.csv(densdat, paste('2015_Greifensee_density_dynamics_', format(Sys.time(), '%Y_%m_%d'),
                         '.csv', sep = ''), row.names = FALSE)
