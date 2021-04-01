#######################
### data subsetting ###
#######################
# march 2020

set.seed(42)
library(data.table)

#specify input output directories
indir<-"Q:/Querprojekte/aquaprobe/2.Lab_cytobuoy data/2020_beer_coasters/processed_sample_data"
        # indir<-"Q:/Querprojekte/aquaprobe/2.Lab_cytobuoy data/2020_beer_coasters/blanks_processed"
outdir<-"Q:/Querprojekte/aquaprobe/2.Lab_cytobuoy data/2020_beer_coasters"

filenames <-list.files(indir,pattern=c('Listmode'),recursive=FALSE)

#Make empty data table to populate with values from each data file
numberoffiles<-length(filenames)
numberofpoints<-100000L
pointsperfile<-as.integer(numberofpoints/numberoffiles)
numberofparameters<-94L #raw

m<-matrix(1,nrow=(numberoffiles*pointsperfile),ncol=numberofparameters)
dt <- as.data.table(m)

#Make the first row (id) an integer to avoid warning messages
dt <- dt[, V1:=as.integer(V1)]
# dt <- dt[, V98:=as.integer(V98)]

#Removing the matrix
rm(m)

#Create a list of row IDs where each new file's values will start at
idseq<-seq(1,(numberoffiles*pointsperfile),pointsperfile)

#Loop through all files to extract the data and populate the data table dt
setwd(indir)

for (index in 1:numberoffiles){
  print(index)
  # fread is a fast read function in the data.table package.
  # It does not read the dropped columns into R, increasing file reading speed considerably.
  # I have dropped the most unhelpful columns during the data processing stage, but if you want to exclude more, 
  # you will need to change the 'numberofparameters' value. 
  # Other ones that  you might want to drop: 'id', which is not really useful except to identify the point in the original file once again, 
  # and 'SWS Time of Arrival', which is probably only  useful if many samples were measured in combination as a single file with time lags between samples.

#Read data file
  dat<-fread(filenames[index], drop=c('V1'))
  #get rid of rare cases where total fluorescence=0, since they are obviously junk and we will log transform the data later
  dat2=dat[dat$`FL Red Total`>0]
  #Impose other needed conditions here. E.g. below excludes points with SWS.Length<2. Good practice to change the variable name (so that it is no longer dat2) if you are going to do this.
  # dat2=dat2[dat$SWS.Length>=2]
  #take a random subset from the file being read
  dat3=dat2[sample(.N, pointsperfile)]

  #populate existing data table with values
  set(dt, i=idseq[index]:(idseq[index]+pointsperfile-1), j=1:(numberofparameters), dat3)

}

rm(index)
#Renames the (uninformative) column names in dt to those from the last data file
setnames(dt, names(dt), names(dat))

# Write data subset
setwd(outdir)
write.csv(dt,paste0('alldatasubset_',numberofpoints,'beercoasters.csv'))

        # write.csv(dt,paste0('alldatasubset_',numberofpoints,'beercoasterblanks.csv'))

### END ###