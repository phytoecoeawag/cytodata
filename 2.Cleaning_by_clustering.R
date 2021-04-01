###############################################
### cleaning cytobuoy samples by clustering ###
###############################################
# march 2020

      #Install and load required packages
      
      list.of.packages <- c('scales', 'data.table','randomForest', 'plot3Drgl')
      new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
      if(length(new.packages)) install.packages(new.packages)
      
      #This package has a different source and needs to be installed separately
      source("http://bioconductor.org/biocLite.R")
      biocLite("flowPeaks")
      
      
          if (!requireNamespace("BiocManager", quietly = TRUE))
            install.packages("BiocManager")
          BiocManager::install(version = "3.10")
          BiocManager::install("flowPeaks")    
    
#Load packages

library(flowPeaks); library(scales);library(data.table);library(randomForest);library(plot3Drgl)

#Set seed for repeatability
set.seed(42)

#Set your working directory, which contains both the data subset and the files to be cleaned (the latter can be in subfolders)
indir<-"Q:/Querprojekte/aquaprobe/2.Lab_cytobuoy data/2020_beer_coasters"
outdir<-"Q:/Querprojekte/aquaprobe/2.Lab_cytobuoy data/2020_beer_coasters/20190305_cleaned_data"

#Read in the data subset to be used for clustering
setwd(indir)
# df<-read.csv('alldatasubset_100000beercoasters.csv')

        # I decided to use blanks to define junk
        df<-read.csv('alldatasubset_100000beercoasterblanks.csv')

# Define colour palette for clusters in 3D plots
palette1 <- c('black','red','green3','blue','cyan','magenta','yellow','gray',
              'brown','purple','pink','orange','darkgreen', 'goldenrod')

df <- data.table(df)
##Define new variables to help identify electonic noise
#These essentially quantify how different the first and last values are. 0.1 added because a) there
# are zeroes, b) it is just below the lowest non-zero measurement, c) we need to log transform later
df$FL.Red.Gradient <- abs(df$FL.Red.First - df$FL.Red.Last) + 0.1
df$FL.Orange.Gradient <- abs(df$FL.Orange.First - df$FL.Orange.Last) + 0.1
df$FL.Green.Gradient<- abs(df$FL.Green.First-df$FL.Green.Last)+0.1
df$FL.Yellow.Gradient<- abs(df$FL.Yellow.First-df$FL.Yellow.Last)+0.1

df <- df[df$FWS.Length > 0.2,]

#Calculate new Range variables
df[, FWS.Range := (FWS.Maximum - FWS.Minimum)]
df[, SWS.Range := (SWS.Maximum - SWS.Minimum)]
df[, FL.Orange.Range := (FL.Orange.Maximum - FL.Orange.Minimum)]
df[, FL.Red.Range := (FL.Red.Maximum - FL.Red.Minimum)]
df[, FL.Green.Range := (FL.Green.Maximum - FL.Green.Minimum)]
df[, FL.Yellow.Range := (FL.Yellow.Maximum - FL.Yellow.Minimum)]
df[, Curvature.Range := (Curvature.Maximum - Curvature.Minimum)]

#Note: we do not really use FWS.Length in clustering because size  estimates at small sizes (<5 microns) are dodgy right now. Worth taking a look at, though, so it's used here for visualization only
names(df)

#Select only the variables used for clustering and visualization
dat.clust<-df[,c('FWS.Range','FL.Red.Range','FL.Green.Range', 'FL.Orange.Range',
                 'FL.Yellow.Range', 'FL.Green.Gradient','FL.Red.Gradient','FL.Yellow.Gradient',
                 'FL.Yellow.Maximum','SWS.Length','FWS.Fill.factor')]

#Subset data to remove entirely ridiculous cases that will mess with clustering.
dat.clust <- subset(dat.clust,dat.clust$FWS.Range > 0)
dat.clust <- subset(dat.clust,dat.clust$FL.Red.Range > 0)
dat.clust <- subset(dat.clust,dat.clust$FL.Orange.Range > 0)
dat.clust <- subset(dat.clust,dat.clust$FL.Yellow.Range > 0)
dat.clust <- subset(dat.clust,dat.clust$FL.Green.Range > 0)

#The next two are more of a judgement call. Take a look at plots and clustering before trying them.
dat.clust<-subset(dat.clust,dat.clust$FWS.Range > 0.2)
dat.clust<-subset(dat.clust,dat.clust$SWS.Length > 1.5)

#Take the log of all except for FWS.Fill.factor, which ranges from 0 to 1
dat.clustlog<-data.frame(log10(dat.clust[,-c(11)]), dat.clust[,c(11)])
colnames(dat.clustlog)[11]<- 'FWS.Fill.factor'

#Cluster using the best parameters
names(dat.clustlog)
subflow <- flowPeaks(dat.clustlog[,c(1:5, 7)], tol = 0.5, h0 = 0.05, h = 2)

# examine details if interested. The 'weight' column indicatest the proportion of the data belonging to that cluster
summary(subflow)

#Examine any plots of interest
scatter3Drgl(dat.clustlog$FWS.Range, dat.clustlog$FL.Red.Range, dat.clustlog$FL.Orange.Range,
             colvar = subflow$peaks.cluster, col = palette1[1:max(subflow$peaks.cluster)],
             xlab = 'FWS.Range', ylab = 'FL.Red.Range', zlab = 'FL.Orange.Range')

scatter3Drgl(dat.clustlog$FL.Red.Range, dat.clustlog$FL.Orange.Range, dat.clustlog$FL.Red.Gradient,
             colvar = subflow$peaks.cluster, col = palette1[1:max(subflow$peaks.cluster)],
             xlab = 'FL.Red.Range', ylab = 'FL.Orange.Range', zlab = 'FL.Red.Gradient')

scatter3Drgl(dat.clustlog$FL.Red.Range, dat.clustlog$FL.Yellow.Range, dat.clustlog$FL.Orange.Range,
             colvar = subflow$peaks.cluster, col = palette1[1:max(subflow$peaks.cluster)],
             xlab = 'FL.Red.Range', ylab = 'FL.Yellow.Range', zlab = 'FL.Orange.Range')


#Visually identify junk clusters, which are characterized by low FWS.Average and low fluorescence. 
# High fluorescence on one channel but very low fluorescence on other channels can indicate one type of junk, electronic noise. Typically the first and second largest clusters (black and red colours) are junk, and some much smaller clusters as well. These are specific to each dataset and must be decided based on the plots. Once they are identified, use the cluster numbers to clean the dataset.

#Re-examine plot without the putative junk clusters (this may reveal smaller ones that were not
#visible earlier). Enter the numbers of your junk clusters below

    junk.clusterids <- c(1:2)

dat.clustlog2 <- data.frame(dat.clustlog,subflow$peaks.cluster)
dat.good <- subset(dat.clustlog2, !(dat.clustlog2$subflow.peaks.cluster %in% junk.clusterids))

# Examine plots of the cleaned data to identify more junk clusters
scatter3Drgl(dat.good$FWS.Range, dat.good$FL.Red.Range, dat.good$FL.Orange.Range,
              colvar = dat.good$subflow.peaks.cluster, col = palette1[1:max(subflow$peaks.cluster)],
              xlab = 'FWS.Range', ylab = 'FL.Red.Range', zlab = 'FL.Orange.Range')

# scatter3Drgl(dat.good$FL.Red.Range, dat.good$FL.Orange.Range, dat.good$FL.Yellow.Range,
#               colvar = dat.good$subflow.peaks.cluster, col = palette1[1:max(subflow$peaks.cluster)],
#               xlab = 'FL.Red.Range', ylab = 'FL.Orange.Range', zlab = 'FL.Yellow.Range')
# 
# scatter3Drgl(dat.good$FL.Red.Range, dat.good$FL.Orange.Range, dat.good$FL.Red.Gradient,
#               colvar = dat.good$subflow.peaks.cluster, col = palette1[1:max(subflow$peaks.cluster)],
#               xlab = 'FL.Red.Range', ylab = 'FL.Orange.Range', zlab = 'FL.Red.Gradient')


## Recluster cleaned data to see if additional junk clusters are now evident
# (uses a slightly different set of variables)
    # cleanflow <- flowPeaks(dat.good[,c(1:5, 7)], tol = 0.15, h0 = 0.05, h = 2)
    # 
    # #examine details if interested
    # summary(cleanflow)
    # 
    # # Plot reclustered data, coloured by new cluster identity
    # scatter3Drgl(dat.good$FL.Red.Range, dat.good$FL.Orange.Range, dat.good$FL.Yellow.Range,
    #              colvar = cleanflow$peaks.cluster, col = palette1[1:max(cleanflow$peaks.cluster)],
    #              xlab = 'FL.Red.Range', ylab = 'FL.Red.Range', zlab = 'FL.Yellow.Range')
    # 
    # 
    # # Specify any new junk clusters that emerge; leave blank if there are none
    # junk.clusterids.round2 <- c(6, 13, 14)
    # 
    # #Add cluster IDs to dataset
    # dat.clustlog3<-data.frame(dat.good,cleanflow$peaks.cluster)

# Train a random forest classifier to identify the clusters in both cleaning rounds

# First round of cleaning
# rf_cleaning1 <- randomForest(factor(subflow.peaks.cluster) ~ FWS.Range + FL.Red.Range + FL.Orange.Range +
#                                 FL.Red.Gradient, data = dat.clustlog2, ntree = 999, na.action = na.omit)
# setwd(indir)
# saveRDS(rf_cleaning1, 'Trained_RF_cleaning_2020_03_05.rds')
rf_cleaning1 <- readRDS('Trained_RF_cleaning_2020_03_05.rds')

    # rf_cleaning2 <- randomForest(factor(cleanflow.peaks.cluster) ~ FL.Red.Range + FL.Orange.Range +
    #                                FL.Yellow.Range + FL.Green.Range, data = dat.clustlog3,
    #                              ntree = 1001, na.action = na.omit)
    # saveRDS(rf_cleaning2, 'E:/Cytobuoy_data/Trained_RF_2015_data_cleaning2_2017_02_28.rds')
    # rf_cleaning2 <- readRDS('E:/Cytobuoy_data/Trained_RF_2015_data_cleaning2_2017_02_28.rds')

#Load the list of raw data files
indir<-"Q:/Querprojekte/aquaprobe/2.Lab_cytobuoy data/2020_beer_coasters/processed_sample_data"
outdir<-"Q:/Querprojekte/aquaprobe/2.Lab_cytobuoy data/2020_beer_coasters/20190305_cleaned_data"

setwd(indir)
filenames <-list.files(pattern=c('Listmode.csv'),full.names=FALSE,recursive=FALSE)

# Remove objects that are no longer needed
rm(cleanflow, subflow, dat.clust, dat.clustlog, dat.clustlog2, dat.clustlog3, dat.good,
   df, palette1)

    ### Classification function ###

    classifydat <- function(x){
      #Read in data file, exclude unhelpful columns unless needed
      setwd(indir)
      dat <- fread(x)
      names(dat) <- sub(" ", ".", names(dat))
      names(dat) <- sub(" ", ".", names(dat))
      colnames(dat)[which(colnames(dat) == "Particle.ID")] = "id"
      
      dat = dat[,FL.Red.Gradient := abs(FL.Red.First - FL.Red.Last) + 0.1]
      dat = dat[,FL.Orange.Gradient := abs(FL.Orange.First - FL.Orange.Last) + 0.1]
      
      # Calculate new Range variables
      dat[, FWS.Range := (FWS.Maximum - FWS.Minimum)]
      dat[, SWS.Range := (SWS.Maximum - SWS.Minimum)]
      dat[, FL.Orange.Range := (FL.Orange.Maximum - FL.Orange.Minimum)]
      dat[, FL.Red.Range := (FL.Red.Maximum - FL.Red.Minimum)]
      dat[, FL.Green.Range := (FL.Green.Maximum - FL.Green.Minimum)]
      dat[, FL.Yellow.Range := (FL.Yellow.Maximum - FL.Yellow.Minimum)]
      dat[, Curvature.Range := (Curvature.Maximum - Curvature.Minimum)]
      
      #Impose same conditions as at the beginning of this process, or clustering may go awry.
      dat = subset(dat,FWS.Range > 0)
      dat = subset(dat,FL.Red.Range > 0)
      dat = subset(dat,FL.Orange.Range > 0)
      dat = subset(dat,FL.Yellow.Range > 0)
      dat = subset(dat,FL.Green.Range > 0)
      
      dat = subset(dat, FWS.Length > 0.2)
      
      #Arrange by id. Done because the data.table package requires this in order to perform later
      #functions
      setkey(dat, id)
      
      #Select only the variables used for clustering and log-transform
      dat2 = dat[,list(FWS.Range,FL.Red.Range,FL.Green.Range, FL.Orange.Range,
                       FL.Yellow.Range, FL.Red.Gradient)]
      dat2 = log10(dat2)
      
      #Predict the cluster identity of particles using the random forests classifier (first round)
      dat[,classification1 := predict(rf_cleaning1, dat2)]
      
      #Remove the clusters that are believed to be junk
      dat3 = dat[!(classification1 %in% junk.clusterids)]
      
      ##Cleaning round 2
      #Select only the variables used for clustering and log-transform
      # dat4 = dat3[,list(FL.Red.Range, FL.Orange.Range, FL.Yellow.Range, FL.Green.Range)]
      # dat4 = log10(dat4)
      # 
      # #Add a new column with the predicted cluster identity
      # dat3[,classification2 := predict(rf_cleaning2, dat4)]
      # 
      # #Remove the clusters that are believed to be junk
      # dat5 = dat3[!(classification2 %in% junk.clusterids.round2)]
      # 
      
      #Remove classification columns if desired.
      dat5<-dat3
      dat5[, classification1 := NULL]
      # dat5[, classification2 := NULL]
      
      #Define RedOrange.ratio parameter
      dat5[, RedOrange.ratio := (FL.Red.Range / FL.Orange.Range)]
      dat5[, RedYellow.ratio := (FL.Red.Range / FL.Yellow.Range)]
      dat5[, RedGreen.ratio := (FL.Red.Range / FL.Green.Range)]
      
      #Order by id (i.e. original measurement sequence). Not necessary, but helpful for later error checking
      datfinal<-dat5
      setkey(datfinal, id)
      
      #Write cleaned file to .csv, file name edited by removing 'listmode_' and appending '_cleaned'
      #You might need to alter the digits within the 'substr' commands depending on how your files
      #are named
      
      # Cleaned files are written to different subfolder
      setwd(outdir)
      write.csv(datfinal, paste(substr(x, 1, nchar(x) - 4),'_cleaned','.csv',sep = ''),
                row.names = FALSE)
      
    }

    #Run function on all files. Can be run on multiple processors except on Windows
    mapply(classifydat, x=filenames)
    
### END ###
