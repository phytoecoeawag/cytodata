#Install and load required packages

# list.of.packages <- c('scales', 'data.table','randomForest', 'plot3Drgl')
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if (length(new.packages)) install.packages(new.packages)

#This package has a different source and needs to be installed separately
# source("http://bioconductor.org/biocLite.R")
# biocLite("flowPeaks")

#Load packages

library(flowPeaks)
library(lubridate)
library(data.table)
library(randomForest)
library(plot3Drgl)
library(dplyr)
library(dtplyr)

#Set seed for repeatability
set.seed(42)

# Create palette for cluster plots
# palette1 <- c('black','red','green3','blue','cyan','magenta','yellow','gray','brown','purple',
              # 'goldenrod','pink','orange','darkgreen', 'coral2', 'khaki4')

#Set your working directory, which contains both the data subset and the files to be cleaned (the latter can be in subfolders)
setwd('C:/Users/pomatifr/Dropbox/pomatifr/Documents/projects_docs/Eawag/Ponds/ETH_Lehmann')

#Read in the data subset to be used for clustering
# df <- read.csv('2015_Greifensee_cleandatasubset_100k_2017_03_02.csv')

#Select only the variables used for clustering

# dat.clust <- df[,c('RedOrange.ratio', 'RedGreen.ratio', 'RedYellow.ratio',
#                    'FL.Yellow.Range', 'FL.Orange.Range', 'FL.Red.Range', 'FL.Green.Range',
#                    'FL.Orange.Gradient','FL.Yellow.Fill.factor')]

#Take the log of all except for FWS.Fill.factor, which ranges from 0 to 1
# dat.clustlog <- data.frame(log10(dat.clust))
#Cluster
# subflow <- flowPeaks(dat.clustlog[,c(1:9)], tol = 0.25, h0 = 0.05, h = 2)
#examine details if interested
# summary(subflow)

#Visually identify junk clusters, which are characterized by low FWS.Average, high FWS.Fill.factor, and low fluorescence. High fluorescence on one channel but very low fluorescence on other channels can also indicate junk, specifically electronic noise. Typically the first and second largest clusters (black and red colours) are junk, and some much smaller clusters as well. These are specific to each dataset and must be decided based on the plots. Once they are identified, use the cluster numbers to clean the dataset.

# scatter3Drgl(dat.clustlog$RedOrange.ratio, dat.clustlog$RedYellow.ratio, dat.clustlog$RedGreen.ratio,
#              colvar = subflow$peaks.cluster, col = palette1[1:max(subflow$peaks.cluster)],
#              xlab = 'RedOrange.ratio', ylab = 'RedYellow.ratio', zlab = 'RedGreen.ratio')
#
# scatter3Drgl(dat.clustlog$FL.Red.Range, dat.clustlog$FL.Orange.Range, dat.clustlog$FL.Yellow.Range,
#              colvar = subflow$peaks.cluster, col = palette1[1:max(subflow$peaks.cluster)],
#              xlab = 'FL.Red.Range', ylab = 'FL.Orange.Range', zlab = 'FL.Yellow.Range')
#

# dat.clustlog$subflow.peaks.cluster <- subflow$peaks.cluster

#Train a random forest classifier to identify the clusters
# rf_groups <- randomForest(factor(subflow.peaks.cluster) ~ RedOrange.ratio + RedGreen.ratio + RedYellow.ratio + FL.Yellow.Range + FL.Orange.Range + FL.Red.Range + FL.Green.Range + FL.Orange.Gradient + FL.Yellow.Fill.factor, data = dat.clustlog, ntree = 1001, na.action = na.omit)
#
# saveRDS(rf_groups, '../Trained_RF_2015_grouping_cleaned_data_2017_03_02.rds')
#
#Read random forest for group classification
rf_groups <- readRDS('./Trained_RF_2015_grouping_cleaned_data_2017_03_02.rds')


#Load the list of cleaned data files
cleanfiles <- list.files(path = './20190305_cleaned_data/',
                         pattern = c('cleaned'), recursive = T, full.names = TRUE)

#Generate empty file to be populated with summary data
numberoffiles <- length(cleanfiles)
m <- matrix(1, nrow = (numberoffiles * length(rf_groups$classes)), ncol = 4)
dt <- as.data.table(m)
rm(m)

#Number of clusters
cluster_number <- length(rf_groups$classes)

#Fill in cluster number column
set(dt, i = NULL, j = 1L, rep(1L:cluster_number,length(cleanfiles)))

#Fill in column names
setnames(dt, names(dt), c('cluster','cluster_count','cluster_percent',
                          'sample_name'))

#Set non-numeric column types
dt <- dt[, cluster_count := as.integer(cluster_count)]
dt <- dt[, sample_name := as.character(sample_name)]

#Create list of starting row values for each file (i.e. start after every cluster has a row)
idseq <- seq(1,(numberoffiles * cluster_number), cluster_number)

#Create a list of all clusters for merging later (otherwise missing clusters are omitted)
clust_list <- data.frame(cluster = seq(1:cluster_number))


for (index in 1:numberoffiles) {

  print(index)

  #Read in data file
  dat <- fread(cleanfiles[index])
  names(dat) <- sub(" ", ".", names(dat))
  
  #Define RedOrange.ratio parameter
  dat[, RedOrange.ratio := (FL.Red.Range / FL.Orange.Range)]
  dat[, RedYellow.ratio := (FL.Red.Range / FL.Yellow.Range)]
  dat[, RedGreen.ratio := (FL.Red.Range / FL.Green.Range)]
  
  # select columns for group clustering (parameters based on training of RF for classification)
  dat1 <- dat[, c('RedOrange.ratio', 'RedGreen.ratio', 'RedYellow.ratio',
                                             'FL.Yellow.Range', 'FL.Orange.Range', 'FL.Red.Range',
                   'FL.Green.Range', 'FL.Orange.Gradient','FL.Yellow.Fill.factor','FWS.Length')]
  
  #Log-transform relevant variables and add a new one to classify with random forest
  dat2 = log10(dat1)

  #Predict the cluster identity of particles using the random forests classifier
  preds <- predict(rf_groups, dat2)

  #Add cluster IDs to dataset
  dat3 = dat[, cluster := preds]

  #save number of cells belonging to each cluster
  set(dt, i = idseq[index]:(idseq[index] + cluster_number - 1), j = 2L, summary(preds))

  #Save proportion of cells belonging to each cluster
  set(dt, i = idseq[index]:(idseq[index] + cluster_number - 1), j = 3L,
      summary(preds) / sum(summary(preds)))

  #Save sample information
  set(dt, i = idseq[index]:(idseq[index] + cluster_number - 1), j = 4L,
      rep((as.character(substr(cleanfiles[index], 25, 32))),
          cluster_number))

}

rm(index)

head(dt)

densdat <- read.csv('./beercoasters_density_dynamics_2020_03_07.csv')
head(densdat)

clustdat <- merge(dt, densdat, by = 'sample_name')
clustdat$dens_clust <- clustdat$dens_clean * clustdat$cluster_percent
head(clustdat)

write.csv(clustdat, paste('./beercoasters_group_dynamics_', format(Sys.time(), '%Y_%m_%d'),
                         '.csv', sep = ''), row.names = FALSE)

### ==================== ###
### which group is cyano ###
### ==================== ###

ds<-fread(cleanfiles[5]);head(ds)
### WRONG - i should use a subset of the clean data, or just one large clean file!!
dim(ds)
names(ds) <- sub(" ", ".", names(ds))

ds<-ds[log10(ds$FWS.Length)>0 & log10(ds$FWS.Length)<3, ] # kill outliers 

ds<-ds[,c('RedOrange.ratio', 'RedGreen.ratio', 'RedYellow.ratio',
                  'FL.Yellow.Range', 'FL.Orange.Range', 'FL.Red.Range',
                  'FL.Green.Range', 'FL.Orange.Gradient','FL.Yellow.Fill.factor','FWS.Length')] # input variales for classification, logged
dx<-log10(ds)

preds <- predict(rf_groups, dx)
summary(preds)
# barplot(summary(preds), main = "counts x clusters, with #1 likely cyanos")
## from Mridul's previous plots, it seems that cyanos is cluster 2, instead from below it's likely #1

ds<-data.frame(ds,preds); head(ds)

cyan<-ds[ds$preds == 1, ]
rest1<-ds[ds$preds == 2, ] # most common eukaryotic algal cluster
# rest2<-ds[ds$preds == 4, ] 
# rest3<-ds[ds$preds == 7, ]
# rest4<-ds[ds$preds == 8, ]

dim(ds); dim(cyan); dim(rest)

plot(ds$FL.Red.Range, ds$FL.Orange.Range, pch=NA,
     main = "green=cyanos (#1); brown = most common eukaryotic algae (#2)")
points(cyan$FL.Red.Range, cyan$FL.Orange.Range, col="green")
points(rest1$FL.Red.Range, rest1$FL.Orange.Range, col="brown")
# points(rest2$FL.Red.Range, rest2$FL.Orange.Range, col="red")
# points(rest3$FL.Red.Range, rest3$FL.Orange.Range, col="orange")
# points(rest4$FL.Red.Range, rest4$FL.Orange.Range, col="purple")




### END ###
