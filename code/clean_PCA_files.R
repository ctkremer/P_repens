############################################

#### Preparing PCA axes based on seven environmental variables ####

# This script processes environmental data on seven variables () for locations where P. repens is observed across the region, including the focal 19 populations examined in more detail in this study. The resulting PCA axes are used in subsequent analyses exploring divergence among populations.

# Focal environmental variables are:
#   pptcon
#   MAP
#   tmax_jan_a
#   MAT
#   altitude
#   tmin july
#   summer rainfall

# OUTPUT: sample_PCA_and_covars.csv

# Written by R. Prunier, updated by C.T. Kremer, 3/26/17

############################################

#### Define helper functions ####
standardize <- function(x){
  if (is.numeric(x)) {
    y <- (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
  }else{
    y <- x
  }
  y
}

standardize.data.frame <- function(x) {
  x[] <- lapply(x, standardize)
  return(x)
}


#### Load data ####

# Shultze values for all protea repens in the protea atlas dataset
env<-read.csv("./data/repens_7_env_1997.csv")

# Additional Schultze values for our specific sampled populations
samples<-read.csv("./data/sample_env_goldblatt.csv")


#### Organize data ####

# Select columns for focal populations corresponding to the environmental variables in the larger P. repens environmental data set:
samples_1<-samples[,names(env)]
sample_names<-samples_1$name     # save the names of our focal populations

# Combine data for the sampled locations to the Schultze values from the protea atlas dataset:
env<-rbind(env, samples_1)


#### Standardize the environmental variables ####
# (so that they can be used in the same PCA)
env_standardized<-standardize.data.frame(env[,c(4:10)])   # first set of columns doesn't need to be standardized
env_standardized<-cbind(env[,1:3], env_standardized)    # restore these ID columns


### Remove rows that are missing data points ####
env_complete<-env_standardized[complete.cases(env_standardized==T),]
#nrow(env_standardized)-nrow(env_complete)  # 176 rows are dropped
#summary(env_complete)


#### Run PCA on the seven environmental variables ####
env_pca<-princomp(env_complete[,c(4:10)], na.rm=T,scale=TRUE)


#### Pull out the scores of each sample location on the PCA axes ####
location_values<-env_pca[[6]]
loadings<-env_pca[[2]]
location_values<-as.data.frame(location_values)
#summary(location_values)


#### Extract the axes matching the sampled locations ####
location_values<-data.frame(env_complete[,1:3],location_values)
sample_axes<-location_values[location_values$name %in% sample_names,]
sample_axes<-as.data.frame(sample_axes)


#### Add additional covariates back in to this data set: ####
sample_axes<-data.frame(sample_axes, goldblatt=samples$goldblatt)
names(sample_axes)[1]<-"sample_names"


#### Save this intermediate output, which supports subsequent analyses: ####
#write.table(sample_axes, file="./intermediate_data/sample_PCA_and_covars.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
