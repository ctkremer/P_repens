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


## Shultze values for all protea repens in the protea atlas dataset
env<-read.csv("./data/repens_7_env_1997.csv")

## Additional Schultze values for our sampled populations
samples<-read.csv("./data/sample_env_goldblatt.csv")

### Select columns for focal populations corresponding to the environmental variables in the larger P. repens environmental data set:
samples_1<-samples[,names(env)]
sample_names<-samples_1$name

### Combine data for the sampled locations to the Schultze values from the protea atlas dataset:
env<-rbind(env, samples_1)


# targets:
#pptcon
#MAP
#tmax_jan_a
#MAT
#altitude
#tmin july
#summer rainfall


### Standardize the seven environmental variables so that they can be used in the same PCA
env_standardized<-standardize.data.frame(env[,c(4:10)])
env_standardized<-cbind(env[,1:3], env_standardized)
#head(env_standardized)

### Pull out the names of the focal sample locations
#sample_names<-env_standardized[(nrow(env_standardized)-18):nrow(env_standardized),1]

### Remove rows that are missing data points
env_complete<-env_standardized[complete.cases(env_standardized==T),]
#nrow(env_standardized)-nrow(env_complete)  # 176 rows are dropped
#summary(env_complete)

### Run PCA on the seven environmental variables
env_pca<-princomp(env_complete[,c(4:10)], na.rm=T,scale=TRUE)
#env_pca<-prcomp(env_complete[,c(4:10)], na.rm=T,scale=TRUE)

### Pull out the scores of each location on the PCA axes
location_values<-env_pca[[6]]
loadings<-env_pca[[2]]
location_values<-as.data.frame(location_values)
#summary(location_values)

### Extract the axes matching the sampled locations
location_values<-data.frame(env_complete[,1:3],location_values)
sample_axes<-location_values[location_values$name %in% sample_names,]
sample_axes<-as.data.frame(sample_axes)

### Add additional covariates back in to this data set:
sample_axes<-data.frame(sample_axes, goldblatt=samples$goldblatt)
names(sample_axes)[1]<-"sample_names"

### Save this intermediate output, which supports subsequent analyses:
#write.table(sample_axes, file="sample_PCA_etc.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


### Error checking:

tmp<-read.csv("./data/sample_axes_1997.csv",stringsAsFactors = F)

c1<-tmp
c2<-sample_axes[,names(tmp)]

library(reshape2)

m1<-melt(c1,id.vars = 'sample_names')
m2<-melt(c2,id.vars = 'sample_names')
names(m2)[3]<-'new.value'
m3<-merge(m1,m2,by=c('sample_names','variable'))

plot(new.value~value,data=m3)
abline(0,1)
cor(x = m3$value,y=m3$new.value)
# 0.9884307

head(m3)
dim(m1)
head(m1)

dim(c1)==dim(c2)
i<-3
for(i in 1:ncol(c1)){
  plot(c1[,i],c2[,i])
  print(c1[,i]==c2[,i])
}

c1[,1]==c2[,1]
head(c1)
head(c2)

head(sample_axes)
as.character(tmp$sample_names)==as.character(sample_axes$sample_names)


as.character(tmp$sample_names)==as.character(sample_axes$sample_names)


