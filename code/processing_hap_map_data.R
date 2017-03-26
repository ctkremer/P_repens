######################################

#### This code takes the output of the UNEAK pipeline and processes it for further analyses ####

# Original input: "HapMap.hmn.txt" (the alleles called in each individual) and "HapMap.hmc.txt" (the number of reads on which each call is based)

# OUTPUT: 
#   loci_20_v2.csv
#   all_allele_frequencies_032417.csv

######################################

#### Index of Script ####
#    1. Quality control and identification of markers to be used for further analyses
#    2. Remove outliers for phylogeography analyses (MDS, TREEMIX, STRUCTURE)
#    3. Make Structure Files for the non-outlier loci
#    4. Make Structure Files for the outlier loci
#    5. Calculate the frequency of the alleles in each population 

######################################

#### 1. Quality control and identification of markers to be used for further analyses ####

## Load data:

# This data set tracks which alleles are observed in which individuals across markers.
data<-read.delim("./data/HapMap.hmn.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = FALSE, comment.char = "",stringsAsFactors=F)

# This data set tracks how many reads were collected for each combination of individual and marker
counts=read.delim("./data/HapMap.hmc.txt", header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = FALSE, comment.char = "",stringsAsFactors=F)


## Clean data:

# Remove "empty" column and misc columns
data<-data[,-c(2:11, 242)]
counts<-counts[,-c(232, 719:723)]
#dim(data)
#dim(counts)


## Remove calls from counts for which there are fewer than 5 reads per marker/individual

# This helper function breaks up the terrible data format in count into usable numbers
# **** IT ALSO forces NA's on the results of the string split that are less than 5 ****
strsplitsum<-function(xx){
  res<-sapply(xx,function(x) strsplit(x,"|",fixed=T)[[1]])
  res<-as.numeric(res[1,])+as.numeric(res[2,])
  res<-ifelse(res<5,NA,res)
  res
}
# test it:
# strsplitsum(counts[1:100,2])


## Make a new dataframe with only the counts that are higher than 5

# Extract the number of reads from the counts data set:
newcount<-matrix(NA,nrow(counts),ncol(counts))
for(i in 2:717){
  print(i/715)
  xs<-strsplitsum(counts[,i])
  newcount[,i]<-xs
}
newcount[,1]<-counts[,1]  # fill in locus names in first column
colnames(newcount)<-colnames(counts)
#dim(newcount)
#newcount[1:5,1:4]


## Count up the total number of reads per individual
indiv_counts<-apply(newcount[,2:ncol(newcount)], 2, FUN=sum, na.rm=T)
#head(indiv_counts)
#hist(log(indiv_counts),50,col='blue')


## Remove individuals that do not have 500 reads from both counts and data results in 663 individuals

# These are the names of individuals with more than 500 total reads (plus the index column)
high.names<-unique(c(names(newcount[1]),names(indiv_counts[indiv_counts > 500])))

# Thin the data and count matrices to contain only these focal individuals
newcount<-newcount[,high.names]
newdata<-data[,high.names]
#dim(data)
#dim(newdata)
#ncol(data)-ncol(newdata)  # this has effectively removed 53 individuals.


## Remove calls from DATA for which there are fewer than 5 reads per marker/individual

# This next chunk transfers NA's from newcount into newdata, creating newnewdata
# such NA's arise from the occurence of < 5 reads.
newnewdata<-matrix(NA, nrow(newdata),ncol(newdata))
for (i in 2:ncol(newdata)){
  print(i/664)
  xs<-ifelse(is.na(newcount[,i]),".",as.character(newdata[,i]))
  newnewdata[,i]<-xs
}
newnewdata[,1]<-data$rs
newnewdata<-as.data.frame(newnewdata)
names(newnewdata)<-names(newdata)

# Check how many calls were lost?
helper<-function(x){
  length(x[x!="."])
}

tally1<-apply(newdata[,2:ncol(newdata)],2,helper)
tally2<-apply(newnewdata[,2:ncol(newnewdata)],2,helper)

# This is how many calls were dropped due to the thinning criteria (>= 5 reads)
#sum(tally1)-sum(tally2)
# [1] 2250218



## Separate individuals that have at least one marker scored. 
plate1<-newnewdata[,apply(newnewdata, 2, function(x) sum(x!= ".", na.rm=T)>0)]
#ncol(newnewdata)
#ncol(plate1)


## Count up how many individuals each marker worked in

# note: for each cell, the %in% returns (TRUE/1) if the cell has 0, 1, or 2, and (FALSE/0) for anything else.  Sum then adds those 0's and 1's up, effectively counting the number of individuals for which a marker works
worked<-apply(newnewdata,1,FUN=function(x) sum(x %in% c("0","1","2")))
worked<-as.vector(worked)
names(worked)<-newnewdata[,1]

## Separate out those markers that are represented in 20% of individuals
worked20<-worked[worked>0.2*(ncol(plate1)-1)]
plate120<-plate1[worked>0.2*(ncol(plate1)-1),]
#names(worked20)
#hist(worked20)
#nrow(plate120)


## Tidy up the results:
plate120T<-t(plate120)   #transposing dataframe
colnames(plate120T)<-plate120T[1,]   #naming the columns
plate120T<-plate120T[-which(row.names(plate120T)=="rs"),]   #removing the row with the marker names
#nrow(plate120T)
#dim(plate120T)


## Write out the resulting data file
write.table(plate120T, "./data/loci_20_v2.csv", quote=F, row.names=T, col.names=T, sep=",")





#### 2. Remove outliers for phylogeography analyses (MDS, TREEMIX, STRUCTURE) ####

## Load data
# data<-plate120T   # If continuing script from 1. above
# data<-read.csv("./data/loci_20_v2.csv",stringsAsFactors=FALSE) # If starting script here


## Format data
data<-data.frame(names=row.names(data), data)
row.names(data)<-NULL
#data[1:5,1:10]


## Load names of outliers
outliers<-read.csv("./data/outliers_columnnumbers.csv")
outliers<-data.frame(lapply(outliers, as.character), stringsAsFactors=FALSE)
#str(outliers)


## Extract names and indices of outlier markers

outliers$reallyright<-as.numeric(outliers$right.column)+1

## Subset data to: A) remove all outliers, and B) include only high outliers

# A) with no outliers 
noouts<-data[, -outliers$reallyright]

# save the no outlier data for treemix analysis
#write.table(noouts, "./data/noouts_data.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


# B) high outliers
high_outs<-data[,outliers$reallyright[outliers$level=="high"]]
high_outs<-data.frame(names=data$names, high_outs)

# save high outlier data (if desired)
#write.table(high_outs, "./data/high_outs.csv", quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)



#### 3. Make structure files for nonoutliers.  These will be used both for structure and for arlequin (after being read through PGD spider) ####
 
## Continue from above, or read in:
#noouts<-read.csv("noouts_data.csv",stringsAsFactors=F)
bob<-noouts

## Format data:
row.names(bob)<-bob$name  # set row names
bob<-bob[,-1] # remove name column
#bob[1:5, 1:10]

## Change missing data code from "." to "-18"  
for (i in 1:ncol(bob)) {	
  bob[,i]<-as.numeric(gsub(".", "-18", as.character(bob[,i]), fixed=TRUE))
}
#str(bob)
  
## Make two new dataframes

# part 1:
# 0 = 0
# 1 = 1
# 2 = 1
bobA<-data.frame(ceiling(bob/2))
bobA<-data.frame(name=as.character(row.names(bobA)),bobA)
row.names(bobA)<-NULL
#head(bobA)
  
# part 2:
# 0 = 0
# 1 = 0
# 2 = 1
bobB<-data.frame(floor(bob/2))
bobB<-data.frame(name=as.character(row.names(bobB)),bobB)
row.names(bobB)<-NULL
#nrow(bob)

## Make a column that has two times the row names
  
# glue the two together
bobs<-rbind(bobA,bobB)
bobs<-bobs[order(bobs$name),]
#bobs[1:10,1:20]
  
pop<-substr(bobs$name, 1, 3)
bobs<-cbind(pop, bobs)
bobs$pop<-as.numeric(as.factor(bobs$pop))
bobs<-bobs[,c(2,1,seq(3,ncol(bobs),1))]
#names(bobs)[1:5]
#summary(bobs$pop)
#bobs[1:5,1:10]
#ncol(bobs)

## Save output:  
#write.table(bobs, file="./data/p_repens_noouts_v2.str", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="	")


  
  
 #### 4. Make structure files for the outlier loci.  These will be used for arlequin (after being read through PGD spider) ####
  
## Continue from above, or read in:
#high_outs<-read.csv("high_outs.csv",stringsAsFactors=F)
bob<-high_outs

## Format data:
row.names(bob)<-bob$name
bob<-bob[,-1]
#bob[1:5, 1:10]

## Change missing data code from "." to "-18" 
for (i in 1:ncol(bob)) {	
  bob[,i]<-as.numeric(gsub(".", "-18", as.character(bob[,i]), fixed=TRUE))
}
#str(bob)
  
## Make two new dataframes

# part 1:
# 0 = 0
# 1 = 1
# 2 = 1
bobA<-data.frame(ceiling(bob/2))
bobA<-data.frame(name=as.character(row.names(bobA)),bobA)
row.names(bobA)<-NULL
#head(bobA)
  
# part 2:
# 0 = 0
# 1 = 0
# 2 = 1
bobB<-data.frame(floor(bob/2))
bobB<-data.frame(name=as.character(row.names(bobB)),bobB)
row.names(bobB)<-NULL
#nrow(bob)

## Make a column that has two times the row names
  
# glue the two together
bobs<-rbind(bobA,bobB)
bobs<-bobs[order(bobs$name),]
#bobs[1:10,1:20]
#dim(bobs)

pop<-substr(bobs$name, 1, 3)
bobs<-cbind(pop, bobs)
bobs$pop<-as.numeric(as.factor(bobs$pop))
bobs<-bobs[,c(2,1,seq(3,ncol(bobs),1))]
#names(bobs)[1:5]
#summary(bobs$pop)
#bobs[1:5,1:10]
#ncol(bobs)
  
## Save output:  
#write.table(bobs, file="./data/p_repens_high_outs_v2.str", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="	")
  



#### 5. Calculate the frequency of the alleles in each population.  This creates the imput for the effect size vs FST analysis ####

## Load tools
library(dplyr)
  
## Read in loci data, produced from earlier parts of script:
data<-read.csv("./data/loci_20_v2.csv",stringsAsFactors = F)
ind.id<-row.names(data)

## helper function for cleaning population codes  
splt<-function(x){strsplit(x,"_")[[1]][1]}

## process population codes
pop.id<-as.vector(sapply(ind.id,splt))

## format data  
data<-data.frame(pop.id,data)

## function for calculating allele frequencies
get.freq<-function(x){
  x<-x[x!="."]
  if(length(x)>0){
    out<-round(1-sum(as.numeric(x))/(2*length(x)),5)
  }else{
    out<-0
  }
  return(out)
}

## Calculate frequencies, cycling across populations and markers/loci
upops<-unique(pop.id)
freq.mat<-matrix(NA,nrow=length(upops),ncol=(ncol(data)-1))
for(i in 2:ncol(data)){
  tmp<-data[,c(1,i)]
  names(tmp)[2]<-"loci"
  freq<-tmp %>% group_by(pop.id) %>% summarize(freq=get.freq(loci))
  freq.mat[,i-1]<-freq$freq
}
allele_freq<-data.frame(pop=freq$pop.id,freq.mat)
names(allele_freq)[2:ncol(allele_freq)]<-names(data)[2:ncol(data)]

## Save allele frequency results as intermediate data file:
#write.csv(allele_freq,"./data/all_allele_frequencies_032417.csv",row.names=F)

