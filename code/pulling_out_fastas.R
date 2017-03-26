############################################

#### Preparing fasta table for blasting ####

# This script processes sequence data in HapMap.fas.txt, extracting the sequences of the 2066 markers that passed QC, and assembling them into a new fasta file so these sequences can be blasted.

# Written by R. Prunier, updated by C.T. Kremer, 3/26/17

############################################

# Define helper function
is.odd <- function(x) x %% 2 != 0

#### Load data ####

# Data on prior classification of markers
markers<-read.csv("./data/marker_fst_table.csv")

# See data file on Dryad (too large for github)
fastas<-read.table("HapMap.fas.txt", sep=" ", colClasses="character")


#### Clean up loci names ####
f1<-function(x,n) strsplit(x,"_")[[1]][n]
fastas2<-unlist(sapply(fastas[,1],f1,n=1))
attr(fastas2,"names")<-NULL
#str(fastas2)


#### Make a dataframe from the fasta file ####
# with locus names in the first column and the corresponding sequence the second column
names=fastas2[seq(from =1, to=length(fastas2), by=4)]
seqs=fastas2[seq(from =2, to=length(fastas2), by=4)]
fasta=as.data.frame(cbind(names, seqs))
fasta$names<-as.character(fasta$names)
#head(fasta)


#### Extract the names of the different sets of markers ####

# Classified as low, high, or non-outliers in previous script, see 'marker' table...
high.name<-paste(">", markers$marker[markers$level2=="high"], sep="")
notout.name<-paste(">", markers$marker[markers$level2=="nonoutlier"], sep="")


#### Pull out the seqs corresponding to these sets of names ####

notout.seq<-fasta[fasta$names %in% notout.name,]
highout.seq<-fasta[fasta$names %in% high.name,]


#### Build a new fasta file for non outliers ####

nonout_fast=vector("logical", length=2*nrow(notout.seq))
check=vector("logical", length=2*nrow(notout.seq))

for(i in 1: length(nonout_fast)) {
	check[i]=ifelse(is.odd(i)==TRUE, (i+1)/2, i/2)
	nonout_fast[i]=ifelse(is.odd(i)==TRUE, notout.seq$names[(i+1)/2], as.character(notout.seq$seqs[(i/2)]))
}
#str(nonout_fast)


#### Save non outliers fasta file ####
#write.table(nonout_fast, "fasta_all_nonouts.txt", quote=FALSE, sep="	", row.names=FALSE, col.names=FALSE)


#### Build a new fasta file for high outliers ####

highout_fast=vector("logical", length=2*nrow(highout.seq))
check=vector("logical", length=2*nrow(highout.seq))

for(i in 1: length(highout_fast)) {
  check[i]=ifelse(is.odd(i)==TRUE, (i+1)/2, i/2)
  highout_fast[i]=ifelse(is.odd(i)==TRUE, highout.seq$names[(i+1)/2], as.character(highout.seq$seqs[(i/2)]))
}
#str(highout_fast)


#### Save high outliers fasta file ####
#write.table(highout_fast, "fasta_highouts.txt", quote=FALSE, sep="	", row.names=FALSE, col.names=FALSE)

