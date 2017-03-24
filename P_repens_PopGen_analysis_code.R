#################################################################################

####~~~~~~~~~~~~~~~~  P. repens Population Genetics Analysis ~~~~~~~~~~~~~~~~####

# Code by Colin T. Kremer & Rachel Prunier
# Updated 12/22/16
# Updated along with final revisions 3/24/17

# Written for:
#
#   "Isolation by distance and isolation by environment contribute to population differentiation 
#   in Protea repens (Proteaceae L.), a widespread South African species"
#         American Journal of Botany

# This code details the data processing and analyses supporting sections of the above paper, 
# chiefly the MCMCglmm analyses of the effects of physical and environmental distance on 
# between population genetic differentiation. You can either start at the beginning, or jump
# directly to sections 3, 5, or 7, and load intermediate data sets.

#################################################################################
####                    Index                                                ####
#################################################################################

# (1) Load libraries, data, and declare functions
#
# (2) Create pairwise population comparisons
#     a) genetic distance
#     b) environmental/physical distance
#     c) merge environmental and genetic distance
#
# (3) MCMCglmm analyses of markers
#     a) Set up MCMCglmm analyses
#     b) Run complete set of MCMCglmm analyses
#     c) Process output of MCMCglmm runs
#     d) Combine model coefficients with marker data
#     e) Focus in on a single marker for additional graphs
#
# (4) randomize allele frequencies for 5 null models
#
# (5) run MCMCglmm analyses on null models
#     a) Set up and run simulation #1
#     b) Set up and run simulation #2
#     c) Set up and run simulation #3
#     d) Set up and run simulation #4
#     e) Set up and run simulation #5
#     f) Combine simulation results with original data
#
# (6) fit GAMs and plot results
#     a) run GAM fits
#     b) plot GAM fits
#     c) GAM fit against original data
#     d) GAM fit against simulation1 data
#     e) Show distribution of high outliers & non-outliers vs FST
#     f) Save statistical output from GAM fits


#################################################################################
#### (1) Load libraries, data, and declare functions                         ####
#################################################################################

library(parallel)
library(introgress)
library(MCMCglmm)
library(fossil)
library(reshape2)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(dplyr)

# Load data on allele frequencies in all populations (rows) by all markers (columns):
afreq<-read.csv("all_allele_frequencies.csv",stringsAsFactors = F)
#afreq[1:5,1:6]

# Load data on pair-wise population comparisons, for each marker
#delts<-read.csv("all_marker_deltas_for_analysis.csv")
#delts[1:5,1:15]

# loci_20, markers by individuals
# provide allele frequencies
# Environmental/spatial data for each population (not pairwise)


# Load environmental PCA scores and other covariates for each population's location
covars.dat<-read.csv("./intermediate_data/sample_PCA_and_covars.csv")

# Split this into environmental/PCA variables
envaxes<-covars.dat[,c("sample_names","Comp.1","Comp.2","Comp.3","Comp.4","Comp.5","Comp.6","Comp.7")]

# ... and physical locations/phytogeographic province of populations 
locs<-covars.dat[,c('sample_names','LONDD','LATDD','goldblatt')]
names(locs)[1]<-"name"

# Load table of marker names, fst values, and outlier categories:
fst.tab<-read.csv("marker_fst_table.csv",stringsAsFactors = F)
head(fst.tab)


##############################

# helper function for processing distance matrices quickly:
make.longform<-function(xmat,nms,vname=NA){
  rownames(xmat)<-nms
  colnames(xmat)<-nms
  mt1<-melt(xmat)
  res.mat<-data.frame(pop.pair=paste(mt1$Var1,mt1$Var2),value=mt1$value)
  
  if(!is.na(vname)){
    names(res.mat)[2]<-vname
  }
  
  return(res.mat)
}

# Function for standardizing covariates by their mean and standard deviation
standardize <- function(x) {
  if (is.numeric(x)) {
    y <- (x - mean(x, na.rm=TRUE))/sd(x, na.rm=TRUE)
  } else {
    y <- x
  }
  return(y)
}


# Function for calculating the marginal and conditional R2 values given MCMCglmm objects
# - follows from work of Nakagawa and Schielzeth 2013
get.R2.MCMCglmm<-function(model){
  sm<-summary(model)
  
  pred.vals<-predict(model)
  
  # variance of fixed effects
  v.fixed<-var(pred.vals)
  
  # variance of random effects
  v.rand<-sm$Gcovariances[1]
  
  # residual variance
  v.resid<-sm$Rcovariances[1]
  
  #### Calculate R2's
  
  # marginal R2:
  R2.m<-v.fixed/(v.fixed+v.rand+v.resid)
  
  # conditional R2:
  R2.c<-(v.fixed+v.rand)/(v.fixed+v.rand+v.resid)
  
  vec<-data.frame(R2.m,R2.c)
  
  return(vec)
}

# Function for tidying up the lists produced by markerglmm() after running MCMCglmm models:
acc<-function(x,i){
  if(length(x)<i){
    res<-NA
  }else{
    res<-x[[i]]    
  }
  return(res)
}


#################################################################################
#### (2) Create pairwise population comparisons                              ####
#################################################################################

####      a) Genetic distances                                               ####

# Calculating delta between each pair of populations and for each marker
m<-ncol(afreq)-1  # of markers
n<-nrow(afreq)  # of populations
npairs<-sum(seq(1,n-1,1))  # number of unique pairings of 19 populations
delt<-matrix(nrow=npairs, ncol=m)  # data storage
pop.mat<-matrix(nrow=npairs,ncol=2)   # store population ID
for(d in 1:m){
  k<-0
  #print(d/m)
  for(i in 1:(n-1)){
    for(j in (i+1):n){	
      k<-k+1
      delt[k,d]<-delta(c(afreq[i,d+1], (1-afreq[i,d+1])), c(afreq[j,d+1], (1-afreq[j,d+1])))
      
      if(d==1){
        pop.mat[k,]<-c(afreq[i,1],afreq[j,1])
      }
    }
  }
}
delt<-as.data.frame(delt)
pop.mat<-as.data.frame(pop.mat)
delt<-data.frame(pop1=pop.mat[,1],pop2=pop.mat[,2],delt)
names(delt)[3:ncol(delt)]<-names(afreq)[2:ncol(afreq)]
delt<-data.frame(pop.pair=paste(delt$pop1,delt$pop2),delt)
#delt[1:5,1:5]

####      b) Environmental distances                                               ####

# Calculate environmental distances on each of the first three PCA axes 
# (which together explain 90% of variation)
distpc1<-dist(envaxes[,2], diag=TRUE, upper=TRUE)
distpc1<-as.matrix(distpc1)
pc1.mat<-make.longform(x=distpc1,nms=as.character(envaxes$sample_names),vname='distpc1.a')
pc1.mat$distpc1.a<-standardize(pc1.mat$distpc1.a)  # standardize covariate

distpc2<-dist(envaxes[,3], diag=TRUE, upper=TRUE)
distpc2<-as.matrix(distpc2)
pc2.mat<-make.longform(x=distpc2,nms=as.character(envaxes$sample_names),vname='distpc2.a')
pc2.mat$distpc2.a<-standardize(pc2.mat$distpc2.a)  # standardize covariate

distpc3<-dist(envaxes[,4], diag=TRUE, upper=TRUE)
distpc3<-as.matrix(distpc3)
pc3.mat<-make.longform(x=distpc3,nms=as.character(envaxes$sample_names),vname='distpc3.a')
pc3.mat$distpc3.a<-standardize(pc3.mat$distpc3.a)  # standardize covariate

# Calculate physical distances between populations
locsonly<-data.frame(locs$LONDD, locs$LATDD)
distkm<-earth.dist(locsonly, dist=F)	
dist.mat<-make.longform(x=distkm,nms=as.character(locs$name),vname='physdist')
dist.mat$physdist<-standardize(dist.mat$physdist)  # standardize covariate


# Score whether or not populations belong to same Goldblatt phytogeographic zone:
tmp<-locs[,c('name','goldblatt')]
names(tmp)[1]<-'pop1'
tmp2<-data.frame(pop2=sort(rep(tmp$pop1,nrow(tmp))),pop1=tmp$pop1)
tmp2<-merge(tmp2,tmp)
names(tmp)<-c('pop2','goldblatt2')
tmp2<-merge(tmp2,tmp)
tmp2$goldblattsame<-ifelse(tmp2$goldblatt==tmp2$goldblatt2,1,0)
gb.mat<-data.frame(pop.pair=paste(tmp2$pop1,tmp2$pop2),goldblattsame=tmp2$goldblattsame)


# Combine various distances into a set of covariates:
env.dists<-merge(pc1.mat,pc2.mat)
env.dists<-merge(env.dists,pc3.mat)
env.dists<-merge(env.dists,dist.mat)
env.dists<-merge(env.dists,gb.mat)

dim(env.dists)
head(env.dists)


# Check: Do these need to be tracked for the all outlier, nonoutlier analyses? fst.noout, fst.outlier


####      c) Merge genetic & environmental distances                         ####

delt2<-merge(env.dists,delt,all.y=T,by='pop.pair')
delt2<-delt2[,c(1,7,8,5,2,3,4,6,seq(9,ncol(delt2),1))]
delt2<-delt2[,-which(names(delt2)=="pop.pair")]
#dim(delt2)
#delt[1:5,1:5]
#delt2[1:4,1:12]

# Save the output for faster loading in later sections of analysis:
write.csv(delt2,"pairwise_population_differences.csv",row.names = F)


#################################################################################
#### (3) MCMCglmm analyses of markers                                        ####
#################################################################################

####      a) Set up MCMCglmm analyses                                        ####

# Load the pair-wise genetic distances and covariates for all 2066 markers:
delt<-read.csv("pairwise_population_differences.csv")
#delt[1:5,1:15]

# Goal: run an MCMCglmm model examining variation in genetic distances between populations
# as a function of environmental and physical distances between populations, for each of 2066
# different markers. This is a large amount of computation, so we take advantage of parallel
# processing to complete the analysis in a reasonable amount of time - ~12 hrs running 6 cores
# on a 2.8 GHz Intel Core i7 CPU. 

# To conveniently run analyses in parallel, we need a function that completes a full analysis 
# given the name of a specific marker that needs to be analyzed. Then we pass this function to
# parallel processing commands in R.

# Function for running MCMCglmm model for a given marker (whose index is specificed by 'x'):
markerglmm=function(x, nitt, burnin, thin){
  
  # Wrapper function for MCMCglmm run
  # - note that this containes a random effect tracking the identity of population pairs, to account
  # for a lack of independence in pairwise data sets, forming a 'multiple membership' mixed effects model
  MCMCglmm2<-function(fixed, data=data){
    test<-MCMCglmm(fixed=fixed, random= ~ idv(mult.memb(~ pop1+pop2)), data=data, 
                   pr=TRUE, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
    return(test)
  }
  
  # Pull out the name of marker/locus of focus, given x...
  nm<-names(delts)[x]
  
  # Fit model, invoking main effects of all covariates
  mod.direct<-eval(substitute(
              MCMCglmm2(fixed=tmp~goldblattsame+physdist+distpc1.a+distpc2.a+distpc3.a,data=delts),
              list(tmp=as.name(nm))))
  
  # list of models
  mods<-list(mod.direct)
  
  # Combine model results:
  # this is legacy code from running multiple models/conducting model comparison; the more important
  # section for current analyses follows under 'Model coefficient estimates'
  DICs.rad_outs <- data.frame(array(NA, c(1, 7)))
  names(DICs.rad_outs) <- c('model','outlier','DIC', 'deltaDIC', 'DICweight',
                            'R2.marginal','R2.conditional')
  DICs.rad_outs$model <- c('alldirect')
  DICs.rad_outs$outlier<-names(delts)[x]
  DICs.rad_outs$DIC<-unlist(lapply(mods,FUN=function(x) x$DIC))
  DICs.rad_outs$deltaDIC<- DICs.rad_outs$DIC-min(DICs.rad_outs$DIC)
  DICs.rad_outs$DICweight <- with(data.frame(DICs.rad_outs), exp(-deltaDIC/2) / sum(exp(-deltaDIC/2)))
  DICs.rad_outs[,c('R2.marginal','R2.conditional')]<-do.call(rbind,lapply(mods,FUN=get.R2.MCMCglmm))
  DICs.rad_outs<-DICs.rad_outs[order(DICs.rad_outs$deltaDIC),]
  
  # Model coefficient estimates
  akk<-function(x){
    sm<-summary(x)
    out<-as.data.frame(sm$solutions[,c('post.mean','pMCMC')])
    out$cf<-row.names(out)
    return(out)
  }
  cf.res<-akk(mod.direct)
  cf.res$marker<-names(delts)[x]
  
  # Compile results
  # (Model level diagnostics, coefficient estimates)
  results_list=list(DICs.rad_outs, cf.res)
  
  return(results_list)	
}

## Test out and time this function:

# Target markers:
names(delts)[1:10]  # marker names start in column 9
xs<-seq(8,ncol(delts))  # these are the indices for all markers that need to be analyzed
length(xs)    # 2066 in total

# Try it out for the first two markers, and a very short MCMCglmm run:
bob<-lapply(c(9:10),markerglmm,nitt=300,burnin=50,thin=2)

# Process list returned by markerglmm and lapply:
bob2<-do.call(rbind,lapply(bob,acc,2))
bob2

# Time the execution of full length MCMCglmm analyses on two markers, using 2 cores:
# CAUTION - this takes several minutes.
#system.time(
#  friend<-mclapply(c(9:10), markerglmm, mc.cores=2,nitt=1000000, burnin=500000, thin=750)
#)

####################################################################################

####      b) Run complete set of MCMCglmm analyses                              ####

# Again, make sure the right set of indices, matching to all 2066 markers, is targeted:
#names(delts)[1:10]
xs<-seq(8,ncol(delts))
length(xs)

# RUN Analysis in parallel - 
# CAUTION: this will take ~13 hours and occupy 6 cores of a 2.8 GHz Intel Core i7 CPU. 
markers_friend=mclapply(xs, markerglmm, mc.cores=6,nitt=1000000, burnin=500000, 
                           thin=750,mc.preschedule = F)

#save(markers_friend,file="all_markers_all_direct.RData")



####      c) Process output of MCMCglmm runs                              ####

#load(file="all_markers_all_direct.RData")


#  Extract model level output:
bob1<-do.call(rbind,lapply(markers_friend,acc,1))

# NOTE: for two markers, MCMCglmm fails to run - TP88062 and TP41983
# These are both 'low outliers' based on their FST, and have no variation between populations.
# They are excluded from subsequent analyses.

# [1] "TP88062"
# 1354: "Error in MCMCglmm(fixed = fixed, random = ~idv(mult.memb(~pop1 + pop2)),  : \n  Mixed model equations singular: use a (stronger) prior\n\n"

# [1] "TP41983"
# [811] "Error in MCMCglmm(fixed = fixed, random = ~idv(mult.memb(~pop1 + pop2)),  : \n  Mixed model equations singular: use a (stronger) prior\n\n"


#  Extract model coefficients, accounting for 2 bad runs:
bob2<-do.call(rbind,lapply(markers_friend,acc,2))
bob2<-na.omit(bob2)   # omitting results of TP88062 and TP41983
head(bob2)


####      d) Combine model coefficients with marker data                        ####

# pull in information from fst.tab, loaded at top of document:
cf.data<-merge(bob2,fst.tab,by='marker',all.x=T)
head(cf.data)



####      e) Focus in on a single marker for additional graphs                  ####

cf.data[cf.data$marker=="TP26245",]

# load input data
delt<-read.csv("pairwise_population_differences.csv")

# pull data for focal marker
loc<-which(names(delt)=="TP26245")
tmp.data<-delt[,c(1:7,loc)]
head(tmp.data)

# run model:
tmark<-MCMCglmm(fixed=TP26245~goldblattsame+physdist+distpc1.a+
                    distpc2.a+distpc3.a,
                random= ~ idv(mult.memb(~ pop1+pop2)), data=tmp.data, 
                pr=TRUE, nitt=1000000, burnin=500000, thin=750, verbose=FALSE)
summary(tmark)


# pull out points for partial regression plots:
obs<-melt(tmp.data,id.vars=c('pop1','pop2','TP26245'))
head(obs)

# calculate mean values of each covariate:
means <- obs %>% group_by(variable) %>% summarise(mn=mean(value))
means$variable<-as.character(means$variable)
means<-rbind(means,c('(Intercept)',1))
means$mn<-as.numeric(means$mn)
names(means)[1]<-'cf'

# extract and format coefficients
cfs<-summary(tmark)$solutions[,c('post.mean')]
cfs<-melt(cfs)
cfs$cf<-rownames(cfs)

# merge means and coefficients
ests<-merge(cfs,means)

# calculate mean predictions for partial regression plots:
ests$term<-ests$value*ests$mn

# now use this as the basis for making predictions
uvars<-unique(obs$variable)
i<-1
full.res<-NULL
for(i in 1:5){
  uvars[i]
  xs<-obs$value[obs$variable==uvars[i]]
  
  loc<-which(ests$cf==uvars[i])
  pred.xs<-ests$value[loc]*xs
  
  int<-rep(ests$term[1],171)
  pc1<-rep(ests$term[2],171)
  pc2<-rep(ests$term[3],171)
  pc3<-rep(ests$term[4],171)
  gbt<-rep(ests$term[5],171)
  phys<-rep(ests$term[6],171)
  
  contribs<-data.frame(int,pc1,pc2,pc3,gbt,phys)
  contribs[,loc]<-pred.xs
  
  ys<-apply(contribs,1,sum)  
  
  res<-data.frame(variable=uvars[i],value=xs,pred=ys)
  
  full.res<-rbind(full.res,res)    
}
head(full.res)


p1<-ggplot(obs,aes(x=value))+
  geom_point(aes(y=TP26245),alpha=0.3)+
  geom_line(data=full.res,aes(y=pred),colour='red')+
  facet_wrap(~variable)+
  theme_bw()+
  theme(panel.grid = element_blank())

ggsave("mcmcglmm_fit_TP26245.pdf",p1,width=6.5,height=5)



#################################################################################
#### (4) randomize allele frequencies for 5 null models                      ####
#################################################################################

# Depends on afreq data set of allele frequencies, loaded at top of document

# Shuffle allele frequencies 5 times: 
# - and save output as new afreq files
for(p in 1:5){
  afreq.R <- afreq
  inds<-seq(1,nrow(afreq.R))
  
  for(i in 2:ncol(afreq.R)){
    shuffled.inds<-sample(inds,size = length(inds),replace = F)
    afreq.R[,i]<-afreq.R[shuffled.inds,i]
  }
  
  # re-calculate delta for each marker
  m<-ncol(afreq)-1  # of markers
  n<-nrow(afreq)  # of populations
  npairs<-sum(seq(1,n-1,1))  # number of unique pairings of 19 populations
  delt3<-matrix(nrow=npairs, ncol=m)  # storage for data
  pop.mat<-matrix(nrow=npairs,ncol=2)   # store population ID
  for(d in 1:m){
    k<-0
    for(i in 1:(n-1)){
      for(j in (i+1):n){	
        k<-k+1
        delt3[k,d]<-delta(c(afreq.R[i,d+1], (1-afreq.R[i,d+1])), c(afreq.R[j,d+1], (1-afreq.R[j,d+1])))
        
        if(d==1){
          pop.mat[k,]<-c(afreq.R[i,1],afreq.R[j,1])
        }
      }
    }
  }
  # process data
  pop.mat<-as.data.frame(pop.mat)
  delt3<-as.data.frame(delt3)
  delt3<-data.frame(pop1=pop.mat[,1],pop2=pop.mat[,2],delt3)
  names(delt3)[3:ncol(delt3)]<-names(afreq.R)[2:ncol(afreq.R)]

  delt3<-data.frame(pop.pair=paste(delt3$pop1,delt3$pop2),delt3)
  
  # merge with environmental data
  #head(env.dists)    # comes from part 2b
  delt3<-merge(env.dists,delt3,all.y=T,by='pop.pair')
  delt3<-delt3[,-which(names(delt3)=='pop.pair')] # remove organizing column
  nms<-c('pop1','pop2','physdist','distpc1.a','distpc2.a','distpc3.a','goldblattsame',names(delt3)[8:ncol(delt3)])
  delt3<-delt3[,nms]  # re-arrange columns
  
  # save output:
  write.csv(delt3,paste("randomized_all_marker_deltas_for_analysis_v",p,".csv",sep=""),row.names=F)
}




#################################################################################
#### (5) run MCMCglmm analyses on null models                                ####
#################################################################################


####  a) Set up and run simulation #1                                        ####

#   - for details on code, see above section (3)

# Read individual marker distances based on randomized allele frequencies
delts<-read.csv("randomized_all_marker_deltas_for_analysis_v1.csv")

# Re-initialize function for running model sets, with new values for delts:
markerglmm=function(x, nitt, burnin, thin){
  
  # Wrapper function for MCMCglmm run
  MCMCglmm2<-function(fixed, data=data){
    test<-MCMCglmm(fixed=fixed, random= ~ idv(mult.memb(~ pop1+pop2)), data=data, 
                   pr=TRUE, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
    return(test)
  }
  
  # Pull out the name of marker/locus of focus, given x...
  nm<-names(delts)[x]
  
  # Fit model, invoking main effects of all covariates
  mod.direct<-eval(substitute(
    MCMCglmm2(fixed=tmp~goldblattsame+physdist+distpc1.a+distpc2.a+distpc3.a,data=delts),
    list(tmp=as.name(nm))))
  
  # list of models
  mods<-list(mod.direct)
  
  # Combine model results:
  # this is legacy code from running multiple models/conducting model comparison; the more important
  # section for current analyses follows under 'Model coefficient estimates'
  DICs.rad_outs <- data.frame(array(NA, c(1, 7)))
  names(DICs.rad_outs) <- c('model','outlier','DIC', 'deltaDIC', 'DICweight',
                            'R2.marginal','R2.conditional')
  DICs.rad_outs$model <- c('alldirect')
  DICs.rad_outs$outlier<-names(delts)[x]
  DICs.rad_outs$DIC<-unlist(lapply(mods,FUN=function(x) x$DIC))
  DICs.rad_outs$deltaDIC<- DICs.rad_outs$DIC-min(DICs.rad_outs$DIC)
  DICs.rad_outs$DICweight <- with(data.frame(DICs.rad_outs), exp(-deltaDIC/2) / sum(exp(-deltaDIC/2)))
  DICs.rad_outs[,c('R2.marginal','R2.conditional')]<-do.call(rbind,lapply(mods,FUN=get.R2.MCMCglmm))
  DICs.rad_outs<-DICs.rad_outs[order(DICs.rad_outs$deltaDIC),]
  
  # Model coefficient estimates
  akk<-function(x){
    sm<-summary(x)
    out<-as.data.frame(sm$solutions[,c('post.mean','pMCMC')])
    out$cf<-row.names(out)
    return(out)
  }
  cf.res<-akk(mod.direct)
  cf.res$marker<-names(delts)[x]
  
  # Compile results
  # (Model level diagnostics, coefficient estimates)
  results_list=list(DICs.rad_outs, cf.res)
  
  return(results_list)	
}

# Testing:
#friend=mclapply(c(9:12), markerglmm, mc.cores=2,nitt=100, burnin=5, thin=10)
#bob2<-do.call(rbind,lapply(friend,acc,2))
#bob2<-na.omit(bob2)
#head(bob2)

#system.time(friend<-mclapply(c(9:10), markerglmm, mc.cores=2,nitt=1000000, burnin=500000, thin=750))

#### Run MCMCglmms

# column indicies for markers:
xs<-seq(8,ncol(delts))
length(xs)

# CAUTION - 13+ hour run time:
markers_randomized_v1=mclapply(xs, markerglmm, mc.cores=6,nitt=1000000, burnin=500000, 
                                  thin=750,mc.preschedule = F)

#save(markers_randomized_v1,file="randomized_all_markers_all_direct_v1.RData")
#load(file="randomized_all_markers_all_direct_v1.RData")

#  extract sim1 model coefficients
b1<-do.call(rbind,lapply(markers_randomized_v1,acc,2))
b1<-na.omit(b1)
row.names(b1)<-NULL

# pull in information from fst.tab, loaded at top of document:
sim1.data<-merge(b1,fst.tab,by='marker',all.x=T)
sim1.data<-data.frame(run='sim1',sim1.data)
head(sim1.data)




####  b) Set up and run simulation #2                                        ####

#   - for details on code, see above section (3)

# Read individual marker distances based on randomized allele frequencies
delts<-read.csv("randomized_all_marker_deltas_for_analysis_v2.csv")

# Re-initialize function for running model sets, with new values for delts:
markerglmm=function(x, nitt, burnin, thin){
  
  # Wrapper function for MCMCglmm run
  MCMCglmm2<-function(fixed, data=data){
    test<-MCMCglmm(fixed=fixed, random= ~ idv(mult.memb(~ pop1+pop2)), data=data, 
                   pr=TRUE, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
    return(test)
  }
  
  # Pull out the name of marker/locus of focus, given x...
  nm<-names(delts)[x]
  
  # Fit model, invoking main effects of all covariates
  mod.direct<-eval(substitute(
    MCMCglmm2(fixed=tmp~goldblattsame+physdist+distpc1.a+distpc2.a+distpc3.a,data=delts),
    list(tmp=as.name(nm))))
  
  # list of models
  mods<-list(mod.direct)
  
  # Combine model results:
  # this is legacy code from running multiple models/conducting model comparison; the more important
  # section for current analyses follows under 'Model coefficient estimates'
  DICs.rad_outs <- data.frame(array(NA, c(1, 7)))
  names(DICs.rad_outs) <- c('model','outlier','DIC', 'deltaDIC', 'DICweight',
                            'R2.marginal','R2.conditional')
  DICs.rad_outs$model <- c('alldirect')
  DICs.rad_outs$outlier<-names(delts)[x]
  DICs.rad_outs$DIC<-unlist(lapply(mods,FUN=function(x) x$DIC))
  DICs.rad_outs$deltaDIC<- DICs.rad_outs$DIC-min(DICs.rad_outs$DIC)
  DICs.rad_outs$DICweight <- with(data.frame(DICs.rad_outs), exp(-deltaDIC/2) / sum(exp(-deltaDIC/2)))
  DICs.rad_outs[,c('R2.marginal','R2.conditional')]<-do.call(rbind,lapply(mods,FUN=get.R2.MCMCglmm))
  DICs.rad_outs<-DICs.rad_outs[order(DICs.rad_outs$deltaDIC),]
  
  # Model coefficient estimates
  akk<-function(x){
    sm<-summary(x)
    out<-as.data.frame(sm$solutions[,c('post.mean','pMCMC')])
    out$cf<-row.names(out)
    return(out)
  }
  cf.res<-akk(mod.direct)
  cf.res$marker<-names(delts)[x]
  
  # Compile results
  # (Model level diagnostics, coefficient estimates)
  results_list=list(DICs.rad_outs, cf.res)
  
  return(results_list)	
}

# Testing:
#friend=mclapply(c(9:12), markerglmm, mc.cores=2,nitt=100, burnin=5, thin=10)
#bob2<-do.call(rbind,lapply(friend,acc,2))
#bob2<-na.omit(bob2)
#head(bob2)

#system.time(friend<-mclapply(c(9:10), markerglmm, mc.cores=2,nitt=1000000, burnin=500000, thin=750))

#### Run MCMCglmms

# column indicies for markers:
xs<-seq(8,ncol(delts))
length(xs)

# CAUTION - 13+ hour run time:
markers_randomized_v2=mclapply(xs, markerglmm, mc.cores=6,nitt=1000000, burnin=500000, 
                               thin=750,mc.preschedule = F)

#save(markers_randomized_v2,file="randomized_all_markers_all_direct_v1.RData")
#load(file="randomized_all_markers_all_direct_v2.RData")

#  extract sim1 model coefficients
b1<-do.call(rbind,lapply(markers_randomized_v2,acc,2))
b1<-na.omit(b1)
row.names(b1)<-NULL

# pull in information from fst.tab, loaded at top of document:
sim2.data<-merge(b1,fst.tab,by='marker',all.x=T)
sim2.data<-data.frame(run='sim2',sim2.data)
head(sim2.data)




####  c) Set up and run simulation #3                                        ####

#   - for details on code, see above section (3)

# Read individual marker distances based on randomized allele frequencies
delts<-read.csv("randomized_all_marker_deltas_for_analysis_v3.csv")

# Re-initialize function for running model sets, with new values for delts:
markerglmm=function(x, nitt, burnin, thin){
  
  # Wrapper function for MCMCglmm run
  MCMCglmm2<-function(fixed, data=data){
    test<-MCMCglmm(fixed=fixed, random= ~ idv(mult.memb(~ pop1+pop2)), data=data, 
                   pr=TRUE, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
    return(test)
  }
  
  # Pull out the name of marker/locus of focus, given x...
  nm<-names(delts)[x]
  
  # Fit model, invoking main effects of all covariates
  mod.direct<-eval(substitute(
    MCMCglmm2(fixed=tmp~goldblattsame+physdist+distpc1.a+distpc2.a+distpc3.a,data=delts),
    list(tmp=as.name(nm))))
  
  # list of models
  mods<-list(mod.direct)
  
  # Combine model results:
  # this is legacy code from running multiple models/conducting model comparison; the more important
  # section for current analyses follows under 'Model coefficient estimates'
  DICs.rad_outs <- data.frame(array(NA, c(1, 7)))
  names(DICs.rad_outs) <- c('model','outlier','DIC', 'deltaDIC', 'DICweight',
                            'R2.marginal','R2.conditional')
  DICs.rad_outs$model <- c('alldirect')
  DICs.rad_outs$outlier<-names(delts)[x]
  DICs.rad_outs$DIC<-unlist(lapply(mods,FUN=function(x) x$DIC))
  DICs.rad_outs$deltaDIC<- DICs.rad_outs$DIC-min(DICs.rad_outs$DIC)
  DICs.rad_outs$DICweight <- with(data.frame(DICs.rad_outs), exp(-deltaDIC/2) / sum(exp(-deltaDIC/2)))
  DICs.rad_outs[,c('R2.marginal','R2.conditional')]<-do.call(rbind,lapply(mods,FUN=get.R2.MCMCglmm))
  DICs.rad_outs<-DICs.rad_outs[order(DICs.rad_outs$deltaDIC),]
  
  # Model coefficient estimates
  akk<-function(x){
    sm<-summary(x)
    out<-as.data.frame(sm$solutions[,c('post.mean','pMCMC')])
    out$cf<-row.names(out)
    return(out)
  }
  cf.res<-akk(mod.direct)
  cf.res$marker<-names(delts)[x]
  
  # Compile results
  # (Model level diagnostics, coefficient estimates)
  results_list=list(DICs.rad_outs, cf.res)
  
  return(results_list)	
}

# Testing:
#friend=mclapply(c(9:12), markerglmm, mc.cores=2,nitt=100, burnin=5, thin=10)
#bob2<-do.call(rbind,lapply(friend,acc,2))
#bob2<-na.omit(bob2)
#head(bob2)

#system.time(friend<-mclapply(c(9:10), markerglmm, mc.cores=2,nitt=1000000, burnin=500000, thin=750))

#### Run MCMCglmms

# column indicies for markers:
xs<-seq(8,ncol(delts))
length(xs)

# CAUTION - 13+ hour run time:
markers_randomized_v3=mclapply(xs, markerglmm, mc.cores=6,nitt=1000000, burnin=500000, 
                               thin=750,mc.preschedule = F)

#save(markers_randomized_v3,file="randomized_all_markers_all_direct_v1.RData")
#load(file="randomized_all_markers_all_direct_v3.RData")

#  extract sim1 model coefficients
b1<-do.call(rbind,lapply(markers_randomized_v3,acc,2))
b1<-na.omit(b1)
row.names(b1)<-NULL

# pull in information from fst.tab, loaded at top of document:
sim3.data<-merge(b1,fst.tab,by='marker',all.x=T)
sim3.data<-data.frame(run='sim3',sim3.data)
head(sim3.data)




####  d) Set up and run simulation #4                                        ####

#   - for details on code, see above section (3)

# Read individual marker distances based on randomized allele frequencies
delts<-read.csv("randomized_all_marker_deltas_for_analysis_v4.csv")

# Re-initialize function for running model sets, with new values for delts:
markerglmm=function(x, nitt, burnin, thin){
  
  # Wrapper function for MCMCglmm run
  MCMCglmm2<-function(fixed, data=data){
    test<-MCMCglmm(fixed=fixed, random= ~ idv(mult.memb(~ pop1+pop2)), data=data, 
                   pr=TRUE, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
    return(test)
  }
  
  # Pull out the name of marker/locus of focus, given x...
  nm<-names(delts)[x]
  
  # Fit model, invoking main effects of all covariates
  mod.direct<-eval(substitute(
    MCMCglmm2(fixed=tmp~goldblattsame+physdist+distpc1.a+distpc2.a+distpc3.a,data=delts),
    list(tmp=as.name(nm))))
  
  # list of models
  mods<-list(mod.direct)
  
  # Combine model results:
  # this is legacy code from running multiple models/conducting model comparison; the more important
  # section for current analyses follows under 'Model coefficient estimates'
  DICs.rad_outs <- data.frame(array(NA, c(1, 7)))
  names(DICs.rad_outs) <- c('model','outlier','DIC', 'deltaDIC', 'DICweight',
                            'R2.marginal','R2.conditional')
  DICs.rad_outs$model <- c('alldirect')
  DICs.rad_outs$outlier<-names(delts)[x]
  DICs.rad_outs$DIC<-unlist(lapply(mods,FUN=function(x) x$DIC))
  DICs.rad_outs$deltaDIC<- DICs.rad_outs$DIC-min(DICs.rad_outs$DIC)
  DICs.rad_outs$DICweight <- with(data.frame(DICs.rad_outs), exp(-deltaDIC/2) / sum(exp(-deltaDIC/2)))
  DICs.rad_outs[,c('R2.marginal','R2.conditional')]<-do.call(rbind,lapply(mods,FUN=get.R2.MCMCglmm))
  DICs.rad_outs<-DICs.rad_outs[order(DICs.rad_outs$deltaDIC),]
  
  # Model coefficient estimates
  akk<-function(x){
    sm<-summary(x)
    out<-as.data.frame(sm$solutions[,c('post.mean','pMCMC')])
    out$cf<-row.names(out)
    return(out)
  }
  cf.res<-akk(mod.direct)
  cf.res$marker<-names(delts)[x]
  
  # Compile results
  # (Model level diagnostics, coefficient estimates)
  results_list=list(DICs.rad_outs, cf.res)
  
  return(results_list)	
}

# Testing:
#friend=mclapply(c(9:12), markerglmm, mc.cores=2,nitt=100, burnin=5, thin=10)
#bob2<-do.call(rbind,lapply(friend,acc,2))
#bob2<-na.omit(bob2)
#head(bob2)

#system.time(friend<-mclapply(c(9:10), markerglmm, mc.cores=2,nitt=1000000, burnin=500000, thin=750))

#### Run MCMCglmms

# column indicies for markers:
xs<-seq(8,ncol(delts))
length(xs)

# CAUTION - 13+ hour run time:
markers_randomized_v4=mclapply(xs, markerglmm, mc.cores=6,nitt=1000000, burnin=500000, 
                               thin=750,mc.preschedule = F)

#save(markers_randomized_v4,file="randomized_all_markers_all_direct_v1.RData")
#load(file="randomized_all_markers_all_direct_v4.RData")

#  extract sim1 model coefficients
b1<-do.call(rbind,lapply(markers_randomized_v4,acc,2))
b1<-na.omit(b1)
row.names(b1)<-NULL

# pull in information from fst.tab, loaded at top of document:
sim4.data<-merge(b1,fst.tab,by='marker',all.x=T)
sim4.data<-data.frame(run='sim4',sim4.data)
head(sim4.data)




####  e) Set up and run simulation #5                                        ####

#   - for details on code, see above section (3)

# Read individual marker distances based on randomized allele frequencies
delts<-read.csv("randomized_all_marker_deltas_for_analysis_v5.csv")

# Re-initialize function for running model sets, with new values for delts:
markerglmm=function(x, nitt, burnin, thin){
  
  # Wrapper function for MCMCglmm run
  MCMCglmm2<-function(fixed, data=data){
    test<-MCMCglmm(fixed=fixed, random= ~ idv(mult.memb(~ pop1+pop2)), data=data, 
                   pr=TRUE, nitt=nitt, burnin=burnin, thin=thin, verbose=FALSE)
    return(test)
  }
  
  # Pull out the name of marker/locus of focus, given x...
  nm<-names(delts)[x]
  
  # Fit model, invoking main effects of all covariates
  mod.direct<-eval(substitute(
    MCMCglmm2(fixed=tmp~goldblattsame+physdist+distpc1.a+distpc2.a+distpc3.a,data=delts),
    list(tmp=as.name(nm))))
  
  # list of models
  mods<-list(mod.direct)
  
  # Combine model results:
  # this is legacy code from running multiple models/conducting model comparison; the more important
  # section for current analyses follows under 'Model coefficient estimates'
  DICs.rad_outs <- data.frame(array(NA, c(1, 7)))
  names(DICs.rad_outs) <- c('model','outlier','DIC', 'deltaDIC', 'DICweight',
                            'R2.marginal','R2.conditional')
  DICs.rad_outs$model <- c('alldirect')
  DICs.rad_outs$outlier<-names(delts)[x]
  DICs.rad_outs$DIC<-unlist(lapply(mods,FUN=function(x) x$DIC))
  DICs.rad_outs$deltaDIC<- DICs.rad_outs$DIC-min(DICs.rad_outs$DIC)
  DICs.rad_outs$DICweight <- with(data.frame(DICs.rad_outs), exp(-deltaDIC/2) / sum(exp(-deltaDIC/2)))
  DICs.rad_outs[,c('R2.marginal','R2.conditional')]<-do.call(rbind,lapply(mods,FUN=get.R2.MCMCglmm))
  DICs.rad_outs<-DICs.rad_outs[order(DICs.rad_outs$deltaDIC),]
  
  # Model coefficient estimates
  akk<-function(x){
    sm<-summary(x)
    out<-as.data.frame(sm$solutions[,c('post.mean','pMCMC')])
    out$cf<-row.names(out)
    return(out)
  }
  cf.res<-akk(mod.direct)
  cf.res$marker<-names(delts)[x]
  
  # Compile results
  # (Model level diagnostics, coefficient estimates)
  results_list=list(DICs.rad_outs, cf.res)
  
  return(results_list)	
}

# Testing:
#friend=mclapply(c(9:12), markerglmm, mc.cores=2,nitt=100, burnin=5, thin=10)
#bob2<-do.call(rbind,lapply(friend,acc,2))
#bob2<-na.omit(bob2)
#head(bob2)

#system.time(friend<-mclapply(c(9:10), markerglmm, mc.cores=2,nitt=1000000, burnin=500000, thin=750))

#### Run MCMCglmms

# column indicies for markers:
xs<-seq(8,ncol(delts))
length(xs)

# CAUTION - 13+ hour run time:
markers_randomized_v5=mclapply(xs, markerglmm, mc.cores=6,nitt=1000000, burnin=500000, 
                               thin=750,mc.preschedule = F)

#save(markers_randomized_v5,file="randomized_all_markers_all_direct_v1.RData")
#load(file="randomized_all_markers_all_direct_v5.RData")

#  extract sim1 model coefficients
b1<-do.call(rbind,lapply(markers_randomized_v5,acc,2))
b1<-na.omit(b1)
row.names(b1)<-NULL

# pull in information from fst.tab, loaded at top of document:
sim5.data<-merge(b1,fst.tab,by='marker',all.x=T)
sim5.data<-data.frame(run='sim5',sim5.data)
head(sim5.data)



####    f) Combine simulation results with original data:   #####

cf.data<-data.frame(run='original',cf.data)

mod.results<-rbind(cf.data,sim1.data,sim2.data,sim3.data,sim4.data,sim5.data)

# save output:
write.csv(mod.results,"all_coefficient_estimates_combined.csv",row.names=F)




#################################################################################
#### (6) fit GAM models and plot results                                     ####
#################################################################################


####    a) run GAM fits                                                     #####

# load MCMCglmm model results
cfs<-read.csv("all_coefficient_estimates_combined.csv")
head(cfs)

# remove non-outliers and observation with extremely high fst:
cfs2<-cfs[cfs$fst<=0.29 & cfs$level!='low', ]

# change the sign of the 'goldblattsame' coefficient; this simply flips the interpretation from examining the effect of being in the same phytogeographic province to being in different provinces.
cfs2$post.mean<-ifelse(cfs2$cf=='goldblattsame',-1*cfs2$post.mean,cfs2$post.mean)


# run GAM fits:
cfid<-unique(cfs2$cf)
runid<-unique(cfs2$run)
res<-c()
stat.mat<-c()
for(i in 1:length(runid)){
  for(j in 1:length(cfid)){
    tmp<-cfs2[cfs2$cf==cfid[j] & cfs2$run==runid[i],]
    gm1<-gam(post.mean~s(fst,k=3),sp=c(-1),data=tmp)
    
    smg<-summary(gm1)
    
    pd1<-predict(gm1,se.fit = T)
    tmp$pred<-pd1$fit
    tmp$se<-pd1$se.fit
    
    # save statistical outputs
    stat.mat<-rbind(stat.mat,c(runid[i],as.character(cfid[j]),smg$p.table[,c(1:4)],smg$edf,smg$s.table[,3:4],smg$n,smg$r.sq,smg$dev.expl))
    
    # GAM fit diagnostics:
    # gam.check(gm1)
    # plot(fitted(gm1),residuals(gm1))
    # qq.gam(gm1,rep=100)
    # plot(residuals(gm1)~tmp$fst)
    # plot(gm1)
    
    res<-rbind(res,tmp)    
  }
}
stat.mat<-as.data.frame(stat.mat)

# save archive version of the results:
resX<-res


####    b) plot GAM fits                                                     #####

# drop intercept term - exists on a much larger scale and is not the focus of the analysis
res<-res[res$cf!="(Intercept)",]

# 95% Confidence band
res$ci.se<- 1.96*res$se

# Change panel labeling options:
res$cf2<-factor(res$cf,levels=c('physdist','distpc1.a','distpc2.a','distpc3.a','goldblattsame'),labels=c('Distance','PC1','PC2','PC3','Phytogeography'))

res.sim<-res[res$run!='original',]
res.orig<-res[res$run=='original',]

dummy<-data.frame(fst=0.015,pred=0.065,txt=c('A.','B.','C.','D.','E.'),cf2=c('Distance','PC1','PC2','PC3','Phytogeography'))

# Comparison plot
g9<-ggplot(res.sim,aes(x=fst,y=pred,cf2))+
  geom_ribbon(aes(group=run,ymin=pred-ci.se,ymax=pred+ci.se),
              colour=NA,fill='gray',alpha=0.4)+
  geom_line(aes(group=run),colour=gray(0.7))+
  geom_ribbon(data=res.orig,aes(ymin=pred-ci.se,ymax=pred+ci.se),
              colour=NA,fill='blue',alpha=0.4)+
  geom_line(data=res.orig,colour='blue')+
  geom_hline(yintercept = 0, colour=gray(0.1), size=0.2)+
  geom_text(data=dummy,aes(label=txt))+
  facet_wrap(~cf2)+
  scale_x_continuous("Marker Fst",limits=c(0,0.25))+
  scale_y_continuous("Effect size")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_rect(colour='white',fill='white'),
        strip.text = element_blank())
g9

ggsave("GAM_fits_comparison_FigXXX_v4.pdf",g9,width=6.5,height=5)




####    c) GAM fit against original data                                     #####

# confidence band:
resX$ci.se<- 1.96*resX$se
resX$pred<-as.numeric(resX$pred)

# various panel labelling options:
resX$cf2<-factor(resX$cf,levels=c('physdist','distpc1.a','distpc2.a','distpc3.a','goldblattsame','(Intercept)'),labels=c('Distance','PC1','PC2','PC3','Phytogeography','Intercept'))

dummy2<-data.frame(fst=0.015,pred=c(rep(0.29,5),0.58),txt=c('A.','B.','C.','D.','E.','F.'),cf2=c('Distance','PC1','PC2','PC3','Phytogeography','Intercept'),cf=c("physdist","distpc1.a","distpc2.a","distpc3.a","goldblattsame","(Intercept)"))

cfid<-c("physdist","distpc1.a","distpc2.a","distpc3.a","goldblattsame","(Intercept)")
plts<-NULL
for(i in 1:length(cfid)){
  if(cfid[i]=="(Intercept)"){
    ylims<-c(0,0.6)
  }else{
    ylims<-c(-0.31,0.31)
  }
  
  dt<-resX[resX$run=="original" & resX$fst<0.25 & resX$cf==cfid[i],]
  plts[[i]]<-ggplot(data=dt,aes(x=fst,y=pred,cf))+
    geom_point(aes(y=post.mean,colour=level),alpha=0.2)+
    geom_ribbon(aes(ymin=pred-ci.se,ymax=pred+ci.se),
                colour=NA,fill='gray',alpha=0.7)+
    geom_line(aes(y=pred))+
    geom_hline(yintercept=0,colour=gray(0.1), size=0.2)+
    geom_text(data=dummy2[dummy2$cf==cfid[i],],aes(label=txt))+
    scale_x_continuous("Marker Fst",limits=c(0,0.25))+
    scale_y_continuous("Effect size",limits=ylims)+
    scale_color_discrete(guide=F)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_rect(colour='white',fill='white'),
          strip.text = element_blank())
}

# Show together:
grid.arrange(plts[[1]],plts[[2]],plts[[3]],plts[[4]],plts[[5]],plts[[6]],nrow=2)

# format for export
ag1<-arrangeGrob(plts[[1]],plts[[2]],plts[[3]],plts[[4]],plts[[5]],plts[[6]],nrow=2)

ggsave("GAM_fit_w_points_original_FigXXX_v3.pdf",ag1,width=6.5,height=5)




####    d) GAM fit against simulation1 data                                     #####

# Simulation:

# confidence band:
resX$ci.se<- 1.96*resX$se
resX$pred<-as.numeric(resX$pred)

# various panel labelling options:
resX$cf2<-factor(resX$cf,levels=c('physdist','distpc1.a','distpc2.a','distpc3.a','goldblattsame','(Intercept)'),labels=c('Distance','PC1','PC2','PC3','Phytogeography','Intercept'))

dummy2<-data.frame(fst=0.015,pred=c(rep(0.29,5),0.58),txt=c('A.','B.','C.','D.','E.','F.'),cf2=c('Distance','PC1','PC2','PC3','Phytogeography','Intercept'),cf=c("physdist","distpc1.a","distpc2.a","distpc3.a","goldblattsame","(Intercept)"))

cfid<-c("physdist","distpc1.a","distpc2.a","distpc3.a","goldblattsame","(Intercept)")
plts2<-NULL
for(i in 1:length(cfid)){
  if(cfid[i]=="(Intercept)"){
    ylims<-c(0,0.6)
  }else{
    ylims<-c(-0.31,0.31)
  }
  
  dt<-resX[resX$run=="sim1" & resX$fst<0.25 & resX$cf==cfid[i],]
  plts2[[i]]<-ggplot(data=dt,aes(x=fst,y=pred,cf))+
    geom_point(aes(y=post.mean,colour=level),alpha=0.2)+
    geom_ribbon(aes(ymin=pred-ci.se,ymax=pred+ci.se),
                colour=NA,fill='gray',alpha=0.7)+
    geom_line(aes(y=pred))+
    geom_hline(yintercept=0,colour=gray(0.1), size=0.2)+
    geom_text(data=dummy2[dummy2$cf==cfid[i],],aes(label=txt))+
    scale_x_continuous("Marker Fst",limits=c(0,0.25))+
    scale_y_continuous("Effect size",limits=ylims)+
    scale_color_discrete(guide=F)+
    theme_bw()+
    theme(panel.grid = element_blank(),
          strip.background = element_rect(colour='white',fill='white'),
          strip.text = element_blank())
}

# Show together:
grid.arrange(plts2[[1]],plts2[[2]],plts2[[3]],plts2[[4]],plts2[[5]],plts2[[6]],nrow=2)

ag2<-arrangeGrob(plts2[[1]],plts2[[2]],plts2[[3]],plts2[[4]],plts2[[5]],plts2[[6]],nrow=2)

ggsave("GAM_fit_w_points_sim1_FigXXX_v3.pdf",ag2,width=6.5,height=5)



####    e) Show distribution of high outliers & non-outliers vs FST          #####

# Extract subset of data for plotting
tmp<-cfs[cfs$run=='original' & cfs$cf=='physdist' & cfs$fst <0.29,]

# Pull up ggplot default colours
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Show plot
g6<-ggplot(tmp,aes(level2))+
  geom_hline(yintercept=0,colour=gray(0.1), size=0.2)+
  geom_freqpoly(aes(fst,colour=level),breaks=seq(0.005,0.25,0.005))+
  scale_x_continuous('Marker Fst',expand=c(0,0),limits=c(0,0.26))+
  scale_y_continuous('Frequency',expand=c(0,0),limits=c(-1,151))+
  scale_color_manual('Marker type',
                     values=c(gg_color_hue(2)[1],'purple',gg_color_hue(2)[2]))+
  theme_bw()+
  theme(panel.grid = element_blank())
g6  

# Save it:
ggsave("outlier_status_frequency_plot_v3.pdf",g6,width=5,height=4.1)




####    f) Save statistical output from GAM fits                            #####

# Organize statistical output
for(j in 3:ncol(stat.mat)){
  stat.mat[,j]<-as.numeric(as.character(stat.mat[,j]))
}
names(stat.mat)[1:2]<-c("run","cf")
names(stat.mat)[7]<-"edf"
names(stat.mat)[10:12]<-c('n','r2','dev.expl')
head(stat.mat)

# save output:
write.csv(stat.mat,"gam_fit_stats_output.csv",row.names=F)


