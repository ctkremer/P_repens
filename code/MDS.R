######################################

#### This code runs the MDS analyses ####

# Written by K. Holsinger

######################################

## Load required tools ##
require(ggplot2)
require(maps)
require(mapdata)
require(plyr)
require(labeling)


#### 1. Set up functions for MDS analysis ####

rm(list=ls())

strip <- function(x) {
  y <- vector(mode="character", length=length(x))
  for (i in 1:length(x)) {
    y[i] <- strsplit(x[i], "_")[[1]][1]
  }
  y
}

read.marker.data <- function(filename) {
  markers <- read.csv(filename, header=TRUE, na.strings=".")

  markers$pop <- as.factor(strip(as.character(markers$names)))
  markers <- subset(markers, pop!="EMPTY")
  markers$names <- NULL
  n.pops <-length(unique(markers$pop))
  n.loci <- ncol(markers)-1
  list(markers=markers, n.pops=n.pops, n.loci=n.loci, n.indivs=nrow(markers))
}

get.freqs <- function(x) {
  n.loci <- ncol(x)
  p <- numeric(n.loci)
  for (i in 1:n.loci) {
    k <- sum(x[,i], na.rm=TRUE)
    n <- nrow(x) - length(x[is.na(x[,2]),2])
    p[i] <- k/(2*n)
  }
  p
}

get.ibd.dist <- function(x, y) {
  n.loci <- length(x)
  f <- 0.0
  k <- 0
  for (i in 1:n.loci) {
    if (!is.na(x[i]) & !is.na(y[i])) {
      p.1 <- x[i]
      q.1 <- 1.0 - p.1
      p.2 <- y[i]
      q.2 <- 1.0 - p.2
      f <- f + (p.1*p.2 + q.1*q.2)
      k <- k + 1
    }
  }
  f <- f/k
  f <- (1-f)/f
  f
}

set.colors <- function(x, min.alpha=0.1, max.alpha=0.9) {
  min.x <- min(x)
  max.x <- max(x)
  alpha <- (max.alpha-min.alpha)*(x - min.x)/(max.x - min.x) + min.alpha
  color <- rgb(1.0, 0.0, 0.0, alpha)
  color
}

set.colors.rb <- function(x, y) {
  min.x <- min(x)
  max.x <- max(x)
  min.y <- min(y)
  max.y <- max(y)
  red <- 1.0 - (x - min.x)/(max.x - min.x)
  blue <- 1.0 - (y - min.y)/(max.y - min.y)
  color <- rgb(red, 0.0, blue)
  color
}

#### 2. performing MDS analysis on FST of non-outlier markers (figure 2) ####

## Load data
dat <- read.marker.data("./data/noouts_data.csv")

## Process it

# First individual allele frequencies at each locus
p <- get.freqs(dat$markers[,1:dat$n.loci])
p.indiv <- matrix(nrow=dat$n.indivs, ncol=dat$n.loci)
cat("Calculating individual allele frequencies at each locus...\n")
markers <- as.matrix(dat$markers)
for (i in 1:dat$n.indivs) {
  for (j in 1:dat$n.loci) {
    p.indiv[i,j] <- as.numeric(markers[i,j])/2.0
  }
}
#p.indiv[1:10,1:8]

# Calculate pairwise ibd distances
f.dist <- matrix(nrow=dat$n.indivs, ncol=dat$n.indivs)
cat("Calculating pairwise ibd distances...")
for (i in 1:(dat$n.indivs-1)) {
  f.dist[i,i] <- 0.0
  cat("\n from individual ", i, sep="")
  for (j in (i+1):dat$n.indivs) {
    if ((j %% 50) == 0) {
      cat(".", sep="")
    }
    f.dist[i,j] <- get.ibd.dist(p.indiv[i,],
                                p.indiv[j,])
    f.dist[j,i] <- f.dist[i,j]
  }
}
f.dist[dat$n.indivs,dat$n.indivs] <- 0.0
cat("\n")
#f.dist[1:10,1:8]


#### 3. Set up for exploring the results graphically ####

f.mds <- cmdscale(f.dist)

for.plot <- data.frame(pop=dat$markers$pop,
                       x=f.mds[,1],
                       y=f.mds[,2],
                       data="Individual trait")
for.plot$pop <- factor(for.plot$pop,
                       levels=c("ALC",
                                "LOE",
                                "KAR",
                                "BAV",
                                "UNI",
                                "SWA",
                                "GAR",
                                "KSW",
                                "ANY",
                                "POT",
                                "BRD",
                                "RND",
                                "MGU",
                                "KLM",
                                "CER",
                                "RIV",
                                "BAN",
                                "CDB",
                                "VAN"))
means <- ddply(for.plot, c("pop"), summarize, X=mean(x), Y=mean(y))
means$pop<- factor(means$pop,
                   levels=c("ALC",
                            "LOE",
                            "KAR",
                            "BAV",
                            "UNI",
                            "SWA",
                            "GAR",
                            "KSW",
                            "ANY",
                            "POT",
                            "BRD",
                            "RND",
                            "MGU",
                            "KLM",
                            "CER",
                            "RIV",
                            "BAN",
                            "CDB",
                            "VAN"))
means$colors.1 <- set.colors(means$X, min.alpha=0.2, max.alpha=1.0)
means$colors.2 <- set.colors(means$Y, min.alpha=0.2, max.alpha=1.0)
means$colors.3 <- set.colors.rb(means$X, means$Y)
means$data <- "Population mean"

colors <- numeric(length(means$pop))
for (i in 1:length(means$pop)) {
  colors[i] <- rgb(i-1, 0, 9, maxColorValue=length(means$pop)-1)
}




#### 4. Make figure 2,  green/yellow plot ####
dev.new()

p <- ggplot(indivcolors, aes(x=x, y=y, color= pop, size=data)) +
	theme_bw() +
     theme(panel.grid.major=element_blank())+
     theme(panel.grid.minor=element_blank())+
     theme(legend.key = element_blank())+
     geom_hline(yintercept=0, colour="gray")+
     geom_vline(xintercept=0, colour="gray")+
     geom_point(alpha=0.8) +
     geom_point(data=means1, aes(x=X, y=Y, color= pop, size=data)) +
     xlab("MDS f-dist Axis 1") +
     ylab("MDS f-dist Axis 2") +
     scale_color_manual(values=means1$newcol) +
     guides(fill=guide_legend(reverse=TRUE))
     
print(p)
ggsave("MDS-f-dist-repens_greenv2.pdf")




