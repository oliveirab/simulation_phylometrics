
################################################################
#### Simulate phylogenetic trees and calculate phylometrics
#### Brunno Oliveira, 2014
#### Universidade Federal do Rio Grande do Norte - Brasil
#### Stony Brook University - USA
################################################################

#### List of phylometrics

# Phylogenetic diversity (PD)  Sum of all branch lengths in a tree.	Faith (1992)
# Mean phylogenetic distance (MPD)	Mean of phylogenetic distances between all species pairs within a tree.	Webb (2000)
# Mean nearest taxon distance (MNTD)	Mean phylogenetic distance between each species and its nearest neighbor in a tree.	Webb (2000)
# Phylogenetic species variability (PSV)	Expected variance among species in a neutrally evolving trait.	Helmus et al. (2007) 
# Mean root distance (MRD)	The mean number of nodes that separating species from the root of their tree. Alternatively, when time trees are available MRD can be calculated as the mean distance between each species nodes and the root of its phylogenetic tree (MRD.time). 	Kerr and Currie (1999); Hawkins et al. (2012)
# Species evolutionary distinctiveness (ED)  After dividing each branch length by the number of species subtending that branch, the obtained values are summed across all branches from which a species descends.  	Redding and Mooers (2006); Isaac et al. (2007)
# Species-level diversification rate (DivRate)	Inverse of the evolutionary distinctness (ED). Species rapidly diversifying will have short edge lengths shared among many species.	Jetz et al. (2012), Kembel et al. (2010)
# Relative branch length (RBL)	Approximates time between speciation events across a tree by calculating the mean ((stem age – crown age) / crown age).	Marin and Blair, In prep.
# Diversification rates (DR)	Species accumulation through time (log(clade diversity)/clade age).	Magallon and Sanderson (2001); Nee (2001); Phillimore et al. (2006)
# Gamma statistics (GAM)	The distribution of branching events throughout the tree. This phylometric distinguished between trees with relatively long inter-nodal distances towards the tips (tippy trees) and trees with relatively longer inter-nodal distances towards the root of the phylogeny (stemmy trees). 	Pybus and Harvey (2000)
# Colless	Measures the degree of tree imbalance. Here, we calculated standardized values of Colless according to the Yule and PDA models.	Colless (1982)


#### load packages
library(geiger)
library(ape)
library(TreeSim)
library(picante)
library(reshape)
library(diversitree)
library(e1071)
library(phytools)
library(apTreeshape)
library(plyr)

rm(list=ls())

### Function range data between zero and one
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

### PSV function by Brunno Oliveira
# Modified from psv function in {picante}
psv.b<-function (samp, tree) {
  Cmatrix <- vcv.phylo(tree, corr = F)
  SR <- rowSums(data.frame(samp))
  nlocations <- 1
  nspecies <- length(tree$tip.label)
  index <- seq(1:nrow(Cmatrix))
  n <- length(index)
  C <- Cmatrix[index, index]
  PSV <- (n * sum(diag(as.matrix(C))) - sum(C))/(n * (n - 1))
  PSV
}

### MRD function by Brunno Oliveira
# Modified from ELIOT MILLER (in http://www.umsl.edu/~emmq7/Menu/Rphylo/MRD.R)
MRD.b <- function(phy) 
{
  phylo.bl1 <- compute.brlen(phy, 1)
  all.dist <- dist.nodes(phylo.bl1)
  root.dist <- all.dist[length(phy$tip.label)+1, 1:length(phy$tip.label)]
  tips.to.root <- data.frame(tipnames=phy$tip.label,root.dist)
  mrd <- mean(tips.to.root$root.dist)
  return(mrd)
}

# Phyvars
phyvars <- c("SPD","PD","age","MRD","MRD.time","ED","ED2","RBL",
       "DR","MPD","MNTD","PSV","GAM","IMY","IMP","DivRates","DivRates2")

#### Set WD
setwd('~/Dropbox/Doutorado Brunno/Manuscritos/Chap1 Age and FD/TreeSim')

#### load image
load("simTree.RData")

################################################################

#### SIMULATE PHYLOGENETIC TREES

#1) EFFECT OF TIME

# set parameters
ages <- rnorm(1000, 100, 20)
N <- 5000 # number of simulation
n <- 100 # number of taxa (tips)

# sample ages
tsamp<-sample(ages,N,replace = T)

# Create table for storing results
model1<-data.frame(matrix(data=NA,nrow=N,ncol=length(phyvars)+2))
colnames(model1)<-c("age","run",phyvars)


# Simulate N trees under a uniform birth-death process
for(i in 1:N){
  cat("\r",i,"of", N)
  time<-tsamp[i]
  trx <- sim.bd.taxa.age(n=100, numbsim=1, lambda=log(n/2)/time, mu=0, age=time, mrca=TRUE)
  trx <- trx[[1]]
  
  # get age for each species
  spp.ages<-data.frame(matrix(data=NA,nrow=length(trx$tip.label),ncol=2))
  colnames(spp.ages)<-c("Species","age")
  max<-max(branching.times(trx))
  for (j in 1:length(trx$tip.label)){
    Spp<-trx$tip.label[j]
    spp.ages$age[j] <- trx$edge.length[which.edge(trx,Spp)]
    spp.ages$Species[j]<-Spp
  }
  
  # get ED for each species
  spp.ED <- evol.distinct(trx, type = c("fair.proportion"), ### Species evolutionary distinctiveness (ED)
              scale = TRUE, use.branch.lengths = TRUE)
  
  # needed for MPD, MNTD and PSV  
  h<-data.frame(rep(1,length(trx$tip.label)))
  rownames(h)<-trx$tip.label
  h<-t(h)
  
  model1$run[i] <- i
  model1$SPD[i] <- length(trx$tip.label) ### Species richness
  model1$PD[i] <- sum(trx$edge.length) ### Phylogenetic diversity (PD)
  model1$MPD[i] <- mpd(h,cophenetic(trx)) ### Mean phylogenetic distance (MPD)
  model1$MNTD[i] <- mntd(h,cophenetic(trx)) ### Mean nearest taxon distance (MNTD)
  model1$PSV[i] <- psv.b(h,trx) ### Phylogenetic species variability (PSV)

  model1$MRD[i] <- MRD.b(trx) ### Mean root distance (MRD)  The mean number of nodes that separating species from the root of their tree. 
  model1$MRD.time[i] <- mean(max(branching.times(trx))-spp.ages$age) ### Mean root distance (MRD.time)  How far from the base of the tree species arise.
  
  model1$ED[i] <- mean(spp.ED[,2]) ### Species evolutionary distinctiveness (ED) 
  model1$ED2[i] <- mean(spp.ED[,2]) ### Needed to repeat - matching tables of results
  model1$DivRate[i] <- mean(1/spp.ED[,2]) ### Species-level diversification rate (DivRate)
  model1$DivRate2[i] <- mean(1/spp.ED[,2]) ### Needed to repeat - matching tables of results
  
  model1$RBL[i] <- mean(trx$edge.length/(as.numeric(branching.times(trx)[1])))  ### Relative branch length (RBL)
  
  model1$DR[i] <- log(model1$SPD[i])/time ### Diversification rates (DR)
  model1$GAM[i] <- gamStat(branching.times(trx),return.list=FALSE) ### Gamma statistics (GAM)  
  model1$IMY[i] <- colless(as.treeshape(trx),norm="yule") ### Colless Yule model (IMY) 
  model1$IMP[i] <- colless(as.treeshape(trx),norm="pda") ### Colless PDA model (IMP) 
  model1$age[i] <- time
  }


################################################################
#2) EFFECT OF DIVERSIFICATION RATES

# set parameters
lamb = rexp(10000,rate=2)
lamb = range01(lamb)
lamb = lamb+0.0001
N <- 5000 # number of simulation
n <- 100 # number of taxa (tips)

#simulate div values
tsamp=sample(lamb,N,replace = T)

# Create table for storing results
model2<-data.frame(matrix(data=NA,nrow=N,ncol=length(phyvars)+2))
colnames(model2)<-c("age","run",phyvars)

# Simulate N trees under a uniform birth-death process

rate<-rep(NA,N)
for(i in 1:N){
  
  cat("\r",i,"of", N)
  ts<-tsamp[i]
  trx <- sim.bd.taxa.age(n=100, numbsim=1, lambda=ts, mu=0, age=40, mrca=TRUE)
  trx <- trx[[1]]
  
  # get age for species
  spp.ages<-data.frame(matrix(data=NA,nrow=length(trx$tip.label),ncol=2))
  colnames(spp.ages)<-c("Species","age")
  max<-max(branching.times(trx))
  for (j in 1:length(trx$tip.label)){
  Spp<-trx$tip.label[j]
  spp.ages$age[j] <- trx$edge.length[which.edge(trx,Spp)]
  spp.ages$Species[j]<-Spp
  }
  
  # get ED for each species
  spp.ED <- evol.distinct(trx, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
              scale = TRUE, use.branch.lengths = TRUE)
  
    
  h<-data.frame(rep(1,length(trx$tip.label)))
  rownames(h)<-trx$tip.label
  h<-t(h)
  
  rate[i] <- ts
  model2$run[i] <- i
  model2$SPD[i] <- length(trx$tip.label) ### Species richness
  model2$PD[i] <- sum(trx$edge.length) ### Phylogenetic diversity (PD)
  model2$MPD[i] <- mpd(h,cophenetic(trx)) ### Mean phylogenetic distance (MPD)
  model2$MNTD[i] <- mntd(h,cophenetic(trx))  ### Mean nearest taxon distance (MNTD)
  model2$PSV[i] <- psv.b(h,trx) ### Phylogenetic species variability (PSV)
  
  model2$MRD[i] <- MRD.b(trx) ### Mean root distance (MRD)  The mean number of nodes that separating species from the root of their tree. 
  model2$MRD.time[i] <- mean(max(branching.times(trx))-spp.ages$age) ### Mean root distance (MRD.time)  How far from the base of the tree species arise.
  
  model2$ED[i] <- mean(spp.ED[,2]) ### Species evolutionary distinctiveness (ED) 
  model2$ED2[i] <- mean(spp.ED[,2]) ### Needed to repeat - matching tables of results
  model2$DivRate[i] <- mean(1/spp.ED[,2]) ### Species-level diversification rate (DivRate)
  model2$DivRate2[i] <- mean(1/spp.ED[,2]) ### Needed to repeat - matching tables of results
  
  model2$RBL[i] <- mean(trx$edge.length/(as.numeric(branching.times(trx)[1]))) ### Relative branch length (RBL)
  
  model2$DR[i] <- log(model2$SPD[i])/max(branching.times(trx)) ### Diversification rates (DR)
  model2$GAM[i] <- gamStat(branching.times(trx),return.list=FALSE) ### Gamma statistics (GAM)  
  model2$IMY[i] <- colless(as.treeshape(trx),norm="yule") ### Colless Yule model (IMY) 
  model2$IMP[i] <- colless(as.treeshape(trx),norm="pda") ### Colless PDA model (IMP) 
  
  model2$age[i] <- as.numeric(branching.times(trx)[1]) ### max lineage age 
}


################################################################
#3) EFFECT OF SPECIES RICHNESS
# We supose that the greater amount of prunning in a phylogenetic tree greater the greater 
# probability of the resulting tree be composed by species with deep evolutionary relationships. 
# Conversely, lower levels of tree prunning is likely to comprise more complete clades'
# histories and concentrate nodes towards the tips (high div and young assemblages). 
# The objective here is to find a phylometric that is not affected by the level of 
# tree prunning (species richness).

# set parameters
N <- 5000 # number of simulations
n <- rnorm(10000,mean=100,sd=38) # number of taxa (tips) * Parameters taken from the observed distribution of species richness values considering 1x1 degree resolution worldwide gridded data for mammals.
n = n+(-1*min(n))
n = n +10 # We consider that it is better to have at least 10 species to calculate phymetrics.

# Simulate a tree in parameters to the actual mammalian tree (Hedges et al., 2015)
trxs <- sim.bd.taxa.age(n=5000, numbsim=1, lambda=0.2, mu=0.14, age=180, mrca=T)
trxs <- trxs[[1]]

# get ED for species
spp.ED <- evol.distinct(trxs, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
            scale = TRUE, use.branch.lengths = TRUE)

# get age for each species (complete tree)
spp.ages.c<-data.frame(matrix(data=NA,nrow=length(trxs$tip.label),ncol=2))
colnames(spp.ages.c)<-c("Species","age")
max<-max(branching.times(trxs))
for (j in 1:length(trxs$tip.label)){
  Spp<-trxs$tip.label[j]
  spp.ages.c$age[j] <- trxs$edge.length[which.edge(trxs,Spp)]
  spp.ages.c$Species[j]<-Spp
}

# Create table for storing results
model3<-data.frame(matrix(data=NA,nrow=N,ncol=length(phyvars)+2))
colnames(model3)<-c("age","run",phyvars)

# Sample the tree to obtain assemblages with different richness
for(i in 1:N){
  
  cat("\r",i,"of", N)
  
  k<- as.integer(sample(n,1,replace = T)) #n species
  
  trx<-drop.tip(trxs,sample(1:5000,5000-k)) # sample tips to obtain assemblages of 100 tips
  
  # get age for each species (pruned tree)
  spp.ages.p<-data.frame(matrix(data=NA,nrow=length(trx$tip.label),ncol=2))
  colnames(spp.ages.p)<-c("Species","age")
  max<-max(branching.times(trx))
  for (j in 1:length(trx$tip.label)){
    Spp<-trx$tip.label[j]
    spp.ages.p$age[j] <- trx$edge.length[which.edge(trx,Spp)]
    spp.ages.p$Species[j]<-Spp
  }
  
  subED<-subset(spp.ED,Species%in%trx$tip.label)
  subAge<-subset(spp.ages.c,Species%in%trx$tip.label)
  
  h<-data.frame(rep(1,length(trx$tip.label)))
  rownames(h)<-trx$tip.label
  h<-t(h)
  
  model3$run[i] <- i
  model3$SPD[i] <- length(trx$tip.label) ### Species richness
  model3$PD[i] <- sum(trx$edge.length) ### Phylogenetic diversity (PD)
  model3$MPD[i] <- mpd(h,cophenetic(trx)) ### Mean phylogenetic distance (MPD)
  model3$MNTD[i] <- mntd(h,cophenetic(trx)) ### Mean nearest taxon distance (MNTD)
  model3$PSV[i] <- psv.b(h,trx) ### Phylogenetic species variability (PSV)
  
  model3$MRD[i] <- MRD.b(trx) ### Mean root distance (MRD)  The mean number of nodes that separating species from the root of their tree. 
  model3$MRD.time[i] <- mean(max(branching.times(trx))-spp.ages.p$age) ### Mean root distance (MRD.time)  How far from the base of the tree species arise.
  
  ### Species evolutionary distinctiveness (ED) 
  #1) Prune the tree first and then pass the tree in
  model3$ED[i] <- mean(evol.distinct(trx, type = c("fair.proportion"), scale = TRUE, use.branch.lengths = TRUE)[,2]) ## Pruned tree
  #2) Subset the resulting vector
  model3$ED2[i] <- mean(subED[,2])
  
  ### Species-level diversification rate (DivRate)
  #1) Prune the tree first and then pass the tree in
  model3$DivRate[i] <- 1/model3$ED[i] 
  #2) Subset the resulting vector
  model3$DivRate2[i] <- mean(1/subED[,2])
  
  model3$RBL[i] <- mean(trx$edge.length/(as.numeric(branching.times(trx)[1]))) ### Relative branch length (RBL)
  
  model3$DR[i] <- log(model3$SPD[i])/max(branching.times(trx)) ### Diversification rates (DR)
  model3$GAM[i] <- gamStat(branching.times(trx),return.list=FALSE) ### Gamma statistics (GAM)  
  model3$IMY[i] <- colless(as.treeshape(trx),norm="yule") ### Colless Yule model (IMY) 
  model3$IMP[i] <- colless(as.treeshape(trx),norm="pda") ### Colless PDA model (IMP) 
  
  model3$age[i] <- max(branching.times(trx)) ### max lineage age 
}


save.image("simTree.RData")


#### Create table of results

#p-values
results1<-data.frame(matrix(data=NA,nrow=length(phyvars),ncol=3))
colnames(results1)<-c("Richness","Tree Depth","Div. Rates")
rownames(results1)<-phyvars

model2<-cbind(model2,lamb=rate)

models<-c("SPD","age","lamb")
models_files<-list(model3,model1,model2)

for(j in 1:length(models)){
  for(i in 1:length(phyvars)){
  x<-as.numeric(models_files[[j]][paste(models[j])][[1]])
  y<-as.numeric(models_files[[j]][paste(phyvars[i])][[1]])
  test <- cor.test(x,y,method="spearman")
  results1[i,j]<-ifelse(round(test$p.value,3)<0.001,"<0.001",paste(round(test$p.value,3)))
  }
}

# Pearson's R values
results2<-data.frame(matrix(data=NA,nrow=length(phyvars),ncol=3))
colnames(results2)<-c("Richness","Tree Depth","Div. Rates")
rownames(results2)<-phyvars

models<-c("SPD","age","lamb")
models_files<-list(model3,model1,model2)

for(j in 1:length(models)){
  for(i in 1:length(phyvars)){
  x<-as.numeric(models_files[[j]][paste(models[j])][[1]])
  y<-as.numeric(models_files[[j]][paste(phyvars[i])][[1]])
  test <- cor.test(x,y,method='spearman')
  results2[i,j]<-round(test$estimate,2)
  }
}

results<-data.frame(results2[,1],results1[,1],results2[,2],results1[,2],results2[,3],results1[,3])
rownames(results)<-rownames(results1)
colnames(results)<-c("Richness Pearson's R","Richness P","Time Pearson's R","Time P","Div. Rates Pearson's R","Div. Rates P")


write.table(results,"results_spearman.csv",sep=",")

save.image("simTree.RData")

