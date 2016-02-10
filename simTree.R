
################################################################
#### Simulate phylogenetic trees and calculate phylometrics
#### Brunno Oliveira, 2016
#### Universidade Federal do Rio Grande do Norte - Brasil
#### Stony Brook University - USA
#### Contact brunno.oliveira@me.com for any information.
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
library(TreeSim)

setwd("/home/brunno/Dropbox/Doutorado Brunno/Manuscritos/Chap1 Age and FD/TreeSim/GitHub_simulation_phylometrics")

#### Load functions:

### Function range data between zero and one
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

### phylometrics function
source("phylometrics_function.R")

#### Simulate phylogenetic trees:

# 1) Manipulate richness:

# Phyvars
phyvars <- c("SR","PD","MRD","MRD.time","ED","RBL",
             "DR","MPD","MNTD","PSV","GAM","IMY","IMP","DivRate",
             "maxage", "ED2", "DivRate2", "Sp.Ages2")


# We supose that the greater amount of prunning in a phylogenetic tree greater the greater 
# probability of the resulting tree be composed by species with deep evolutionary relationships. 
# Conversely, lower levels of tree prunning is likely to comprise more complete clades'
# histories and concentrate nodes towards the tips (high div and young assemblages). 
# The objective here is to find a phylometric that is not affected by the level of 
# tree prunning (species richness).

# set parameters
N <- 50 # number of simulations
n <- rnorm(10000,mean=100,sd=38) # number of taxa (tips) * Parameters taken from the observed distribution of species richness values considering 1x1 degree resolution worldwide gridded data for mammals.
n = n+(-1*min(n))
n = n +10 # We consider that it is better to have at least 10 species to calculate phymetrics.

# Simulate a tree in parameters to the actual mammalian tree (Hedges et al., 2015)
trxs <- sim.bd.taxa.age(n=500, numbsim=1, lambda=0.2, mu=0.14, age=180, mrca=T)
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


# Sample the tree to obtain assemblages with different richness
for(i in 1:N) {  cat("\r",i,"of", N)
  
  k<- as.integer(sample(n,1,replace = T)) #n species
  
  trx<-drop.tip(trxs,sample(1:5000,5000-k)) # sample tips to obtain assemblages of n tips
  
  subED<-subset(spp.ED,Species%in%trx$tip.label)
  subAge<-subset(spp.ages.c,Species%in%trx$tip.label)
  
  run <- i
  
  ### Manipulated parameter
  par <- length(trx$tip.label)
  
  ### Species ages (Sp.Ages2) Subset the resulting vector
  Sp.Ages2 <- mean(subAge[,2])
  
  ### Species evolutionary distinctiveness (ED2) Subset the resulting vector
  ED2 <- mean(subED[,2])
  
  ### Species-level diversification rate (DivRate) Subset the resulting vector
  DivRate2 <- mean(1/subED[,2])
  
  ### calculate phylometrics
  phylo.res <- phylometrics(trx)
  
  ### save
  if(i==1){
    model1 <- c(run=run, par=par, phylo.res, ED2=ED2, DivRate2=DivRate2, Sp.Ages2=Sp.Ages2)
  }
  else{
    model1 <- rbind(model1, 
                    c(run=run, par=par, phylo.res, ED2=ED2, DivRate2=DivRate2, Sp.Ages2=Sp.Ages2))
  }
}

model1 <- data.frame(model1)

# 2) Manipulate time (age):

# set parameters
ages <- rnorm(1000, 100, 20)
N <- 50 # number of simulation
n <- 100 # number of taxa (tips)

# sample ages
tsamp<-sample(ages,N,replace = T)

# Simulate N trees under a uniform birth-death process
for(i in 1:N){ cat("\r",i,"of", N)
  time<-tsamp[i]
  trx <- sim.bd.taxa.age(n=100, numbsim=1, lambda=log(n/2)/time, mu=0, age=time, mrca=TRUE)
  trx <- trx[[1]]
  
  subED<-subset(spp.ED,Species%in%trx$tip.label)
  subAge<-subset(spp.ages.c,Species%in%trx$tip.label)
  
  run <- i
  
  ### Manipulated parameter
  par <- time
  
  ### Species ages (Sp.Ages2) Subset the resulting vector
  Sp.Ages2 <- NA
  
  ### Species evolutionary distinctiveness (ED2) Subset the resulting vector
  ED2 <- NA
  
  ### Species-level diversification rate (DivRate) Subset the resulting vector
  DivRate2 <- NA
  
  ### calculate phylometrics
  phylo.res <- phylometrics(trx)
  
  ### save
  if(i==1){
    model2 <- c(run=run, par=par, phylo.res, ED2=ED2, DivRate2=DivRate2, Sp.Ages2=Sp.Ages2)
  }
  else{
    model2 <- rbind(model2, 
                    c(run=run, par=par, phylo.res, ED2=ED2, DivRate2=DivRate2, Sp.Ages2=Sp.Ages2))
  }
}

model2 <- data.frame(model2)

# 3) Manipulate diversification rate:

# set parameters
lamb = rexp(10000, rate=2)
lamb = range01(lamb)
lamb = lamb+0.0001
N <- 50 # number of simulation
n <- 100 # number of taxa (tips)

#simulate div values
tsamp=sample(lamb,N,replace = T)

# Create table for storing results
model3<-data.frame(matrix(data=NA,nrow=N,ncol=length(phyvars)+2))
colnames(model3)<-c("age","run",phyvars)

# Simulate N trees under a uniform birth-death process

for(i in 1:N){ cat("\r",i,"of", N)
  
  ts<-tsamp[i]
  trx <- sim.bd.taxa.age(n=100, numbsim=1, lambda=ts, mu=0.14, age=40, mrca=TRUE) # change mu to 0.14
  trx <- trx[[1]]
  
  subED<-subset(spp.ED,Species%in%trx$tip.label)
  subAge<-subset(spp.ages.c,Species%in%trx$tip.label)
  
  run <- i
  
  ### Manipulated parameter
  par <- ts
  
  ### Species ages (Sp.Ages2) Subset the resulting vector
  Sp.Ages2 <- NA
  
  ### Species evolutionary distinctiveness (ED2) Subset the resulting vector
  ED2 <- NA
  
  ### Species-level diversification rate (DivRate) Subset the resulting vector
  DivRate2 <- NA
  
  ### calculate phylometrics
  phylo.res <- phylometrics(trx)
  
  ### save
  if(i==1){
    model3 <- c(run=run, par=par, phylo.res, ED2=ED2, DivRate2=DivRate2, Sp.Ages2=Sp.Ages2)
  }
  else{
    model3 <- rbind(model3, 
                    c(run=run, par=par, phylo.res, ED2=ED2, DivRate2=DivRate2, Sp.Ages2=Sp.Ages2))
  }
}

model3 <- data.frame(model3)



#### Table of results:

#p-values
results1<-data.frame(matrix(data=NA,nrow=length(phyvars),ncol=3))
colnames(results1)<-c("Richness","Tree Depth","Div. Rates")
rownames(results1)<-phyvars

models_files<-list(model1,model2,model3)

for(j in 1:length(models_files)){
  for(i in 1:length(phyvars)){
    x<-as.numeric(models_files[[j]]$par)
    y<-as.numeric(models_files[[j]][paste(phyvars[i])][[1]])
    test <- cor.test(x,y,method="spearman", exact=FALSE)
    results1[i,j]<-ifelse(round(test$p.value,3)<0.001,"<0.001",paste(round(test$p.value,3)))
  }
}

# Pearson's R values
results2<-data.frame(matrix(data=NA,nrow=length(phyvars),ncol=3))
colnames(results2)<-c("Richness","Tree Depth","Div. Rates")
rownames(results2)<-phyvars

for(j in 1:length(models_files)){
  for(i in 1:length(phyvars)){
    x<-as.numeric(models_files[[j]]$par)
    y<-as.numeric(models_files[[j]][paste(phyvars[i])][[1]])
    test <- cor.test(x,y,method='spearman')
    results2[i,j]<-round(test$estimate,2)
  }
}

results<-data.frame(results2[,1],results1[,1],results2[,2],results1[,2],results2[,3],results1[,3])
rownames(results)<-rownames(results1)
colnames(results)<-c("Richness - rho","Richness - P-value","Time - rho","Time - P-value","Div. Rate - rho","Div. Rate - P-value")
results <- results[-1,]

write.table(results,"results_spearman.csv",sep=",")

save.image("simTree.RData")

