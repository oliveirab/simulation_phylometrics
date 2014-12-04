
################################################################
#### Simulate phylogenetic trees and calculate phylometrics
#### Brunno Oliveira, 2014
#### Universidade Federal do Rio Grande do Norte - Brasil
#### Stony Brook University - USA
################################################################

#### load packages
library(geiger)
library(ape)
library(TreeSim)
library(picante)
library(reshape)
library(diversitree)
library(e1071)
library(phytools)

rm(list=ls())

### Function range data between zero and one
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

### Function PSV
psv.b<-function (samp, tree, compute.var = TRUE) {
  Cmatrix <- vcv.phylo(tree, corr = TRUE)
  SR <- rowSums(data.frame(samp))
  nlocations <- 1
  nspecies <- length(tree$tip.label)
  index <- seq(1:nrow(Cmatrix))
  n <- length(index)
  C <- Cmatrix[index, index]
  PSV <- (n * sum(diag(as.matrix(C))) - sum(C))/(n * (n - 1))
  PSV
}

### Function to pairs correlation
# P-value and R coefficient
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y,method="spearman")
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y,method="spearman")
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt))
  text(.5, .75, Signif)
}
# Apply smoth regression line
panel.smooth<-function (x, y, col = "black", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}
# Add histogram to the diagonal
panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="white", ...)
}

### Tip lengths function
Tip.lengths<-function(phy){
  phy<-phy
  tip.length<-as.data.frame(phy$edge.length)
  tip.length<-cbind(rownames(tip.length),tip.length)
  names(tip.length)<-c('parental.node','length')
  rownames(tip.length)<-NULL
  
  ed<-as.data.frame(phy$edge) 
  names(ed)<-c("parental.node","daughter.node")
  
  #select only the nodes from tips
  ed<-subset(ed, daughter.node%in%1:length(phy$tip.label))
  #add age to parental nodes of tips
  tip.length<-merge(ed,tip.length,by.y="parental.node")
  tip.length<-tip.length$length
  
  tip.length
}

#### Set WD
setwd('C:/Users/Brunno/Dropbox/Doutorado Brunno/Manuscritos/Chap3 Age and FD/TreeSim')

#### save image
load("simTree.RData")

################################################################

#### RUN SIMULATION OF PHYLOGENETIC TREES

#### PHYLOMETRICS:

#1) EFFECT OF AGE

# set parameters
ages <- rnorm(1000, 100, 20)
N <- 5000 # number of simulation
n <- 100 # number of taxa (tips)

# Create table for storing results
phyvars <- c("SPD","PHD","AGX","clades","clade.age","Maxclade.age","AGM","MRD","ED","RBL",
             "DR","SAG","MaxSAG","SAC","MPD","MNTD","PSV")

model1<-data.frame(matrix(data=NA,nrow=N,ncol=length(phyvars)+2))
colnames(model1)<-c("age","run",phyvars)


# Simulate N trees under a uniform birth-death process
for(i in 1:N){
    
    cat("\r",i,"of", N)
    t=sample(ages,1)
    trx <- sim.bd.taxa.age(n=100, numbsim=1, lambda=log(n/2)/t, mu=0, age=t, mrca=TRUE)
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
    
    # get DR for species
    spp.ED <- evol.distinct(trx, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
                            scale = TRUE, use.branch.lengths = TRUE)
    
        
    sub<-subset(spp.ages,Species%in%trx$tip.label)
    subED<-subset(spp.ED,Species%in%trx$tip.label)
    
    model1$run[i] <- i
    model1$SPD[i] <- length(trx$tip.label)
    model1$AGX[i] <- max(branching.times(trx)) ### max lineage age 
    model1$clades[i] <- length(getCladesofSize(trx, clade.size=2)) ### clade rich
    model1$clade.age[i] <- mean(sapply(1:model1$clades[i],function(i) 
      max(branching.times(getCladesofSize(trx, clade.size=2)[[i]])))) ### clade age
    model1$Maxclade.age[i] <- max(sapply(1:model1$clades[i],function(i) 
      max(branching.times(getCladesofSize(trx, clade.size=2)[[i]])))) ### clade age
    model1$PHD[i] <- sum(trx$edge.length) ### phylogenetic diversity
    model1$AGM[i] <- median(branching.times(trx)) ### median lineage age
    model1$MRD[i] <- median(max(branching.times(trx))-branching.times(trx)) ### median root distance - mean distance from each node to the root of the trx
    model1$ED[i] <- median(evol.distinct(trx, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
                                         scale = FALSE, use.branch.lengths = TRUE)[,2])
    model1$DR[i] <- mean(1/subED[,2]) # Diversification rates in Jetz et al.2012
    model1$SAG[i] <- mean(sub$age)
    model1$MaxSAG[i] <- max(sub$age)
    
    model1$SAC[i] <- log(model1$SPD[i])/t ### Species accumulation - Diversification rates
    model1$RBL[i] <- median((max(branching.times(trx))-branching.times(trx))/max(branching.times(trx))) # relative tip length - small relative tips represent high DIV rates
    
    h<-data.frame(rep(1,length(trx$tip.label)))
    rownames(h)<-trx$tip.label
    h<-t(h)
    
    model1$MPD[i] <- mpd(h,cophenetic(trx))
    model1$MNTD[i] <- mntd(h,cophenetic(trx))
    model1$PSV[i] <- psv.b(h,trx,compute.var=TRUE)
}


### Spearman correlations

pairs(data.frame(AGX=model1$AGX,clades=model1$clades,cladeages=model1$clade.age,
                 ED=model1$ED,PHD=model1$PHD,AGM=model1$AGM,MRD=model1$MRD,MaxSAG=model1$MaxSAG,
                 RBL=model1$RBL,SAC=model1$SAC,DR=model1$DR,SAG=model1$SAG,
                 MPD=model1$MPD,MNTD=model1$MNTD,PSV=model1$PSV),
      lower.panel=panel.cor, upper.panel=panel.smooth,diag.panel=panel.hist)

################################################################
#2) EFFECT OF DIVERSIFICATION RATES

# set parameters
lamb = rexp(10000,rate=2)
lamb = range01(lamb)
lamb = lamb+0.0001
N <- 5000 # number of simulation
n <- 100 # number of taxa (tips)

# Create table for storing results
phyvars <- c("SPD","PHD","AGX","clades","clade.age","Maxclade.age","AGM","MRD","ED","RBL",
             "DR","SAG","MaxSAG","SAC","MPD","MNTD","PSV")

model2<-data.frame(matrix(data=NA,nrow=N,ncol=length(phyvars)+2))
colnames(model2)<-c("age","run",phyvars)

# Simulate N trees under a uniform birth-death process


for(i in 1:N){
  
  cat("\r",i,"of", N)
  t=sample(lamb,1)
  trx <- sim.bd.taxa.age(n=100, numbsim=1, lambda=t, mu=0, age=40, mrca=TRUE)
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
  
  # get DR for species
  spp.ED <- evol.distinct(trx, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
                          scale = TRUE, use.branch.lengths = TRUE)
  
  
  sub<-subset(spp.ages,Species%in%trx$tip.label)
  subED<-subset(spp.ED,Species%in%trx$tip.label)
  
  model2$run[i] <- i
  model2$SPD[i] <- length(trx$tip.label)
  model2$lamb[i] <- t #lambda
  model2$clades[i] <- length(getCladesofSize(trx, clade.size=2)) ### clade rich
  model2$clade.age[i] <- mean(sapply(1:model2$clades[i],function(i) 
    max(branching.times(getCladesofSize(trx, clade.size=2)[[i]])))) ### clade age
  model2$Maxclade.age[i] <- max(sapply(1:model2$clades[i],function(i) 
    max(branching.times(getCladesofSize(trx, clade.size=2)[[i]])))) ### clade age
  model2$PHD[i] <- sum(trx$edge.length) ### phylogenetic diversity
  model2$AGM[i] <- median(branching.times(trx)) ### median lineage age
  model2$MRD[i] <- median(max(branching.times(trx))-branching.times(trx)) ### median root distance - mean distance from each node to the root of the trx
  model2$ED[i] <- median(evol.distinct(trx, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
                                       scale = FALSE, use.branch.lengths = TRUE)[,2])
  model2$DR[i] <- mean(1/subED[,2]) # Diversification rates in Jetz et al.2012
  model2$SAG[i] <- mean(sub$age)
  model2$MaxSAG[i] <- max(sub$age)
  
  model2$SAC[i] <- log(model2$SPD[i])/t ### Species accumulation - Diversification rates
  model2$RBL[i] <- median((max(branching.times(trx))-branching.times(trx))/max(branching.times(trx))) # relative tip length - small relative tips represent high DIV rates
  
  h<-data.frame(rep(1,length(trx$tip.label)))
  rownames(h)<-trx$tip.label
  h<-t(h)
  
  model2$MPD[i] <- mpd(h,cophenetic(trx))
  model2$MNTD[i] <- mntd(h,cophenetic(trx))
  model2$PSV[i] <- psv.b(h,trx,compute.var=TRUE)
}


### Spearman correlations

pairs(data.frame(lamb=model2$lamb,clades=model2$clades,cladeages=model2$clade.age,
                 ED=model2$ED,PHD=model2$PHD,AGM=model2$AGM,MRD=model2$MRD,MaxSAG=model2$MaxSAG,
                 RBL=model2$RBL,SAC=model2$SAC,DR=model2$DR,SAG=model2$SAG,
                 MPD=model2$MPD,MNTD=model2$MNTD,PSV=model2$PSV),
      lower.panel=panel.cor, upper.panel=panel.smooth,diag.panel=panel.hist)


################################################################
#3) EFFECT OF CLADE DIVERSITY
# Phylometrics could be affected by the tree richness given its affect on the straight 
# of tree prunning. The greater amount of prunning in a phylogenetic tree greater the
# greater probability of the resulting tree be composed by species with deep evolutionary
# relationships. Conversely, lower levels of tree prunning is likely to comprise more
# complete clades and concentrate nodes towards the tips (high div and young assemblages). 
# The objective here is to find a phylometric that is not affected by the level of 
# tree prunning (species richness).

# set parameters
N <- 5000 # number of simulation
n <- rnorm(10000,mean=100,sd=38) # number of taxa (tips)
n = n+(-1*min(n))
n = n +10 # It is better to have at least 10 species to run phymetrics.

# Simulate a tree in parameters to the actual mammal timetree
trxs <- sim.bd.taxa.age(n=5000, numbsim=1, lambda=0.2, mu=0.14, age=180, mrca=T)
trxs <- trxs[[1]]

# get age for species
spp.ages<-data.frame(matrix(data=NA,nrow=length(trxs$tip.label),ncol=2))
colnames(spp.ages)<-c("Species","age")
max<-max(branching.times(trxs))
for (i in 1:length(trxs$tip.label)){
  Spp<-trxs$tip.label[i]
  spp.ages$age[i] <- trxs$edge.length[which.edge(trxs,Spp)]
  spp.ages$Species[i]<-Spp
}

# get DR for species
spp.ED <- evol.distinct(trxs, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
                        scale = TRUE, use.branch.lengths = TRUE)

# Create table for storing results
phyvars <- c("SPD","PHD","AGX","clades","clade.age","Maxclade.age","AGM","MRD","ED","RBL",
             "DR","SAG","MaxSAG","SAC","MPD","MNTD","PSV")

model3<-data.frame(matrix(data=NA,nrow=N,ncol=length(phyvars)+2))
colnames(model3)<-c("age","run",phyvars)

# Sample the tree to obtain assemblages with different richness
for(i in 1:N){
  
  cat("\r",i,"of", N)
  
  k<- as.integer(sample(n,1)) #n species
  trx<-drop.tip(trxs,sample(1:5000,5000-k)) # sample tips to obtain assemblages of 100 tips
  
  sub<-subset(spp.ages,Species%in%trx$tip.label)
  subED<-subset(spp.ED,Species%in%trx$tip.label)
  
  model3$run[i] <- i
  model3$SPD[i] <- length(trx$tip.label)
  model3$AGX[i] <- max(branching.times(trx)) ### max lineage age 
  model3$clades[i] <- length(getCladesofSize(trx, clade.size=2)) ### clade rich
  model3$clade.age[i] <- mean(sapply(1:model3$clades[i],function(i) 
    max(branching.times(getCladesofSize(trx, clade.size=2)[[i]])))) ### clade age
  model3$Maxclade.age[i] <- max(sapply(1:model3$clades[i],function(i) 
    max(branching.times(getCladesofSize(trx, clade.size=2)[[i]])))) ### clade age
  model3$PHD[i] <- sum(trx$edge.length) ### phylogenetic diversity
  model3$AGM[i] <- median(branching.times(trx)) ### median lineage age
  model3$MRD[i] <- median(max(branching.times(trx))-branching.times(trx)) ### median root distance - mean distance from each node to the root of the trx
  model3$ED[i] <- median(evol.distinct(trx, type = c("fair.proportion"), ### Species' evolutionary distinctiveness
                                       scale = FALSE, use.branch.lengths = TRUE)[,2])
  model3$DR[i] <- mean(1/subED[,2]) # Diversification rates in Jetz et al.2012
  model3$SAG[i] <- mean(sub$age)
  model3$MaxSAG[i] <- max(sub$age)
  
  model3$SAC[i] <- log(model3$SPD[i])/t ### Species accumulation - Diversification rates
  model3$RBL[i] <- median((max(branching.times(trx))-branching.times(trx))/max(branching.times(trx))) # relative tip length - small relative tips represent high DIV rates
  
  h<-data.frame(rep(1,length(trx$tip.label)))
  rownames(h)<-trx$tip.label
  h<-t(h)
  
  model3$MPD[i] <- mpd(h,cophenetic(trx))
  model3$MNTD[i] <- mntd(h,cophenetic(trx))
  model3$PSV[i] <- psv.b(h,trx,compute.var=TRUE)
}


### Spearman correlations

pairs(data.frame(SPD=model3$SPD,clades=model3$clades,cladeages=model3$clade.age,
                 ED=model3$ED,PHD=model3$PHD,AGM=model3$AGM,MRD=model3$MRD,MaxSAG=model3$MaxSAG,
                 RBL=model3$RBL,SAC=model3$SAC,DR=model3$DR,SAG=model3$SAG,
                 MPD=model3$MPD,MNTD=model3$MNTD,PSV=model3$PSV),
      lower.panel=panel.cor, upper.panel=panel.smooth,diag.panel=panel.hist)

save.image("simTree.RData")


#### Create table of results
phyvars <- c("PHD","clades","clade.age","Maxclade.age","AGM","MRD","ED","RBL",
             "DR","SAG","MaxSAG","SAC","MPD","MNTD","PSV")

#p-values
results1<-data.frame(matrix(data=NA,nrow=length(phyvars),ncol=3))
colnames(results1)<-c("Tree Depth","Div. Rates","Richness")
rownames(results1)<-phyvars

models<-c("AGX","lamb","SPD")
models_files<-list(model1,model2,model3)

for(j in 1:length(models)){
  for(i in 1:length(phyvars)){
    x<-as.numeric(models_files[[j]][paste(models[j])][[1]])
    y<-as.numeric(models_files[[j]][paste(phyvars[i])][[1]])
    test <- cor.test(x,y,method="spearman")
    results1[i,j]<-ifelse(round(test$p.value,3)<0.001,"<0.001",paste(round(test$p.value,3)))
  }
}

# rho values
results2<-data.frame(matrix(data=NA,nrow=length(phyvars),ncol=3))
colnames(results2)<-c("Tree Depth","Div. Rates","Richness")
rownames(results2)<-phyvars

models<-c("AGX","lamb","SPD")
models_files<-list(model1,model2,model3)

for(j in 1:length(models)){
  for(i in 1:length(phyvars)){
    x<-as.numeric(models_files[[j]][paste(models[j])][[1]])
    y<-as.numeric(models_files[[j]][paste(phyvars[i])][[1]])
    test <- cor.test(x,y,method="spearman")
    results2[i,j]<-round(test$estimate,2)
  }
}

results<-data.frame(results2[,1],results1[,1],results2[,2],results1[,2],results2[,3],results1[,3])
rownames(results)<-rownames(results1)
colnames(results)<-c("Tree Depth Rho","Tree Depth P","Div. Rates Rho","Div. Rates P","Richness Rho","Richness P")


write.table(results,"results.csv",sep=",")

save.image("simTree.RData")

