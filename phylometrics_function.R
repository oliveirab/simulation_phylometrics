
################################################################
#### Simulate phylogenetic trees and calculate phylometrics
#### Brunno Oliveira, 2016
#### Universidade Federal do Rio Grande do Norte - Brasil
#### Stony Brook University - USA
#### Contact brunno.oliveira@me.com for any information.
################################################################

phylometrics <- function(phy){
  
  if(!is.element('picante', installed.packages()[,1])) { stop("Package not found. Install picante") }
  if(!is.element('apTreeshape', installed.packages()[,1])) { stop("Package not found. Install apTreeshape") }
  if(!is.element('laser', installed.packages()[,1])) { stop("Package not found. Install laser") }
  
  require("picante")
  require("apTreeshape")
  require("laser")
  
  ### PSV function by Brunno Oliveira .Modified from psv function in {picante}
  psv.b<-function (phy) {
    h<-data.frame(rep(1,length(phy$tip.label)))
    rownames(h)<-phy$tip.label
    h<-t(h)
    Cmatrix <- vcv.phylo(phy, corr = F)
    SR <- rowSums(data.frame(h))
    nlocations <- 1
    nspecies <- length(phy$tip.label)
    index <- seq(1:nrow(Cmatrix))
    n <- length(index)
    C <- Cmatrix[index, index]
    PSV <- (n * sum(diag(as.matrix(C))) - sum(C))/(n * (n - 1))
    return(PSV)
  }
  
  ### MRD function by Brunno Oliveira. Modified from ELIOT MILLER (downloaded from http://www.umsl.edu/~emmq7/Menu/Rphylo/MRD.R)
  MRD.b <- function(phy) {
    phylo.bl1 <- compute.brlen(phy, 1)
    all.dist <- dist.nodes(phylo.bl1)
    root.dist <- all.dist[length(phy$tip.label)+1, 1:length(phy$tip.label)]
    tips.to.root <- data.frame(tipnames=phy$tip.label,root.dist)
    mrd <- mean(tips.to.root$root.dist)
    return(mrd)
  }
  
  ### Species ages. The age at the split.
  Spp.Ages <- function(phy) {
    spp.ages <- data.frame(matrix(data=NA,nrow=length(phy$tip.label),ncol=2))
    colnames(spp.ages)<-c("Species","age")
    max<-max(branching.times(phy))
    for (j in 1:length(phy$tip.label)){
      Spp<-phy$tip.label[j]
      spp.ages$age[j] <- phy$edge.length[which.edge(phy,Spp)]
      spp.ages$Species[j]<-Spp
    }
    mean(spp.ages$age)
  }
  
  ### Evolutionary distinctiveness. Borrowed from picante
  ED <- function(phy){
    mean(evol.distinct(phy, type = c("fair.proportion"), scale = TRUE, use.branch.lengths = TRUE)[,2])
  }
  
  ### Species richness
  SR <- function(phy){
    length(phy$tip.label)
  }
  
  ### Phylogenetic diversity (Faith's PD)
  PD <- function(phy){
    sum(phy$edge.length)
  }
  
  ### Mean phylogenetic distance (MPD). Borrowed from picante
  MPD <- function(phy){
    h<-data.frame(rep(1,length(phy$tip.label)))
    rownames(h)<-phy$tip.label
    h<-t(h)
    mpd(h,cophenetic(phy))
  }
  
  ### Mean nearest taxon distance (MNTD)
  MNTD <- function(phy){
    h<-data.frame(rep(1,length(phy$tip.label)))
    rownames(h)<-phy$tip.label
    h<-t(h)
    mntd(h,cophenetic(phy)) 
  }
  
  ### Mean root distance (MRD.time)  How far from the base of the tree species arise.
  MRD.time <- function(phy){
    spp.ages <- data.frame(matrix(data=NA,nrow=length(phy$tip.label),ncol=2))
    colnames(spp.ages)<-c("Species","age")
    max<-max(branching.times(phy))
    for (j in 1:length(phy$tip.label)){
      Spp<-phy$tip.label[j]
      spp.ages$age[j] <- phy$edge.length[which.edge(phy,Spp)]
      spp.ages$Species[j]<-Spp
    }
    mean(max(branching.times(phy))-spp.ages$age) 
  }
  
  ### Species-level diversification rate (DivRate) (Jetz et al. 2012, Nature)
  DivRate <- function(phy){
    1/mean(evol.distinct(phy, type = c("fair.proportion"), scale = TRUE, use.branch.lengths = TRUE)[,2])
  }
  
  ### Relative branch length (RBL)
  RBL <- function(phy){
    mean(phy$edge.length/(as.numeric(branching.times(phy)[1])))
  }
  
  ### Diversification rates (DR)
  DR <- function(phy){
    SR <- length(phy$tip.label)
    t <- max(branching.times(phy))
    return(log(SR)/t)
  }
  
  ### Gamma statistics (GAM). Borrowed from laser package
  GAM <- function(phy){
    gamStat(branching.times(phy),return.list=FALSE)
  }
  
  ### Colless Yule model (IMY). Borrowed from apTreeshape package
  IMY <- function(phy){
    colless(as.treeshape(phy),norm="yule")
  }
  
  ### Colless PDA model (IMP). Borrowed from apTreeshape package
  IMP <- function(phy){
    colless(as.treeshape(phy),norm="pda")
  }
  
  ### Maximum lineage age
  maxage <- function(phy){
    max(branching.times(phy))
  }
  
  ### write results
  return(c(SR=SR(phy), PD=PD(phy), Sp.Ages=Spp.Ages(phy), MRD=MRD.b(phy), MRD.time=MRD.time(phy),
           ED=ED(phy), RBL=RBL(phy), DR=DR(phy), MPD=MPD(phy), MNTD=MNTD(phy), 
           PSV=psv.b(phy), GAM=GAM(phy), IMY=IMY(phy), IMP=IMP(phy),
           DivRate=DivRate(phy), maxage=maxage(phy)))
}
