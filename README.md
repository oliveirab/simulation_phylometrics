Simulates phylogenetic trees with varying degrees of tree prunning (=species richness), ages and diversification rates. 

1. For the pdf version of Supporting Information click [here](https://github.com/oliveirab/simulation_phylometrics/blob/master/simTree.pdf)    
2. For access to the [source code](https://github.com/oliveirab/simulation_phylometrics/blob/master/simTree.R)  
3. For access the code to calculate phylometric click [here](https://github.com/oliveirab/simulation_phylometrics/blob/master/phylometrics_function.R)


Uses the simulated phylogenies to calculate the following phylometrics:

- Phylogenetic diversity (PD)  Sum of all branch lengths in a tree.	Faith (1992)
- Mean phylogenetic distance (MPD)	Mean of phylogenetic distances between all species pairs within a tree.	Webb (2000)
- Mean nearest taxon distance (MNTD)	Mean phylogenetic distance between each species and its nearest neighbor in a tree.	Webb (2000)
- Phylogenetic species variability (PSV)	Expected variance among species in a neutrally evolving trait.	Helmus et al. (2007) 
- Mean root distance (MRD)	The mean number of nodes that separating species from the root of their tree. Alternatively, when time trees are available MRD can be calculated as the mean distance between each species nodes and the root of its phylogenetic tree (MRD.time). 	Kerr and Currie (1999); Hawkins et al. (2012)
- Species evolutionary distinctiveness (ED)  After dividing each branch length by the number of species subtending that branch, the obtained values are summed across all branches from which a species descends.  	Redding and Mooers (2006); Isaac et al. (2007)
- Species-level diversification rate (DivRate)	Inverse of the evolutionary distinctness (ED). Species rapidly diversifying will have short edge lengths shared among many species.	Jetz et al. (2012), Kembel et al. (2010)
- Relative branch length (RBL)	Approximates time between speciation events across a tree by calculating the mean ((stem age – crown age) / crown age).	Marin and Blair, In prep.
- Diversification rates (DR)	Species accumulation through time (log(clade diversity)/clade age).	Magallon and Sanderson (2001); Nee (2001); Phillimore et al. (2006)
- Gamma statistics (GAM)	The distribution of branching events throughout the tree. This phylometric distinguished between trees with relatively long inter-nodal distances towards the tips (tippy trees) and trees with relatively longer inter-nodal distances towards the root of the phylogeny (stemmy trees). 	Pybus and Harvey (2000)
- Colless	Measures the degree of tree imbalance. Here, we calculated standardized values of Colless according to the Yule and PDA models.	Colless (1982)




***
Contact brunno.oliveira@me.com for any further information.  

* This repository follows the principles of reproducible research ([Peng, 2011](http://www.sciencemag.org/content/334/6060/1226)).
