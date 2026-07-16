#Andre E. Moncrieff
#4 May 2026
#Renaming tips and collapsing nodes

#load necessary libraries
#install.packages("phylotools")
library(ape)
library(phylotools)


# clear R's brains
rm(list = ls())


setwd("/Users/andremoncrieff/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/Xipho_revision/IQTREE2/output")

ntree <- read.tree(file = "xipho.min4.phy.varsites.phy.treefile")


#Following after code here: http://evoslav.blogspot.com/2015/01/how-to-collapse-unsupported-branches-of.html

Badnodes <- which(as.numeric(ntree$node.label) < 50) + length(ntree$tip.label)

Badnodes_indexes <- c()
for(node in Badnodes){
  Badnodes_indexes <- c(Badnodes_indexes, which(ntree$edge[,2] == node))
}

ntree$edge.length[Badnodes_indexes] <- 0

tree_multi <- di2multi(ntree)

# Remove outgroup individual
tree_multi <- drop.tip(tree_multi, "XELEGANS_MPEG75162_Mus")

# 
write.tree(tree_multi, file = "xipho_tree_polytomies_noXELEGANS.tre")







