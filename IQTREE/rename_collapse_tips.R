#Andre E. Moncrieff
#4 Feb 2025
#Renaming tips and collapsing nodes

#load necessary libraries
install.packages("phylotools")
library(ape)
library(phylotools)



# clear R's brains
rm(list = ls())


setwd("/Users/andremoncrieff/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/IQTREE/")

ntree <- read.tree(file = "xipho.min4.phy.varsites.phy.treefile")


#tips <- tree$tip.label
#newtips <- read.table("newnames.txt", sep = "\t", header = TRUE)
#oldtips <- read.table("oldnames.txt", sep = "\t", header = TRUE)
#typeof(newtips)
#typeof(oldtips)
#dat <- data.frame(oldtips, newtips)
#ntree <- sub.taxa.label(tree, dat)
#write.tree(ntree, "newtree.tre")


#Following after code here: http://evoslav.blogspot.com/2015/01/how-to-collapse-unsupported-branches-of.html

Badnodes <- which(as.numeric(ntree$node.label) < 50) + length(ntree$tip.label)

Badnodes_indexes <- c()
for(node in Badnodes){
  Badnodes_indexes <- c(Badnodes_indexes, which(ntree$edge[,2] == node))
}

ntree$edge.length[Badnodes_indexes] <- 0

tree_multi <- di2multi(ntree)

write.tree(tree_multi, file = "xipho_tree_polytomies.tre")





