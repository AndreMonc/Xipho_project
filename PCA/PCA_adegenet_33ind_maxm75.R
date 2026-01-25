#Andre E. Moncrieff

#Adegenet,PCA
# PCA for Xiphorhynchus dataset
#Beginning on 6 April 2025

#Super helpful info here: https://tomjenkins.netlify.app/tutorials/r-popgen-getting-started/

#load necessary libraries
library("ape")
library("pegas")
library("seqinr")
library("ggplot2")
library("adegenet")
library("vcfR")
library("stats")
library("ade4")
library("RColorBrewer")
library("dplyr")
#display.brewer.all(colorblindFriendly = TRUE)
#display.brewer.pal(n=12, name = 'Paired')
#brewer.pal(n=10, name = 'Paired')



# clear R's brains
rm(list = ls())

setwd("/Users/andremoncrieff/Dropbox/Work/Postdoc/Manuscript--Xiphorhynchus/Adegenet_PCA")


vcf <- read.vcfR("cluster_vcf_final.recode.vcf", verbose = FALSE )

my_genind <- vcfR2genind(vcf)

data(my_genind)

X <- tab(my_genind, freq = TRUE, NA.method = "mean")

# Perform PCA
pca1 <- dudi.pca(df = X, scale = FALSE, scannf = FALSE, nf = 3)
barplot(pca1$eig[1:50], main = "PCA eigenvalues", col = heat.colors(50))

# Analyse how much percent of genetic variance is explained by each axis
# Plot saved as pdf
percent = pca1$eig/sum(pca1$eig)*100
barplot(percent, ylab = "Genetic variance explained by eigenvectors (%)", ylim = c(0,20),
        names.arg = round(percent, 1))

# Create a data.frame containing individual coordinates
ind_coords = as.data.frame(pca1$li)

# Rename columns of dataframe
colnames(ind_coords) = c("Axis1","Axis2","Axis3")

# Add a column containing individuals
ind_coords$Ind = indNames(my_genind)

# Access file with pop assignments
Pop_file = read.delim(file = "pop_assignments.txt", header = TRUE, sep = "\t")

# Add a column with the cluster IDs
ind_coords$Cluster = Pop_file$Population_assignment

# Calculate centroid (average) position for each population
centroid = aggregate(cbind(Axis1, Axis2, Axis3) ~ Cluster, data = ind_coords, FUN = mean)

# Add centroid coordinates to ind_coords dataframe
ind_coords = left_join(ind_coords, centroid, by = "Cluster", suffix = c("",".cen"))

# Define colour palette
cols = brewer.pal(5, "Set1")

# Custom x and y labels
xlab = paste("Axis 1 (", format(round(percent[1], 1), nsmall=1)," %)", sep="")
ylab = paste("Axis 2 (", format(round(percent[2], 1), nsmall=1)," %)", sep="")

# Custom theme for ggplot2
ggtheme = theme(axis.text.y = element_text(colour="black", size=12),
                axis.text.x = element_text(colour="black", size=12),
                axis.title = element_text(colour="black", size=12),
                panel.border = element_rect(colour="black", fill=NA, size=1),
                panel.background = element_blank(),
                plot.title = element_text(hjust=0.5, size=15) 
)

# Scatter plot axis 1 vs. 2
ggplot(data = ind_coords, aes(x = Axis1, y = Axis2))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  
# points
geom_point(aes(fill = Cluster), shape = 21, size = 3, show.legend = FALSE)+

# centroids
geom_label(data = centroid, aes(label = Cluster, fill = Cluster), size = 4, show.legend = FALSE)+

# colouring
scale_fill_manual(values = cols)+
scale_colour_manual(values = cols)+
  
# custom labels
labs(x = xlab, y = ylab)+
ggtitle("PCA 33 ind, 75 maxm")+
# custom theme
ggtheme


