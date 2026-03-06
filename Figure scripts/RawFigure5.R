################################################################################
#
#       Title:    Plotting drosomycin against Insect phylogeny (Fig5)
#       Project:  Drosomycin
#       Year:     2025
#
################################################################################

#-------------------------------------------------------------------------------
#----------- Introduction ------------------------------------------------------
#-------------------------------------------------------------------------------

#----------- Library installation -----------------------------------------------


#install.packages(c('tidyverse','ape','ggimage','patchwork','phytools'))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#  BiocManager::install("ggtree")


#----------- Library importation -----------------------------------------------

library(tidyverse)
library(ape)
library(ggimage)
library(patchwork)
library(phytools)
library(ggtree)
library(geiger)
library(ggplot2)

#----------- Set work directory ------------------------------------------------

setwd("~/Projects/Drosomycin")

#----------- Tree importation --------------------------------------------------


tree = read.tree("Figure6/Tree-Fig6.tre")


#----------- Dataset importation -----------------------------------------------

data = read.csv("Figure6/Insect_dataset_fig6.csv",
                header = T, 
                sep = ",",dec = ".",
                stringsAsFactors = T) #turn str into factors

str(data)
summary(data)


#----------- Check concordance tree-data ---------------------------------------

check <- name.check(phy = tree, data = data, 
                    data.names = data$Clades)

check

#----------- Figure ggtree -----------------------------------------------------


p <- ggtree(tree) %<+% data + xlim(-.10, 21)+
  geom_tippoint(aes(colour=Drs))+
  geom_tiplab(size = 2.4,offset = 0.15, fontface = 3)+
  guides(colour=guide_legend(title = "Drosomycin"))+
  scale_color_manual(values = c("red2","cyan2"),
                     labels=c("Not found", "Present"))

p

#NB depending on the colour you wish, you have to choose a
#different scale_color (or scale_fill if filling is needed).
# you can find differente scales here
# https://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually

