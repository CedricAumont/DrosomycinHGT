################################################################################
#
#       Title:    Plotting Drs against Blattodea/Mantodea phylogeny (Fig2)
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

setwd("~/Projects/Drosomycin/Figure2")

#----------- Tree importation --------------------------------------------------


tree = read.tree("Figure2-Blattodea_Mantodea_tree.tre")


#----------- Data set importation ----------------------------------------------

data = read.csv("data-file-Figure2.csv",
                header = T, 
                sep = ",",dec = ".",
                stringsAsFactors = T) #turn strings into factors

str(data)
summary(data)


#----------- Check concordance tree-data ---------------------------------------

check <- name.check(phy = tree, data = data, 
                    data.names = data$Species_name)

check

#----------- Figure ggtree -----------------------------------------------------


p <- ggtree(tree) %<+% data + xlim(-.10, 21)+
  geom_tippoint(aes(colour=Drs_occurence, size = BUSCO..C.,shape = Data_type))+
  geom_tiplab(size = 2.4,offset = 0.3, fontface = 3)+
  guides(colour=guide_legend(title = "Drosomycin"))+
  scale_color_manual(values = c("red2","cyan2","grey60"),
                     labels=c( "Absent","Present", "Uncertain"))+
  scale_size_continuous(range = c(.5,5), name = "BUSCO score",
                        guide = guide_legend(reverse = T))+ 
  geom_strip('Nasutitermes_lujae', 'Cryptocercus_punctulatus', barsize=1, color='turquoise2', 
             label = "Isoptera\n+ Cryptocercidae", offset=4.2,offset.text=.4)+ 
  geom_strip('Blatta_orientalis', 'Hebardina_concinna', barsize=1, color='turquoise2', 
             label = "Blattidae", offset=4.2,offset.text = .4, extend = .5)+ 
  geom_strip('Nauphoeta_cinerea', 'Rhabdoblatta_sp', barsize=1, color='deepskyblue3', 
             label = "Blaberidae", offset=4.2,offset.text = .4, extend = .2)+ 
  geom_strip('Hemithyrsocera_sp', 'Margattea_sp', barsize=1, color='deepskyblue3', 
             label = "Ectobiidae", offset=4.2,offset.text = .4, extend = .1)+ 
  geom_strip('Eupolyphaga_sp', 'Eupolyphaga_sp', barsize=1,
             color='maroon2', 
             label = "Corydiidae", offset=4.2,offset.text = .4, extend = .25)+
  geom_strip('Tenodera_sinensis', 'Metallyticus_violacea', barsize=1,
           color='darkred', 
           label = "Mantodea", offset=4.2,offset.text = .4, extend = .0)

p

#END#
