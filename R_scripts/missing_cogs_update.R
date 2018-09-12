#these are the cogs that were described in the science 2010 paper... there not exactly cogs

library(ggplot2)
library(gplots)
library(reshape)
# source("https://bioconductor.org/biocLite.R")
# biocLite("ComplexHeatmap")
# biocLite("ggtree")
library(ComplexHeatmap)
library(tidyverse)
library(ggtree)
#devtools::install_github("guiastrennec/ggplus")
library(ggplus)

cogs <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/10_missing_cogs/combined_table_for_R.csv.2', sep = ' ', header = F)
colnames(cogs) <- c('species','yhom','count')
cogs.c <- cast(cogs, species ~ yhom)
cogs.c.m <- as.matrix(cogs.c, dimnames = list(cogs.c$species))
cogs.c.m <- cogs.c.m[-1,]
colnames(cogs.c.m) <- colnames(cogs.c)[-1]
cogs.c.m.t <- t(cogs.c.m)


# Functional grouping of the genes

annot <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/10_missing_cogs/gene_list_r.txt.atrrib.sorted.2.1',sep = '\t', header = F)

storage.mode(cogs.c.m.t) <- 'integer'

#remove bgh v3
cogs.c.m.t <- cogs.c.m.t[,-1]
#rename the columns to avoid more figure editing later
colnames(cogs.c.m.t) <- c('B. graminis fsp hordeii', "B. graminis fsp tritici","Botrytis cinerea","Colletotrichum graminicola","Erysiphe necator",
                          "Erysiphe pulchra","Fusarium oxysporum","Glarea lozoyensis","Golovinomyces cichoracearum", "Magnaporthe oryzae","Marssonina brunnea","Neurospora crassa",
                          "Oidiodendron maius","Oidium neolycopersici", "Parauncinula septata", "Penicillium digitatum", "Phialocephala scopiformis",
                          "Phialocephala subalpina", "Pseudogymnoascus destructans", "Pseudogymnoascus verrucosus",
                          "Rhynchosporium agropyri", "Rhynchosporium commune", "Rhynchosporium secalis","Sclerotinia borealis", "Sclerotinia sclerotiorum","Verticillium dahliae", "Zymoseptoria tritici")      

Heatmap(cogs.c.m.t,
        col = colorpanel(2, low = 'whitesmoke', high = '#b2df8a'),
        cluster_rows = F,
        #cluster_columns = F,
        split =  annot$V2,
        gap = unit(1, "mm"),
        row_names_gp = gpar(fontsize = 6),
        show_heatmap_legend = T,
        row_title_rot = 0,
        row_title_gp = gpar(fontsize = 9),
        column_names_gp = gpar(fontsize = 8),
        name = 'CAGs'
                )



# This ok, but I want to visuallize all the trees for these Parau homologs to make sure that this make sense and that they are not contaminations from 
# other genomes

tree1 <- read.newick('/Users/lfrantzeskakis/projects/03_jp_mildews/10_missing_cogs/yeast_parau_homologs/OG0000054_tree.txt')

tree2 <- read.newick('/Users/lfrantzeskakis/projects/03_jp_mildews/10_missing_cogs/yeast_parau_homologs/alltrees.txt', keep.multi = T)
tree2 <- read.newick('/Users/lfrantzeskakis/projects/03_jp_mildews/10_missing_cogs/yeast_parau_homologs/alltrees_clean_labels.txt', keep.multi = T)

treenames <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/10_missing_cogs/yeast_parau_homologs/tree_names_yeast.csv', sep = ' ', header = F)
class(tree2) <- "multiPhylo"
names(tree2) <- treenames$V3
#names(tree2) <- paste0(rep('tree',57),c(1:57))



ggtree(tree2, branch.length="none", layout="daylight")+
  geom_tippoint(aes(color=label))+
  geom_tiplab(size =1)+
  geom_tippoint(aes(color=grepl('botrytis|PARAU|sclerotinia_scl',label),shape=grepl('PARAU',label)),size=3)+
  #scale_color_brewer(type='div', palette=2)+
  #theme(legend.position="right")+
  facet_wrap(~.id, scale="free")
  

pdf("Rplot12.pdf")

ggtree(tree2, branch.length="none", layout="slanted")+
  geom_tippoint(aes(color=label))+
  geom_tiplab(size =2, align = T, linesize = 0.1, hjust = 0.5)+
  geom_tippoint(aes(color=grepl('botrytis|PARAU|sclerotinia_scl',label),shape=grepl('PARAU',label)),size=1)+
  #facet_wrap(~.id, scale="free")
  #scale_color_brewer(type='div', palette=2)+
  #theme(legend.position="right")+
  facet_wrap_paginate(~.id, ncol = 2, nrow =2, page = 1, scale="free")+
  #coord_flip()+
  theme_minimal()

dev.off()

pdf("Rplot12.pdf")

ggtree(tree2, branch.length="none", layout="slanted")+
  geom_tippoint(aes(color=grepl('botrytis|PARAU|sclerotinia_scl',label),shape=grepl('PARAU',label)),size=1)+
  scale_color_brewer(palette='Paired')+
  facet_wrap(~.id, scale="free")+
  coord_flip()+
  theme_minimal()+
  theme(legend.title=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank())
  

dev.off()


###


ggtree(tree2, branch.length="none", layout="slanted")+
  facet_wrap(~.id, scale="free")+
  geom_tippoint(aes(color=grepl('botrytis|PARAU|sclerotinia',label),shape=grepl('PARAU',label)),size=3)+
  theme_minimal(legend.position="right")


p<- ggtree(tree2, branch.length="none", layout="rectangular")+
  facet_wrap(~.id, scale="free")+
  geom_tippoint(aes(color=grepl('botrytis|PARAU|sclerotinia',label),shape=grepl('PARAU',label)),size=1.5)
  
