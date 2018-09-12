#Libraries
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicFeatures")
library(GenomicFeatures) 
library('ggplot2')
options( scipen=999)

# Suppl. Fig 1 - QC of the assembly and filtering

## bash part:
#blasting:
#/home/lf216591/miniconda2/bin/tblastn -query /hpcwork/lf216591/15_parauncinula_assembly_cleanup/quercus_prots/GCF_002906115.1_CorkOak1.0_protein.faa -db /work/lf216591/07_pleo_para_annot/01_parauncinula/01_paraunc_annotation/05_selection_by_bgh_homology/scaffolds.draft2.500bp.nobact.fa -outfmt 6 -evalue 10e-5 -max_target_seqs 1 -num_threads 12 > quercus_hits.tsv
#cut -f2,3,11 scaffolds.draft2.500bp.nobact.fa.tblastnenec12.txt | sed 's/_/ /g' |  cut -f4,6,7,8 > size_coverage_evalue.enec.tsv
#cut -f2,3,11 scaffolds.draft2.500bp.nobact.fa.tblastnbgh12.txt | sed 's/_/  /g' |  cut -f4,6,7,8 > size_coverage_evalue.bgh.tsv


sceb <- read.csv('size_coverage_evalue.bgh.tsv', sep = '\t', header = F)
scee <- read.csv('size_coverage_evalue.enec.tsv', sep = '\t', header = F)
sceq <- read.csv('size_coverage_evalue.quercus.tsv', sep = '\t', header = F)
sceb <- cbind(sceb,'B. graminis f.sp. hordeii') # bgh 
scee <- cbind(scee,'E. necator') # necator
sceq <- cbind(sceq,'Q. suber') # plant hits

colnames(sceb) <- c('Size (bp)','Coverage','ID (%)','evalue','species')
colnames(scee) <- c('Size (bp)','Coverage','ID (%)','evalue','species')
colnames(sceq) <- c('Size (bp)','Coverage','ID (%)','evalue','species')
sce<- rbind(sceb,scee,sceq)

## Plotting coverage, size and evalue for the tblastn BGH hits to the assembly
ggplot(sce)+
  geom_point(aes(size,coverage, colour = id), alpha=0.75)+
  scale_colour_gradient2(midpoint = 50, low = 'tomato', high = 'black')+
  scale_x_log10(breaks = c(0,1000,10000,100000,1000000))+
  scale_y_log10(breaks = c(seq(0,100,10),1,5,200,500))+
  facet_wrap(~species)

ggplot(sce)+
  geom_point(aes(`Size (bp)`,`Coverage`, colour =`ID (%)`), alpha=0.5)+
  scale_colour_gradient2( low = 'blue',mid = 'white', high = 'red', midpoint = 50)+
  scale_x_log10(breaks = c(0,1000,10000,100000,1000000))+
  scale_y_log10(breaks = c(seq(0,100,10),1,5,200,500))+
  facet_wrap(~species)+
  theme_minimal()

## In way what you see is that the Quercus hits are in contings <5 depth, and in the >30 its low identity hits. While in Bgh/Enec is kind of the oposite.
## To me it is ok enough to say that I cut out whatever is <20 as plant related and keep >20 as fungal/powdery mildew. Together with BUSCO it should be convincing but we'll see.



