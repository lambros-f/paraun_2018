# Libraries, working paths
setwd("~/Documents/GitHub/jap_pm_2018")
#devtools::install_github('cttobin/ggthemr')
options(scipen=999)
library(devtools)
library(ggplot2)
library(ggridges)
library(ggthemr)

#### Results Figures ####

#### Figure 1 - Gene density/Genome expansion of Bgh ####

#### Figure 1A - Gene density/Intergenic distance plots ####

ds <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/down.txt', sep = '\t', header = F)
us <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/up.txt', sep = '\t', header = F)
sps <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/m.res.per_contigs.high_confidence.gene.lst.fa.tmhmm.notmm', sep = '\t', header = F)

ds <- cbind(ds,'Downstream')
us <- cbind(us,'Upstream')
all <- merge(ds,us, by = 'V4')

colnames(all) <- c('gene1','scaffold1','start1','stop1','scaffold2','start2','stop2','gene2','distance_size_up','ds','scaffold1','start1','stop1','scaffold2','start2','stop2','gene2','distance_size_down','us')

all$distance_size_up <- abs(all$distance_size_up)
all$distance_size_down <- abs(all$distance_size_down)

all$group <- "Non SP"
all$group[all$gene1 %in% sps$V1] <- "SP"



ggplot(all, aes(x= distance_size_up,y=distance_size_down)) +
  #geom_hex()+
  stat_bin2d(binwidth=c(0.06, 0.06))+
  scale_fill_distiller(palette = "Spectral", name="Gene\ncount")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("5' prime intergenic length (bp)") +
  xlab("3' prime intergenic length (bp)") +
  geom_point(data=subset(all, group == "SP"), color = 'black', fill = 'white',shape = 21,alpha = 0.8, size = 2)+
  theme_minimal()

#### Figure 1B - TE distance from start/stop ####

#In bash, remove sub-categorization from TEs, for example make Tad1-12 to Tad1 in order to summarize better the results
#Also remove the entries that are at the end of scaffolds and denoted as '.'
#sed  -e 's/Gypsy4-.*       /Gypsy  /g' -e 's/Copia-.*        /Copia  /g' -e 's/Gypsy-.*      /Gypsy  /g' -e 's/HaTad1-.*     /HaTad1 /g' -e 's/Mariner-.*    /Mariner        /g' -e 's/Mariner5_AO.*        /Mariner       /g' -e 's/Tad1-.*       /Tad1   /g' down_te.txt | grep -v '       .       '  > down_te_fixed.txt
#sed  -e 's/Gypsy4-.*       /Gypsy  /g' -e 's/Copia-.*        /Copia  /g' -e 's/Gypsy-.*      /Gypsy  /g' -e 's/HaTad1-.*     /HaTad1 /g' -e 's/Mariner-.*    /Mariner        /g' -e 's/Mariner5_AO.*        /Mariner       /g' -e 's/Tad1-.*       /Tad1   /g' up_te.txt | | grep -v '       .       ' > up_te_fixed.txt


dste <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/down_te_fixed.txt', sep = '\t', header = F)
uste <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/up_te_fixed.txt', sep = '\t', header = F)
sps <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/m.res.per_contigs.high_confidence.gene.lst.fa.tmhmm.notmm', sep = '\t', header = F)

dste$group <- 'Downstream'
uste$group <- 'Upstream'
allte <- rbind(uste,dste)

colnames(allte) <- c('scaffold1','start1','stop1','gene','scaffold2','start2','stop2','te_type','distance','group')

allte$distance <- abs(allte$distance)

allte$secr <- "Non SP"
allte$secr[allte$gene %in% sps$V1] <- "SP"

# ggplot(allte) +
#   geom_freqpoly(aes(x=distance,colour=te_type),binwidth = 1000)+
#   facet_wrap(~group+secr, scales = 'free')+
#   ggtitle("Types repetive elements upstream and downstream of the genes" ) +
#   theme_minimal()  +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5)) +
#   labs(x= "Distance (bp)", y = "Frequency") +
#   coord_cartesian(xlim = c(-10000, 10000))+
#   scale_colour_brewer(palette = "Paired", name="Types")+
#   coord_cartesian(xlim = c(0, 50000))

ggthemr('dust')
ggplot(allte) +
  geom_point(aes(x=te_type,y=distance),alpha = 0.04, position = 'jitter')+
  geom_violin(aes(x=te_type,y=distance, fill=te_type))+
  facet_wrap(~group+secr)+
  scale_fill_brewer(palette = "Paired", name="Types")+
 # theme_minimal()  +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5)) +
  scale_y_log10()+
  ylab("Distance for start/stop codon (bp)") +
  xlab("Repetitive element type (bp)") 
  

#### Figure 1C - Intergenic distance in Bgh and Parauncinula ####

#make the table for Parauncinula
ds <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/down.txt', sep = '\t', header = F)
us <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/up.txt', sep = '\t', header = F)
sps <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/m.res.per_contigs.high_confidence.gene.lst.fa.tmhmm.notmm', sep = '\t', header = F)

ds$group <- 'Downstream'
us$group <- 'Upstream'
parau_interg <- rbind(us,ds)

colnames(parau_interg) <- c('scaffold1','start1','stop1','gene','scaffold2','start2','stop2','gene2','distance','group')

parau_interg$distance <- abs(parau_interg$distance)

parau_interg$secr <- "Non SP"
parau_interg$secr[parau_interg$gene %in% sps$V1] <- "SP"
parau_interg_final <- cbind.data.frame(parau_interg$gene,parau_interg$distance,parau_interg$group,parau_interg$secr)
colnames(parau_interg_final) <- c('gene','distance','group','secreted')
parau_interg_final$species <- 'Parauncinula septata'

#make the table for Bgh
bgh_interg <- read.csv('/Users/lfrantzeskakis/Documents/GitHub/blumeria_2017/useful_files/intergenic_spaces.bed2R', sep = '\t', header = F)
colnames(bgh_interg) <- c('Scaffold','Gene_start','Gene_end','Gene_name','Scaffold_inter','Intergen_start','Intergen_stop','UpDown')
table_up <- subset.data.frame(bgh_interg, UpDown > 0)
table_up$UpDown <- 'Upstream'
table_up$distance <- abs(table_up$Intergen_start - table_up$Intergen_stop)
table_up$secreted <- 'Non SP'

table_down <- subset.data.frame(bgh_interg, UpDown < 0)
table_down$UpDown <- 'Downstream'
table_down$distance <- abs(table_down$Intergen_start - table_down$Intergen_stop)
table_down$secreted <- 'Non SP'

bgh_interg.2 <- rbind(table_up, table_down)
bgh_interg_final <- cbind.data.frame(bgh_interg.2$Gene_name,bgh_interg.2$distance,bgh_interg.2$UpDown,bgh_interg.2$secreted)
colnames(bgh_interg_final) <- c('gene','distance','group','secreted')
bgh_interg_final$species <- 'Blumeria graminis fsp hordei'

bghsps <- read.csv('/Users/lfrantzeskakis/Documents/GitHub/blumeria_2017/useful_files/cseps_805.lst', sep = '\t', header = F)

bgh_interg_final$secreted <- as.character(bgh_interg_final$secreted)
bgh_interg_final$secreted[bgh_interg_final$gene %in% bghsps$V1] <- "SP"

data <- rbind(parau_interg_final,bgh_interg_final)

ggplot(data)+
  geom_violin(aes(y=distance,x=group))+
  facet_wrap(~species+secreted)+
  scale_y_log10()+
  theme_minimal()
  
ggthemr('dust')
ggplot(data)+
  geom_violin(aes(y=distance,x=group))+
  facet_wrap(~species+secreted, ncol = 4)+
  scale_y_log10()+
  ylab("Distance to the nearest gene (bp)") +
  xlab("") 
  
#### Figure 1B - CAZYmes ####

caz <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_cazy/all_results_for_R_summarized.csv', sep = ' ', header = F)
colnames(caz) <- c('CAZyme','Count','Species')


caz$Count <- as.numeric(caz$Count)
ggthemr('dust')
ggthemr('greyscale')
ggplot(caz, aes(caz$Species,caz$CAZyme))+
  geom_tile(aes(fill = caz$Count))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5))+
  xlab('')+
  ylab('CAZyme group')




#### Figure 2 - RIP ####


### Bimodal distribution of GC/AT

### Occultercut

#Get the data

ocutest.all <-''
for (filez in list.files('/Users/lfrantzeskakis/projects/03_jp_mildews/07_occultercut/', pattern = '*.compositionGC.txt')){
  filez <- paste('/Users/lfrantzeskakis/projects/03_jp_mildews/07_occultercut/',filez,sep = '')
  ocutest <- read.csv(filez, header = F, sep = ' ')
  ocutest$id <- filez
  ocutest.all <- rbind(ocutest.all,ocutest)
  
}

#Cleanup the table
ocutest.all$id <- sub('/Users/lfrantzeskakis/projects/03_jp_mildews/07_occultercut/','',ocutest.all$id, fixed = T)
ocutest.all$id <- sub('.compositionGC.txt','',ocutest.all$id, fixed = T)
ocutest.all$id <- sub('_GC.*','',ocutest.all$id )
ocutest.all$id <- sub('_',' ',ocutest.all$id )
ocutest.all$id <- sub('.fa','',ocutest.all$id )
ocutest.all <- ocutest.all[!apply(ocutest.all == "", 1, all),]
colnames(ocutest.all) <- c('GC','Portion','Species')
ocutest.all$GC <- as.numeric(ocutest.all$GC)
ocutest.all$Portion <- as.numeric(ocutest.all$Portion)

#plot
ggplot(ocutest.all)+
  geom_polygon(aes(GC,Portion), fill = 'tomato')+
  facet_wrap(~Species, scales = 'free_y')+
  theme_minimal()+
  ylab(label = 'Portion of the genome (%)')+
  xlab(label = 'GC content (%)')

#Make a simpler figure for the main figure

ocutest.all.s1 <- subset.data.frame(ocutest.all, Species == 'res.per contigs.high_confidence.lst')
ocutest.all.s2 <- subset.data.frame(ocutest.all, Species == 'bgh dh14_v4')
ocutest.all.s3 <- subset.data.frame(ocutest.all, Species == 'rhynchosporium commune')
ocutest.all.s4 <- subset.data.frame(ocutest.all, Species == 'marssonina brunnea')
ocutest.all.s <- rbind(ocutest.all.s1,ocutest.all.s2,ocutest.all.s3,ocutest.all.s4)
ocutest.all.s$Species <- gsub('bgh dh14_v4','Blumeria graminis f.sp. hordei',ocutest.all.s$Species )
ocutest.all.s$Species <- gsub('rhynchosporium commune','Rhynchosporium commune',ocutest.all.s$Species )
ocutest.all.s$Species <- gsub('marssonina brunnea','Marssonina brunnea',ocutest.all.s$Species )
ocutest.all.s$Species <- gsub('res.per contigs.high_confidence.lst','Parauncinula septata',ocutest.all.s$Species )


ggplot(ocutest.all.s)+
  geom_polygon(aes(GC,Portion), fill = 'tomato')+
  facet_wrap(~Species, scales = 'free_y')+
  theme_minimal()+
  ylab(label = 'Portion of the genome (%)')+
  xlab(label = 'GC content (%)')

### RIP Indexes and dinucleotide frequencies

### RIPCAL

#test with more indexes and dominance
#lets keep it simple and check only the tad1,gypsy,copia

ripidx <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/05_ripcal/01_leotio_dinucl_calc/05_with_dominance/summary.all.noavg.dom', header = T, sep = ';')
colnames(ripidx) <- c("species","family","TA/AT I1","(CpA+TpG)/(ApC+GpT) I2","(CpA+TpG)/TpA","(CpC + GpG)/(TpC + GpA)","CpG/(TpG + CpA)","(CpT + ApG)/(TpT + ApA)","counts_per_family")

ripidx$dom_ca_ta <- ripidx$`(CpA+TpG)/TpA` / (ripidx$`(CpC + GpG)/(TpC + GpA)` + ripidx$`CpG/(TpG + CpA)` + ripidx$`(CpT + ApG)/(TpT + ApA)`)
ripidx$dom_cc_tc <- ripidx$`(CpC + GpG)/(TpC + GpA)`  / ( ripidx$`(CpA+TpG)/TpA` + ripidx$`CpG/(TpG + CpA)` + ripidx$`(CpT + ApG)/(TpT + ApA)`)
ripidx$dom_cg_tg <- ripidx$`CpG/(TpG + CpA)`  / (ripidx$`(CpC + GpG)/(TpC + GpA)` + ripidx$`(CpA+TpG)/TpA` + ripidx$`(CpT + ApG)/(TpT + ApA)`)
ripidx$dom_ct_tt <-  ripidx$`(CpT + ApG)/(TpT + ApA)`  / (ripidx$`(CpC + GpG)/(TpC + GpA)` + ripidx$`CpG/(TpG + CpA)` + ripidx$`(CpA+TpG)/TpA`)


ripidx.m <- melt(ripidx)
ripidx.s <- cbind.data.frame(ripidx$species,ripidx$family,ripidx$`TA/AT I1`,ripidx$`(CpA+TpG)/(ApC+GpT) I2`,ripidx$dom_ca_ta,ripidx$dom_cc_tc,ripidx$dom_cg_tg,ripidx$dom_ct_tt)
colnames(ripidx.s) <- c('species','family','I1','I2','dom_ca','dom_cc','dom_cg','dom_ct')
#ripidx.s <- subset.data.frame(ripidx.s, species == c('neurospora_crassa.res.txt','botrytis_cinirea_b05','rhynchosporium_commune','blumeria_graminis_fsp_hordei','parauncinula_septata'))
ripidx.s1 <- subset.data.frame(ripidx.s, species == 'parauncinula_septata')
ripidx.s2 <- subset.data.frame(ripidx.s, species == 'rhynchosporium_commune')
ripidx.s3 <- subset.data.frame(ripidx.s, species == 'blumeria_graminis_fsp_hordei')
ripidx.s4 <- subset.data.frame(ripidx.s, species == 'marssonina_brunnea')
ripidx.s <- rbind(ripidx.s1,ripidx.s2,ripidx.s3,ripidx.s4)

#Index 1

ripidx.m <- melt(ripidx.s)
ripidx.controls <- subset(ripidx.m, family == 'control')
ripidx.m <- subset(ripidx.m,variable == 'I1')
ripidx.m <- subset(ripidx.m,family == c('Gypsy','Copia','Tad1'))
ripidx.controls <- subset(ripidx.controls, variable == 'I1')


ggplot(ripidx.m)+
  geom_boxplot(aes(x= family , y=value))+
  geom_point(aes(x= family , y=value), colour = 'grey', alpha = 0.4)+
  facet_wrap(~species , scales = 'fixed', ncol =6)+
  geom_hline(data = ripidx.controls, aes(yintercept = value) , linetype="dashed", colour = 'red')+
  geom_hline(yintercept = 0.89, linetype="dashed", colour = 'blue', alpha = 0.5)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5)) +
  theme_minimal()+
  ylab(label = 'TpA/ApT')+
  xlab(label = 'TE family')

#Index 2 

ripidx.m <- melt(ripidx.s)
ripidx.controls <- subset(ripidx.m, family == 'control')
ripidx.m <- subset(ripidx.m,variable == 'I2')
ripidx.m <- subset(ripidx.m,family == c('Gypsy','Copia','Tad1'))
ripidx.controls <- subset(ripidx.controls, variable == 'I2')


ggplot(ripidx.m)+
  geom_boxplot(aes(x= family , y=value))+
  facet_wrap(~species , scales = 'free', ncol =6)+
  geom_point(aes(x= family , y=value), colour = 'grey', alpha = 0.4)+
  geom_hline(data = ripidx.controls, aes(yintercept = value) , linetype="dashed", colour = 'red')+
  geom_hline(yintercept = 1.03, linetype="dashed", colour = 'blue')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5)) +
  theme_minimal()+
  ylab(label = '(CpA + TpG)/(ApC+ GpT)')+
  xlab(label = 'TE family')





# QC Figures

## Suppl. Figure 2A - TBlastn hits on Bgh, Enecator, Q. suber

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

ggplot(sce)+
  geom_point(aes(`Size (bp)`,`Coverage`, colour =`ID (%)`), alpha=0.5)+
  scale_colour_gradient2( low = 'blue',mid = 'white', high = 'red', midpoint = 50)+
  scale_x_log10(breaks = c(0,1000,10000,100000,1000000))+
  scale_y_log10(breaks = c(seq(0,100,10),1,5,200,500))+
  facet_wrap(~species)+
  theme_minimal()

## Suppl. Figure 3A - Density plots for the 10 largest contigs

##for the TE joy plots

genes <- read.table('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/res.per_contigs.high_confidence.lst.gff.bed_top10', sep = '\t',header = F)
tes <- read.table('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/res.per_contigs.high_confidence.lst.fa.out.gff.TE.bed_top10', sep = '\t',header = F)
sps <- read.table('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/res.per_contigs.high_confidence.lst.gff.bed_top10_secreted', sep = '\t',header = F)
genes <- cbind(genes, 'Genes')
tes <- cbind(tes, 'TEs')
sps <- cbind(sps, 'SPs')

colnames(genes) <- c('scaffold','start','stop','id','type')
colnames(tes) <- c('scaffold','start','stop','id','type')
colnames(sps) <- c('scaffold','start','stop','id','type')


dat <- rbind(genes,tes,sps)

ggplot(dat)+ 
  geom_joy(aes(x =start , y = type) ,bandwidth = 5000)+
  labs(title = 'Gene and repetitive element density',
       subtitle = 'Density of elements per 5kb in the ten largest scaffolds') +
  theme_joy(font_size = 8, grid = T) + theme(axis.title.y = element_blank()) +
  scale_x_continuous()+
  facet_wrap(~ scaffold, ncol = 2, scales = 'free_x')













#### OTHER OLDER PART ####


table1 <- read.csv('non-bact.contigs.lst.histcov', sep = '_',header = F)
hist(table1$V1, breaks = 1000, xlim = c(0,10000))
colnames(table1) <- c('Length','Depth')
plot(table1$V1,table1$V2, ylim = c(0,100), xlim = c(0,1000))
table2 <- read.csv('maker_prot_hist.txt', header = F)
hist(table2$V1,xlim = c(0,50),breaks = 1000)


table1.s <- subset(table1, Depth > 25)
interg <- read.csv('intergenic_spaces.filt.bed', sep = '\t', header = F)
hist(interg$V10, breaks = 100, xlim = c(0,20000))

ggplot(interg)+
  geom_histogram(aes(x=interg$V10))+
  facet_wrap(~interg$V1)
  
#parau genedensity

ds <- read.csv('closest.d.50k.txt', sep = '\t', header = F)
us <- read.csv('closest.u.50k.txt', sep = '\t', header = F)

ds <- cbind(ds,'DS')
us <- cbind(us,'US')
sps <- merge(ds,us, by = 'V2')
sps$V3.x <- abs(sps$V3.x)
sps$V3.y <- abs(sps$V3.y)

ggplot(sps, aes(x= V3.x,y=V3.y)) +
  geom_hex()+
  scale_fill_distiller(palette = "Spectral", name="Gene\ncount")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("5' prime intergenic length (bp)") +
  xlab("3' prime intergenic length (bp)") +
  theme_minimal()

bghcont <- read.csv('cov_of_bgh.txt')
hist(bghcont$mummer.bgh, breaks = 100000, xlim = c(0,100))  
