library(ggplot2)
library(reshape2)
library(ggthemr)

ggthemr('dust')
### Bimodal distribution of GC/AT ###

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
#  theme_minimal()+
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
#  theme_minimal()+
  ylab(label = 'Portion of the genome (%)')+
  xlab(label = 'GC content (%)')

### RIP Indexes and dinucleotide frequencies ###

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
#  theme_minimal()+
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
#  theme_minimal()+
  ylab(label = '(CpA + TpG)/(ApC+ GpT)')+
  xlab(label = 'TE family')


#Index 3 - from RIPCAL (shows a depletion of TA?)

ripidx.m <- melt(ripidx)
ripidx.s <- cbind.data.frame(ripidx$species,ripidx$family,ripidx$`(CpA+TpG)/TpA`)
colnames(ripidx.s) <- c('species','family','I3')

ripidx.s1 <- subset.data.frame(ripidx.s, species == 'parauncinula_septata')
ripidx.s2 <- subset.data.frame(ripidx.s, species == 'rhynchosporium_commune')
ripidx.s3 <- subset.data.frame(ripidx.s, species == 'blumeria_graminis_fsp_hordei')
ripidx.s4 <- subset.data.frame(ripidx.s, species == 'marssonina_brunnea')
ripidx.s <- rbind(ripidx.s1,ripidx.s2,ripidx.s3,ripidx.s4)


ripidx.m <- melt(ripidx.s)
ripidx.controls <- subset(ripidx.m, family == 'control')
ripidx.m <- subset(ripidx.m,variable == 'I3')
ripidx.m <- subset(ripidx.m,family == c('Gypsy','Copia','Tad1'))
ripidx.controls <- subset(ripidx.controls, variable == 'I3')

ggplot(ripidx.m)+
  geom_boxplot(aes(x= family , y=value))+
  facet_wrap(~species , scales = 'free', ncol =6)+
  geom_point(aes(x= family , y=value), colour = 'grey', alpha = 0.4)+
  geom_hline(data = ripidx.controls, aes(yintercept = value) , linetype="dashed", colour = 'red')+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5)) +
  theme_minimal()+
  xlab(label = 'TE family')








## SEMI JUNK ##


# This part contains a miscalculation in the raw tables, where averages of summaries are calculated.
#
# ripidx <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/05_ripcal/01_leotio_dinucl_calc/summary.all', header = T, sep = ' ')
# ripidx <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/05_ripcal/01_leotio_dinucl_calc/summary.all.cutoff5count', header = T, sep = ' ')
# 
# colnames(ripidx) <- c('Species','Family','Index1','Index2','Counts')
# ripidx.subset <- subset(ripidx, Family == c('Tad1','Gypsy','Copia'))
# 
# ggplot(ripidx.subset)+
#   geom_point(aes(x= Family, y=Index1 , size = Counts))+
#   facet_wrap(~Species, scales = 'fixed', ncol = 8)+
#   geom_hline(yintercept = 0.89)+
#   geom_point(data=subset(ripidx, Family == 'non_repetitive_control'), aes(x= Family, y=Index1 ), colour="red", size = 6)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 
# ggplot(ripidx.subset)+
#   geom_point(aes(x= Family, y=Index2 , size = Counts))+
#   facet_wrap(~Species, scales = 'fixed', ncol = 8)+
#   geom_hline(yintercept = 1.03)+
#   geom_point(data=subset(ripidx, Family == 'non_repetitive_control'), aes(x= Family, y=Index2 ), colour="red", size = 6)+
#   geom_point(data=subset(ripidx, Family == 'Tad1'), aes(x= Family, y=Index2 ), colour="blue")+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# ### With no averaging - This is also problematic becase subset generates problems in the table
# ripidx <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/05_ripcal/01_leotio_dinucl_calc/summary.all.noavg', header = T, sep = ' ')
# colnames(ripidx) <- c('Species','Family','Index1','Index2','ID')
# 
# ripidx.subset <- subset(ripidx, Family == c('Tad1','Gypsy','Copia'))
# 
# ggplot(ripidx.subset)+
#   geom_violin(aes(x= Family , y=Index2 ))+
#   facet_wrap(~Species, scales = 'fixed', ncol = 8)+
#   geom_hline(yintercept = 1.03)+
#   geom_point(data=subset(ripidx, Family == 'non_repetitive_control'), aes(x= Family, y=Index2 ), colour="red", size = 6)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# 
# ggplot(ripidx.subset)+
#   geom_violin(aes(x= Family , y=Index1 ))+
#   facet_wrap(~Species, scales = 'fixed', ncol = 8)+
#   geom_hline(yintercept = 0.89)+
#   geom_point(data=subset(ripidx, Family == 'non_repetitive_control'), aes(x= Family, y=Index1 ), colour="red", size = 4)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


# ### With no averaging
# # egrep  '\sTad1|Gypsy\s|Copia|Mariner\s|\sHarbinger\s|\sLTR\s|non_repetitive_control' summary.all.noavg > summary.all.noavg.subset
# 
# ripidx <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/05_ripcal/01_leotio_dinucl_calc/summary.all.noavg.subset', header = F, sep = ' ')
# colnames(ripidx) <- c('Species','Family','Index1','Index2','ID')
# 
# #ripidx.subset <- subset(ripidx, Family == c('Tad1','Gypsy','Copia'))
# #keys <- c('Tad1','Gypsy','Copia','Mariner','Harbinger','LTR')
# 
# #ripidx.subset <- grep('Tad1|Gypsy|Copia|Mariner|Harbinger|LTR',ripidx$Family, perl = T)
# 
# #ripidx.subset <- ripidx[ripidx.subset,]
# 
# ggplot(ripidx)+
#   geom_violin(aes(x= Family , y=Index1 ))+
#   facet_wrap(~Species, scales = 'free', ncol = 8)+
#   geom_hline(yintercept = 0.89, linetype="dashed")+
#   geom_point(data=subset(ripidx, Family == 'non_repetitive_control'), aes(x= Family, y=Index1 ), colour="red", size = 3)+
#   #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
#   #scale_y_continuous(limits = c(0, 2))+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5)) 
# 
# ggplot(ripidx)+
#   geom_point(aes(x= Family , y=Index1 ), position = 'jitter')+
#   facet_wrap(~Species, scales = 'fixed', ncol = 8)+
#   geom_hline(yintercept = 0.89, linetype="dashed")+
#   geom_point(data=subset(ripidx, Family == 'non_repetitive_control'), aes(x= Family, y=Index1 ), colour="red", size = 3)+
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# #scale_y_continuous(limits = c(0, 2))+
# 
# 
# 
# ggplot(ripidx)+
#   geom_violin(aes(x= Family , y=Index2 ))+
#   facet_wrap(~Species, scales = 'free', ncol = 8)+
#   geom_hline(yintercept = 1.03, linetype="dashed")+
#   geom_point(data=subset(ripidx, Family == 'non_repetitive_control'), aes(x= Family, y=Index2 ), colour="red", size = 3)+
#   #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
#   #scale_y_continuous(limits = c(0, 2))+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5)) 
# 
