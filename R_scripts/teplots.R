##for the TE plots
setwd("~/Documents/GitHub/jap_pm_2018")
library(ggplot2)
library(ggjoy)


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

##for the kamoun plots

ds <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/down.txt', sep = '\t', header = F)
us <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/up.txt', sep = '\t', header = F)
sps <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_TE_distribution/m.res.per_contigs.high_confidence.gene.lst.fa.tmhmm.notmm', sep = '\t', header = F)

ds <- cbind(ds,'DS')
us <- cbind(us,'US')
all <- merge(ds,us, by = 'V4')
all$V9.x <- abs(all$V9.x)
all$V9.y <- abs(all$V9.y)

all$group <- "Non SP"
all$group[all$V4 %in% sps$V1] <- "SP"



ggplot(all, aes(x= V9.x,y=V9.y)) +
  #geom_hex()+
  stat_bin2d(binwidth=c(0.06, 0.06))+
  scale_fill_distiller(palette = "Spectral", name="Gene\ncount")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("5' prime intergenic length (bp)") +
  xlab("3' prime intergenic length (bp)") +
  geom_point(data=subset(all, group == "SP"), color = 'purple',alpha = 0.5)+
  theme_minimal()



#parau distance from TE


hist(ds$V9)

hist(us$V9)



