library(ggplot2)
library(reshape)
options(scipen=999)

rep_table <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/11_repeat_landscape/res.per_contigs.high_confidence.lst.fa.align.landscape.Div.Rclass.tab', sep = '\t', header = T)
colnames(rep_table) <- c('Rclass',1:50)

rep_table.m <- melt(rep_table)


ggplot(rep_table.m,aes(fill=Rclass,x=variable,weight=value))+ 
  geom_bar() +
  theme_minimal()+
  scale_fill_brewer(palette="Set2") +
  xlab("Divergence (%)")+
  ylab("Sequence (bp)") +
  ggtitle("Landscape plots of transposon divergence" )  +
  theme(axis.text.x = element_text(angle =90, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5)) 
