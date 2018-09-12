##CAZY
library(ggplot2)
library(ggthemr)
library(gplots)
library(reshape2)

library(reshape)

caz <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/06_cazy/all_results_for_R_summarized_final.csv', sep = ',', header = F)
colnames(caz) <- c('CAZyme','Count','Species')


caz$Count <- as.numeric(caz$Count)
ggthemr('dust')
ggthemr('greyscale')



ggplot(caz, aes(caz$Species,caz$CAZyme))+
  geom_tile(aes(fill = caz$Count))+
  geom_text(aes(label = caz$Count ))+
  #scale_fill_continuous(name = "Gene count")+
  scale_fill_gradient(low = "white", high = "red",name = "Gene count")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1), plot.title = element_text(hjust=0.5))+
  xlab('')+
  ylab('CAZyme group')
  
caz.m <- acast(caz,CAZyme~Species, value.var = 'Count')



heatmap.2(caz.m[,-13],
          margins =c(35,5),
          trace="none",
          col = colorpanel(10, low = 'snow',mid = 'indianred', high = 'darkred'),
          dendrogram="column",
          denscol = 'black',
          cexRow = 0.75,
          cexCol = 1,
          key = T,
          key.title = '',
          keysize = 1,
          key.xlab = 'Count',
          key.ylab = '',
          density.info = 'none',
          scale = 'none',
          cellnote=caz.m,
          notecol='black'
                    )
