library(ggplot2)

table_par <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/01_interpro_counts/pfam.parau.frequency.csv', sep = ' ', header = F)
table_bgh <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/01_interpro_counts/pfam.bgh.frequency.csv', sep = ' ', header = F)
colnames(table_par) <- c('PFAM','count_par','freq_par')
colnames(table_bgh) <- c('PFAM','count_bgh','freq_bgh')


table_j <- merge(table_par,table_bgh,by = 'PFAM',all = T)
table_j[is.na(table_j)] <- 0
table_j$fold <- table_j$freq_par/table_j$freq_bgh
write.csv(table_j,'/Users/lfrantzeskakis/projects/03_jp_mildews/01_interpro_counts/table_j_forsuppl.txt')

# subseting to > 2 fold

table_j_sub <- subset(table_j, table_j$fold > 2 & table_j$fold < Inf )
table_j_in <- subset(table_j, table_j$fold == Inf )

write.csv(table_j_sub,'/Users/lfrantzeskakis/projects/03_jp_mildews/01_interpro_counts/table_j_sub.txt')
## in bash: grep -f table_j_sub.txt.lst pfam.tsv | awk 'FS="\t" {print $5,$6}' | sort | uniq > table_j_sub.txt.lst.descr
pfam.descr <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/01_interpro_counts/table_j_sub.txt.lst.descr', sep = '\t', header = F)
colnames(pfam.descr) <- c('PFAM','description')
table_j_sub_m <- merge(table_j_sub,pfam.descr, by = 'PFAM',all = T)

ggplot(table_j_sub_m, aes(x=reorder(description,fold), y=fold)) +
  geom_bar(stat='identity', fill = 'black') +
  coord_flip()+
  theme_minimal()


# subseting to > 5 fold to make a more compact figure

table_j_sub_m5 <- subset(table_j_sub_m, table_j_sub_m$fold > 3 )

ggplot(table_j_sub_m5, aes(x=reorder(description,fold), y=fold)) +
  geom_bar(stat='identity', fill = 'black') +
  coord_flip()+
  xlab('')+
  ylab('Fold difference compared to Bgh')+
  theme_minimal()
