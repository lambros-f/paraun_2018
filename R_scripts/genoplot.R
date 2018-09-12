setwd("~/Documents/GitHub/jap_pm_2018")
library(ggplot2)
library(genoPlotR)

# Example 1
data(three_genes)

comparisons[[1]]$col <- apply_color_scheme(c(0.6, 0.4, 0.5), "grey")

names <- c("Huey", "Dewey", "Louie")

names(dna_segs) <- names

tree <- newick2phylog("(((Huey:4.2,Dewey:3.9):3.1,Louie:7.3):1);")

mid_pos <- middle(dna_segs[[1]])

xlims <- list(c(Inf, -Inf), c(-Inf, Inf), c(1850, 2800))

annot <- annotation(x1=c(mid_pos[1],dna_segs[[1]]$end[2]),
                    x2=c(NA, dna_segs[[1]]$end[3]),
                    text=c(dna_segs[[1]]$name[1], "region1"),
                    rot=c(30, 0), col=c("blue", "black"))

plot_gene_map(dna_segs=dna_segs, comparisons=comparisons,
                annotations=annot, annotation_height=1.3,
                tree=tree, tree_width=2,
                xlims=xlims,
                main="Comparison of Huey, Dewey and Louie")

# Real example

#comparison1 <- read_comparison_from_blast('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/blastn2.out')
#comparison2 <- read_comparison_from_blast('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/blastn.out')
#comparison2b <- read_comparison_from_tab('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/node35-scaf_23.delta.tsv', header = F, gene_type = 'lines')
comparison2b <- read_comparison_from_tab('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/node35-scaf_23.delta.tsv2', header = F, gene_type = 'lines')
comparison3 <- list(comparison2b)

#dna_segs2 <- read_dna_seg_from_fasta('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/scaffold_23.fa')
dna_segs2b <- read_dna_seg_from_embl('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/scaffold_23.fa.embl', tagsToParse=c("gene"))
dna_segs3 <- read_dna_seg_from_embl('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/node35.fa.embl')

annot2 <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/node35.fa.annot', sep = '\t', header = F)
annot2$V3 <- as.character(annot2$V3)
annot2.a <- annotation(annot2$V1,annot2$V2,annot2$V3)
annot2.a <- annotation(annot2$V1,annot2$V2,annot2$V3,col = "red",rot = 90)

annot3 <- read.csv('/Users/lfrantzeskakis/projects/03_jp_mildews/simple_synteny/scaffold_23.fa.annot', sep = '\t', header = F)
annot3$V3 <- as.character(annot3$V3)
annot3.a <- annotation(annot3$V1,annot3$V2,annot3$V3)
annot3.a <- annotation(annot3$V1,annot3$V2,annot3$V3,col = "black",rot = 90)

dna_segs.l <- list(dna_segs3,dna_segs2b)
annot.l <- list(annot2.a,annot3.a)

#gene_types

plot_gene_map(dna_segs = dna_segs.l, 
              comparisons = comparison3,
              annotations = annot.l,
              annotation_cex = 0.1,
              #annotation_height = 1,
              #dna_seg_scale = T,
              #scale = F,
              #gene_type="side_blocks",
              #seg_plot_height = 5,
              #xlims = list(c(100000,250000),c(3500000,4550000)),
              dna_seg_labels = c('node_35','scaffold_23'),
              #gene_type = 'lines',
              global_color_scheme = c('auto','auto','grey',1)
              )
