# Parauncinula septata

Information about the sample:

Parauncinula septata on Quercus serrata
collected: 9 November 2017 in Torimiyama Park, Haibara, Uda-shi, Nara Prefecture, Japan
(DNA extracts were prepared from the same collection)

DNA extraction was performed either by collecting chasmothecia or by acetone/cellulose peeling of the collected leaves. 

DNA was send for sequencing at CEGAT and was sequenced with Illumina HiSeq?, 2x150bp.

This is for later:
Pleochaeta shiraiana on Celtis sinensis
collected: 13 December 2017 in Mie University Campus, Tsu-shi, Mie Prefecture, Japan
(DNA extracts were prepared from powdery mildew on the same host plant individual on earlier dates)
 
Please be careful with Celtis leaves as they can be co-infected with Erysiphe kusanoi besides P. shiraiana (but you can easily differentiate them by chasmothecial morphology and size).
 
I think most probably most of the extracts in tubes contain very tiny amount of DNA (if any at all).
 
## Assembly of the reads

See original script I submitted. This below is indicative of that I run.

```
# First I run bfc to reduce the complexity. I kept only single ended reads, gzip results

bfc -b 32 -k 25 -t 10 out.fastq.gz 2>| input.corr.fastq.gz.bfc.e | seqtk dropse - 2>| input.corr.fastq.gz.seqtk.e | pigz -c - -p 4 -2  1>| input.corr.fastq.gz 2>| input.corr.fastq.gz.pigz.e

# This I run spades, we normally presplit the input  since it is faster.

spades.py \
-m 2000 \
--tmp-dir $SLURM_TMP \
-o spades3 \
--only-assembler \
-k 33,55,77,99,127 \
--meta \
-t 32 \
-1 /global/projectb/scratch/aclum/metagenomes/benchmarking/bfc_params/mock/split/reads1.fasta \
-2 /global/projectb/scratch/aclum/metagenomes/benchmarking/bfc_params/mock/split/reads2.fasta

```

At the end you get 658405 contigs, Total 481 MB, Max 1,8 MB, Min 112bp , N50 of 1,1kb. Not good, but ok.

After removing contigs less than 500bp you get 159829, Total 316 MB, Max 1,8 MB, Min 501bp, N50 of 4,1kb. So we remove a lot of the small contigs with not so much impact in the total sequence content. The small ones are not useful for the downstream analysis anyway.

Ploting contig size and coverage shows that there no correlation, but that there is a population of contigs ~30x coverage and another with less than 2x. This 30x contigs are mostly Parauncinula, which is nice.

## Cleaning up bacterial contaminating sequences

Since the sample is metagenomic, it makes sense to remove as much bacterial contigs as possible. I blasted the contigs to 3837 bacterial genomes that associated with plants which can be found here (http://labs.bio.unc.edu/Dangl/Resources/gfobap_website/index.html , see also Asaf Levy Nat. Genetics paper). 

For the blasting I used:
```
/home/lf216591/utils/ncbi-blast-2.5.0+/bin/blastn \
-query draft_2_meta/scaffolds.draft2.500bp.fasta \
-db bact_genomes/genomes_fna.fa \
-max_hsps 1 \
-max_target_seqs 1 \
-evalue 10e-5 \
-outfmt 6 \
-qcov_hsp_perc 25 \
-num_threads 12 > scaffolds.parau.draft2.fasta.bact.blast.out
```
I also used no qcov filtering, but the qcov can help in keeping contigs with very short hits to bacteria. Results are here:

Stats of removed sequences | With Qcov | No Qcov
----| ----| ----|
Total number | 50270 | 60377
Total length | 64098655 | 129098971
Average length | 1275.09 | 2138.21
Average coverage | 1.30794 | 1.42922
Max size | 67938 | 1808665
Max cov | 2748.177291 | 2748.177291

The cleaned-up assembly is:

Scaffold Stats | Count/Size
----| ----
Seqs | 109,559
Min | 480 
1st Qu. | 639 
Median | 845 
Mean | 2,321 
3rd Qu. | 1,335 
Max | 1,808,665 
Total | 254,390,543
n50 | 7,978 
n90 | 743 
n95 | 611


## Cleaning up fungal and plant contaminating sequences

I started by doing a tblastn of the filtered contigs to the proteomes of Bgh, E. necator and Quercus suber (related to the host plant).

Then I parsed the results in R and generated a plot of size and coverage of contigs with hits from each of the three species.

The figure shows that the majority of the plant hits are in the contigs with coverage lower than 20x. 
For the fungal hits it is very not clear, whether they are predominantly on the lower 20 bracket or the higher 20 bracket. But considering that fungal hits can be on the plant genes in the lower than 20x, I would consider that making a cut-off for contigs lower than 20x.

After making this cut-off the assembly (scaffolds.draft2.500bp.nobact.d20up.fa) looks like:

Scaffold Stats | Count/Size
----| ----
Seqs  |      1,321
Min    |        501 
1st Qu.|      1,774 
Median |      9,327 
Mean   |     45,360 
3rd Qu.|     39,135 
Max    |  1,808,665 
Total  | 59,921,213 
n50    |    176,477 
n90    |     25,733 
n95    |     15,552 

There quite some loss of contigs but not as much sequence.

To verify if this is dataset is meaningful I run BUSCO on the predicted proteome after a MAKER annotation. The MAKER annotation was generated with the Bgh pipeline.

After the annotation I also changed the headers of the fasta and gff to PARAU_XXXXX . These files are: paraud2.maker.gff,paraud2.maker.prots.fa etc and correspond to the scaffolds.draft2.500bp.nobact.d20up.fa assembly.

C:95.4%,S:33.9%,D:61.5%,,F:1.7%,M:2.9%,n:1438

Number|Category
---|---
1372|	Complete BUSCOs (C)
488|	Complete and single-copy BUSCOs (S)
884|	Complete and duplicated BUSCOs (D)
24|	Fragmented BUSCOs (F)
42|	Missing BUSCOs (M)
1438|	Total BUSCO groups searched

Although the results look good, the issue here is the high number of duplicated BUSCOs indicating that this dataset has multiple fungal genomes

See for example the Bgh result:

Number|Category
---|---
1403|	Complete BUSCOs
1266|	Complete and single-copy BUSCOs
137|	Complete and duplicated BUSCOs
21|	Fragmented BUSCOs
14|	Missing BUSCOs
1438|	Total BUSCO groups searched

In order to reduce the dataset and exclude the other fungal genomes I would like to pull out only the ones that have hits to the Leotiomycetes over the majority of the contigs.

Therefore I blastp the predicted proteome of these contigs to the nr. There is a small change to the typical blastp command since I want to get the species name for each of the hits. You also have to make a alias/link to the fungal part of the nr.

```
blastdb_aliastool -gilist refseq_fungi_prot_db -db /hpcwork/rwth0146/db/nr -out refseq_fungi_prot_db_db -title refseq_fungi_prot_db

blastp \
-query paraud2.maker.prots.fa \
-db refseq_fungi_prot_db_db  \
-max_hsps 2 \
-max_target_seqs 2 \
-evalue 10e-5 \
-num_threads 12 \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue sscinames staxids bitscore" > paraud2.maker.prots.fa.nr.REFSEQFUNGAL.blastp.out
```

This returns something like this:

```
PARAU_12488	gi|1069187622|ref|XP_018066279.1|	57.22	783	274	19	1	769	1	736	0.0	Phialocephala scopiformis	149040	771
PARAU_12488	gi|597585795|ref|XP_007295954.1|	56.09	788	266	20	39	771	5	767	0.0	Marssonina brunnea f. sp. 'multigermtubi' MB_m1	1072389	739
PARAU_12346	gi|1069194092|ref|XP_018063050.1|	57.32	717	274	10	1	689	1	713	0.0	Phialocephala scopiformis	149040	752
PARAU_12346	gi|597583617|ref|XP_007294865.1|	60.33	721	248	16	1	692	1	712	0.0	Marssonina brunnea f. sp. 'multigermtubi' MB_m1	1072389	739
PARAU_12345	gi|597587079|ref|XP_007296596.1|	60.57	601	222	8	7	594	15	613	0.0	Marssonina brunnea f. sp. 'multigermtubi' MB_m1	1072389	744
PARAU_12345	gi|1069194032|ref|XP_018063020.1|	60.70	603	218	6	9	595	4	603	0.0	Phialocephala scopiformis	149040	734
PARAU_09993	gi|164422135|ref|YP_001648749.1|	66.01	253	66	5	32	264	311	563	3e-94	Zymoseptoria tritici	1047171	301
PARAU_09993	gi|1236514537|ref|YP_009409399.1|	69.04	197	60	1	27	223	326	521	3e-87	Arthrinium arundinis	335852	283
PARAU_09990	gi|1011058802|ref|YP_009240970.1|	68.54	267	32	3	139	353	1	267	2e-106	Cairneyella variabilis	1802957	325
PARAU_09990	gi|1236517182|ref|YP_009412632.1|	66.42	268	38	3	139	354	1	268	6e-106	Imshaugia aleurites	172621	324
```

Once this is over (~ couple of days), I need to parse the results. Assuming that most of the hits in a correctly assembled contig will come from a single species, I want to get the first, second hit and their frequency along the hits on this contig.

```
#make_blastp_results_table.sh
rm hits.per_contigs.txt
rm res.per_contigs.txt
echo 'Contig Num_Genes Num_Blastp_hits 1st_Frq_hit 2nd_Frq_hit 1st_Frq_count 2nd_Frq_count 1st_Frq_% 2nd_Frq_%' > res.per_contigs.txt
for contig in `cat contig.lst`;do
#echo $contig
grep $contig paraud2.maker.gff | grep -P 'gene\t' | grep -o PARAU_..... | sort | uniq > $contig.genes.lst

	for gene in `cat $contig.genes.lst`;do
#	echo $gene
#	genehit1=$(grep $gene paraud2.maker.prots.fa.nr.REFSEQFUNGAL.blastp.out | head -n 1 | cut -f1 )
	genehit2=$(grep $gene paraud2.maker.prots.fa.nr.REFSEQFUNGAL.blastp.out | sort -nrk3 |  head -n 1 |cut -f12 )
		if [ -z "$genehit2" ]
		then genehit2="N/A"
		fi
	echo $contig"	"$gene"	"$genehit2 >> hits.per_contigs.txt
	done

#echo 'done'
numgenespc=$(wc -l $contig.genes.lst|cut -f1 -d ' ' )
numhitspc=$(grep $contig hits.per_contigs.txt |grep PARAU | grep -v 'N/A'| wc -l | cut -f1 -d ' ')
fsthitN=$(grep $contig hits.per_contigs.txt | grep PARAU | cut -f3 | sort | uniq -c | sort -nrk1 | sed 's/^ *//'|head -n 1 |cut -f1 -d ' ')
fsthitID=$(grep $contig hits.per_contigs.txt| grep PARAU | cut -f3 | sort | uniq -c | sort -nrk1 | sed 's/^ *//'|head -n 1 |cut -f2 -d ' ')
sndhitN=$(grep $contig hits.per_contigs.txt | grep PARAU |grep -v 'N/A'| cut -f3 | sort | uniq -c | sort -nrk1 | sed 's/^ *//'|head -n 2 | tail -n 1 | cut -f1 -d ' ')
sndhitID=$(grep $contig hits.per_contigs.txt| grep PARAU |grep -v 'N/A'| cut -f3 | sort | uniq -c | sort -nrk1 | sed 's/^ *//'|head -n 2 | tail -n 1 | cut -f2 -d ' ')

perc_fst=$(echo $fsthitN/$numgenespc | bc -l)
perc_snd=$(echo $sndhitN/$numgenespc | bc -l)
echo $contig' '$numgenespc' '$numhitspc' '$fsthitID' '$sndhitID' '$fsthitN' '$sndhitN' '$perc_fst' '$perc_snd >> res.per_contigs.txt

done

#High confidence contigs

## here it pulls out the Leotiomycetes species from frequent_species.orders.leot.lst

awk 'NR==FNR{a[$1];next}$4 in a' frequent_species.orders.leot.lst res.per_contigs.txt >1 && awk 'NR==FNR{a[$1];next}$5 in a' frequent_species.orders.leot.lst 1 |column -t > res.per_contigs.high_confidence.txt && rm 1

grep -f res.per_contigs.high_confidence.lst paraud2.maker.gff | grep -P 'gene\t' | grep -o 'PARAU_.....' | sort | uniq > res.per_contigs.high_confidence.gene.lst

sh get_seqs.sh res.per_contigs.high_confidence.gene.lst paraud2.maker.prots.fa res.per_contigs.high_confidence.gene.lst.fa

#Low confidence/Chimeric contigs




grep -f frequent_species.orders.leot.lst res.per_contigs.txt |  awk '{if ($8 + $9 > .5) print }' | column -t | sort -nrk2 | grep -v  'N/A' > res.per_contigs.high_confidence.txt
cut -d ' '  -f1 res.per_contigs.high_confidence.txt  > res.per_contigs.high_confidence.lst


#stats
grep -f frequent_species.orders.leot.lst res.per_contigs.txt |  awk '{if ($8 + $9 > .5) print }' | grep -v  'N/A' | awk 'FS="_"{num+=$4}END{print num}'
grep -f frequent_species.orders.leot.lst res.per_contigs.txt |  \
awk '{if ($8 + $9 > .5) print }' | grep -v  'N/A' | cut -f2 -d ' ' | awk '{num+=$1}END{print num}'
```

This prints a table like this (see res.per_contigs.txt) :

Contig|Num_Genes|Num_Blastp_hits|1st_Frq_hit|2nd_Frq_hit|1st_Frq_count|2nd_Frq_count|1st_Frq_%|2nd_Frq_%
---|---|---|---|---|---|---|---|---|
NODE_1000_length_27779_cov_37.567696 |3 |3 |Phialocephala |Sclerotinia |2 |1 |.66666666666666666666 |.33333333333333333333
NODE_1001_length_27771_cov_29.788901 |4 |4 |Phialocephala |Sclerotinia |2 |1 |.50000000000000000000 |.25000000000000000000
NODE_1030_length_26970_cov_36.415205 |4 |4 |Phialocephala |Marssonina |1 |1 |.25000000000000000000 |.25000000000000000000

Based on the lineages provided by TaxonomyDB (lineages-2018-06-13.csv) you can find out if a genus belongs to the leotiomycetes or not.

But also the sequences for the high confidence contigs (first and and second most frequent hit is Leotiomycetes) and everything else as low confidence contigs. This is a quite strict filter but you can avoid potentially chimeric contigs.

The high and low confidence assemblies look like this:


Scaffold Stats | High confidence | Low confidence
----| ---- | ----
Seqs  | 495|107
Min    | 966| 3204   
1st Qu.| 17297|    
Median | 36067|
Mean   | 56598|
3rd Qu.| 74488|
Max    | 354858| 1808665
Total  | 28016241| 10976294
n50    | 98897|
n90    | 25733|
n95    | 17861|

The rest of the contigs have no fungal hits at all so they are excluded. Of course this is also removing other contigs where genes were not predicted by maker.

The BUSCO results look better, in a sense that the duplicated BUSCOs are significantly reduced:

Number|Category
---|---
1305|Complete BUSCOs
1172|Complete and single-copy BUSCOs
133|Complete and duplicated BUSCOs
57|Fragmented BUSCOs
76|Missing BUSCOs
1438|Total BUSCO groups searched

The completeness is 90% vs 97% of Bgh, which I would say is as high as you can go with this dataset.

So now the final assembly is "res.per_contigs.high_confidence.lst.fa", and based on this result I think we can move on with the analysis.



## Analysis

Note that in the scripts used and the output files the file res.per_contigs.high_confidence.lst.fa is the genome assembly with the high confidence contigs of P. polyspora.

### Genomes/Proteomes used for analysis

The Leotiomycete and Ascomycete datasets used for the downstream analysis see Suppl. Table 8.

### MAKER

For the annotation, the setting files (*.ctl) used for maker are in various_scripts.

### Orthofinder

Orthofinder version used is 1.1.2, and was called with default options e.g. "orthofinder -t 12 -f protein_seq_folder"

The output from the comparison of the 16 Leotiomycete species (proteomes.tar.gz), including the derived species tree is in useful_files/01_Orthofinder_output 

### RepeatMasker

RepeatMasker version used is 4.0.6, and was called by:
```
RepeatMasker -species fungi -pa 12 -excln -gff -no_is genome.fa
```

Output is in useful_files/02_Repeatmasker

### Secretome

For the secretome SignalP v4.1 and TMHMM v2.0c, were used. The output from SignalP was passed to TMHMM to screen for transmembrane domains (tm) and then proteins with tm were removed. The script used was:

```
file1=res.per_contigs.high_confidence.gene.lst.fa

/home/lf216591/utils/signalp-4.1/signalp \
-m m.$file1 \
$file1 > $file1.signalp.out

/home/lf216591/utils/tmhmm-2.0c/bin/tmhmm \
m.$file1 > m.$file1.tmhmm

grep 'Number of predicted TMHs:  0' m.$file1.tmhmm | cut -f2 -d ' ' > m.$file1.tmhmm.notmm

sh get_seqs.sh m.$file1.tmhmm.notmm m.$file1 m.$file1.tmhmm.notmm.fa
```

The get_seqs.sh script is a simple bash script to get sequences from a list of headers. Copied it in various_scripts.




