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

## Cleaning up bacterial contamination

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

