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

