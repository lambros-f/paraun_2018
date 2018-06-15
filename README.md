# paraun_2018
Parauncinula genomics 2018

# Bacterial contamination

bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.out | sort | uniq | awk 'FS="_"{num+=$6;++i}END{print num" "i" "num/i}'
65750.2 50270 1.30794
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.out | sort | uniq | awk 'FS="_"{num+=$4;++i}END{print num" "i" "num/i}'
64098655 50270 1275.09
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.noqcov.out | sort | uniq | awk 'FS="_"{num+=$4;++i}END{print num" "i" "num/i}'
129098971 60377 2138.21
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.noqcov.out | sort | uniq | awk 'FS="_"{num+=$6;++i}END{print num" "i" "num/i}'
86292.3 60377 1.42922

Contigs found, q cov is 25%
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.noqcov.out | sort | uniq | wc -l
60377
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.out | sort | uniq | wc -l
50270

Max size max cov with no qcov
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.noqcov.out | cut -f4,6 -d '_' | sort -nr -k1 -t '_' | head
1808665_34.602153
1608424_34.293083
1288287_33.830732
1091948_33.694573
955146_34.596150
853683_33.525800
836785_33.523844
683076_34.789354
635506_33.547335
615394_35.365817
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.noqcov.out | cut -f4,6 -d '_' | sort -nr -k2 -t '_' | head
1115_2748.177291
1823_2384.869159
716_1629.723967
660_1588.573770
1586_1538.067119
8080_956.982808
6430_613.553885
8477_588.458881
8279_579.071866
25515_568.281845

Max size max cov with qcov
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.out | cut -f4,6 -d '_' | sort -nr -k1 -t '_' | head
67938_3.481033
32314_2.513368
30235_2.779976
27808_3.879843
27178_2.455684
26140_3.466825
25626_2.605605
24414_2.285603
24259_5.762092
24227_3.433820
bash-4.2$ cut -f1 scaffolds.parau.draft2.fasta.bact.blast.out | cut -f4,6 -d '_' | sort -nr -k2 -t '_' | head
1115_2748.177291
1823_2384.869159
716_1629.723967
660_1588.573770
1586_1538.067119
2611_294.599200
511_186.402500
1543_141.188547
6042_104.212612
697_93.957338
