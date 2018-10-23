#!/bin/bash

#BSUB -J iqtree-fast
#BSUB -u lamprinos.frantzeskakis@rwth-aachen.de
#BSUB -o job.%J.%I
#BSUB -W 72:00
#BSUB -M 1800
#BSUB -R "select[hpcwork]"
#BSUB -n 32
#BSUB -a openmpi
#BSUB -N 1

source /home/lf216591/.bash_profile

cd /work/lf216591/07_pleo_para_annot/01_parauncinula/09_rnases/04_rnase_tree

/home/lf216591/utils/iqtree-1.6.beta4-Linux/bin/iqtree -alrt 1000 -bb 1000 -s all_secreted_with_no_pfam.fa.aln -nt 24

