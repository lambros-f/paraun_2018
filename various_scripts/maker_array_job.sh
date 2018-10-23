#!/bin/bash
 
#BSUB -J makpar
#BSUB -B
#BSUB -P rwth0146
#BSUB -N
#BSUB -o arrayjob.%J.%I
#BSUB -W 00:10
#BSUB -M 2000
#BSUB -R select[hpcwork]

export LIBPATH=/home/lf216591/perl5/lib/perl5
export PERL5LIB=/home/lf216591/perl5/lib/perl5
export LD_PRELOAD=/opt/MPI/openmpi-1.10.2/linux/intel_16.0.2.181/lib/libmpi.so

cd /work/lf216591/07_pleo_para_annot/01_parauncinula/01_paraunc_annotation/01_maker

bsub \
-R select[hpcwork] \
-W 48:00 \
-M 2000 \
-P rwth0146 \
-J "mapar[1-104]" run_maker.sh 

