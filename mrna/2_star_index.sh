#!bin/bash

#source ~/miniconda3/bin/activate
#softw=~/miniconda3/bin

STAR=/home1/songlt/.local/bin/STAR

ref_gff=/home1/songlt/ref_genome/genome/hg19/gencode.v33lift37.annotation.gtf
ref_fa=/home1/songlt/ref_genome/genome/hg19/hg19.fa

# 1. star genome  index

$STAR --runThreadN 20 --runMode genomeGenerate \
   --genomeFastaFiles ${ref_fa} \
   --sjdbGTFfile ${ref_gff}  \
   --genomeDir /home1/songlt/ref_genome/index/STAR
