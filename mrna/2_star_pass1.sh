#!bin/bash

ref_gff=/home1/songlt/ref_genome/genome/hg19/gencode.v33lift37.annotation.gtf
ref_fa=/home1/songlt/ref_genome/genome/hg19/hg19.fa
clean_fq=/home1/GENE_proc/SONGLITING/mRNA/clean_fq
STAR=/home1/songlt/.local/bin/STAR

#2. star 1-pass align
mkdir -p /home1/GENE_proc/SONGLITING/mRNA/star_1pass

cat /home1/GENE_proc/SONGLITING/mRNA/samples_info/samples_info_globin.txt | cut -f1 | while read sample


do
$STAR --runThreadN 20 \
  --genomeDir /home1/songlt/ref_genome/index/STAR/ \
  --readFilesCommand zcat \
  --readFilesIn ${clean_fq}/${sample}.gb.R1.fastq.gz   ${clean_fq}/${sample}.gb.R2.fastq.gz \
  --outFileNamePrefix /home1/GENE_proc/SONGLITING/mRNA/star_1pass/${sample}.gb


done 

