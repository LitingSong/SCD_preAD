#!bin/bash

ref_gff=/home1/songlt/ref_genome/genome/hg19/gencode.v33lift37.annotation.gtf
ref_fa=/home1/songlt/ref_genome/genome/hg19/hg19.fa
clean_fq=/home1/GENE_proc/SONGLITING/mRNA/clean_fq
STAR=/home1/songlt/.local/bin/STAR

SJ=`ls /home1/GENE_proc/SONGLITING/mRNA/star_1pass/*SJ.out.tab`

#2. star 1-pass align
mkdir -p /home1/GENE_proc/SONGLITING/mRNA/star_2pass


cat /home1/GENE_proc/SONGLITING/mRNA/samples_info/samples_info.txt |sed -n '1,3p' | cut -f1 | while read sample


do
$STAR --runThreadN 20 \
  --genomeDir /home1/songlt/ref_genome/index/STAR/ \
  --readFilesCommand zcat \
  --readFilesIn ${clean_fq}/${sample}.R1.fastq.gz   ${clean_fq}/${sample}.R2.fastq.gz \
  --sjdbFileChrStartEnd $SJ \
  --outReadsUnmapped Fastx \
  --outSAMtype BAM SortedByCoordinate \
  --outFilterType BySJout \
  --quantMode TranscriptomeSAM \
  --outFilterMultimapNmax 20 \
  --alignSJoverhangMin 8 \
  --alignSJDBoverhangMin 1 \
  --outFilterMismatchNmax 999 \
  --outFilterMismatchNoverReadLmax 0.04 \
  --alignIntronMin 20  \
  --alignIntronMax 1000000 \
  --limitSjdbInsertNsj 2000000 \
  --alignMatesGapMax 1000000 \
  --outFileNamePrefix ${sample}.sort


done 

