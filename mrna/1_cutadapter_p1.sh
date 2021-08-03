#!bin/bash
# cut adapter



mkdir -p /home1/GENE_proc/SONGLITING/mRNA/clean_fq

clean_fq=/home1/GENE_proc/SONGLITING/mRNA/clean_fq
raw_fq=/home2/Cohorts/Genes/NeuroDegenerative/RNA/Project_s905r10001_244Samples_20200506


cat /home1/GENE_proc/SONGLITING/mRNA/samples_info/samples_info.txt  | cut -f1 | while read sample

do
ls ${raw_fq}/*/*/*${sample}*_R1.fastq.gz 
/home1/songlt/miniconda3/bin/cutadapt \
	-a AGATCGGAAGAGC \
	-A AGATCGGAAGAGC \
	-e 0.1 -O 5 -m 50  -n  2  --cores=0 -q 30 \
	-o ${clean_fq}/${sample}.R1.fastq.gz -p ${clean_fq}/${sample}.R2.fastq.gz  ${raw_fq}/*/*/*${sample}*_R1.fastq.gz  ${raw_fq}/*/*/*${sample}*_R2.fastq.gz


done 
 
