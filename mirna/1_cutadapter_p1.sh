#!bin/bash
# cut adapter



mkdir -p /home1/GENE_proc/SONGLITING/miRNA/clean_fq

clean_fq=/home1/GENE_proc/SONGLITING/miRNA/clean_fq
raw_fq=/home2/Cohorts/Genes/NeuroDegenerative/RNA/Project_s905s09001_244Samples_20200430


cat /home1/GENE_proc/SONGLITING/miRNA/samples_info/samples_info.txt | sed -n '1,40p' |cut -f1 | while read sample

do

ls ${raw_fq}/*/*${sample}*_R1.fastq.gz

/home1/songlt/miniconda3/bin/cutadapt \
	-a AACTGTAGGCACCATCAAT --discard-untrimmed \
	-e 0.1 -O 5 --minimum-length=18   --cores=0 -q 30 \
	-o ${clean_fq}/${sample}.R1.fastq.gz   ${raw_fq}/*/*${sample}*_R1.fastq.gz    >& log.txt


done 
 
