#!bin/bash
# cut adapter


source /home1/songlt/miniconda3/bin/activate
mkdir -p /home1/GENE_proc/SONGLITING/mRNA/clean_fq

clean_fq=/home1/GENE_proc/SONGLITING/mRNA/clean_fq
raw_fq=/home2/Cohorts/Genes/NeuroDegenerative/RNA/Project_s905r10001_244Samples_20200506


cat /home1/GENE_proc/SONGLITING/mRNA/samples_info/samples_info.txt  | cut -f1 | while read sample

do
ls ${raw_fq}/*/*/*${sample}*_R1.fastq.gz 

/home1/songlt/miniconda3/bin/trim_galore -j 4 -q 20 --phred33 --stringency 3 --length 50 -e 0.1 \
            --paired ${raw_fq}/*/*/*${sample}*_R1.fastq.gz  ${raw_fq}/*/*/*${sample}*_R2.fastq.gz  \
            -o ${clean_fq}


done 
 
