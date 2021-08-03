#!bin/bash
# cut adapter


source /home1/songlt/miniconda3/bin/activate


mkdir -p /home1/GENE_proc/SONGLITING/miRNA/clean_fq

clean_fq=/home1/GENE_proc/SONGLITING/miRNA/clean_fq
raw_fq=/home2/Cohorts/Genes/NeuroDegenerative/RNA/Project_s905s09001_244Samples_20200430


cat /home1/GENE_proc/SONGLITING/miRNA/samples_info/samples_info.txt |sed -n '1,40p' | cut -f1 | while read sample

do
ls ${raw_fq}/*/*/*${sample}*_R1.fastq.gz 

/home1/songlt/miniconda3/bin/trim_galore --small_rna  -q 20 --phred33 --stringency 3 -j 4 --length 18 -e 0.1 \
            --paired ${raw_fq}/*/*${sample}*_R1.fastq.gz  ${raw_fq}/*/*${sample}*_R2.fastq.gz  \
            -o ${clean_fq}


done 
 
