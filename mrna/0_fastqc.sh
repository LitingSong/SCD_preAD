#!bin/bash

RNA_fq_s=/home2/Cohorts/Genes/NeuroDegenerative/RNA/Project_s905s09001_244Samples_20200430
RNA_fq_r=/home2/Cohorts/Genes/NeuroDegenerative/RNA/Project_s905r10001_244Samples_20200506



mkdir -p Project_s905s09001_244Samples_20200430 
mkdir -p Project_s905r10001_244Samples_20200506 


#ls ${RNA_fq_r}/*/*/*fastq.gz | while read id;do (nohup /home1/songlt/miniconda3/bin/fastqc  -q  -o  Project_s905r10001_244Samples_20200506  $id &);done

ls ./Project_s905s09001_244Samples_20200430/*html|cut -d "_" -f4,5,6 |cut -d '/' -f2 > qc_file

ls ${RNA_fq_s}/*/*fastq.gz > mirna_file

grep -vFf  qc_file mirna_file > dif_mi

head -n 60 dif_mi  | while read id;do (nohup /home1/songlt/miniconda3/bin/fastqc  -t 4   -o  Project_s905s09001_244Samples_20200430  $id &);done


#ls ${RNA_fq_r}/*/*/*fastq.gz > mrna_file
#ls ./Project_s905r10001_244Samples_20200506/*html|cut -d "_" -f4,5,6 |cut -d '/' -f2 > qc_file_mrna

#grep -vFf  qc_file_mrna mrna_file > dif_m


#cat dif_m | while read id;do (nohup /home1/songlt/miniconda3/bin/fastqc  -q  -o  Project_s905r10001_244Samples_20200506  $id &);done




