#!bin/bash

source /home1/songlt/miniconda3/bin/activate
softw=/home1/songlt/miniconda3/bin
ref_fa=/home1/songlt/ref_genome/genome/hg19/hg19.fa

dedup_fq=/home1/GENE_proc/SONGLITING/miRNA/dedup_fq
raw_fq=/home2/Cohorts/Genes/NeuroDegenerative/RNA/Project_s905s09001_244Samples_20200430

cd ${dedup_fq}

cat /home1/GENE_proc/SONGLITING/miRNA/samples_info/samples_info.txt | sed -n '1p' |cut -f1 | while read sample

do

# s1. cut the adapter behind the umi
${softw}/cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  --discard-untrimmed   -e 0.1 -O 5 --minimum-length=18   --cores=0 -q 30 -o ${sample}.s1.fastq.gz  ${raw_fq}/*/*${sample}*_R1.fastq.gz

# s2: extract umi
${softw}/umi_tools extract --stdin=${sample}.s1.fastq.gz  --extract-method=string --bc-pattern=NNNNNNNNNNNN   --3prime  --stdout ${sample}.s2.fastq.gz


# s3. cut the adapter ahead of the umi
${softw}/cutadapt -a AACTGTAGGCACCATCAAT  --discard-untrimmed   -e 0.1 -O 5 --minimum-length=18  --maximum-length=24   --cores=0 -q 30 -o  ${sample}.R1.fastq.gz ${sample}.s2.fastq.gz

# s4. bowtie

# 4.1 construct reference index
####$softw/bowtie-build ${ref_fa} --threads 5 ~/ref_genome/index/bowtie/hg19

# 4.2 biwtie align
gunzip ${sample}.R1.fastq.gz
$softw/bowtie --threads 8 -v 2 -m 10 -a --norc /home1/songlt/ref_genome/index/bowtie/hg19  ${sample}.R1.fastq  --sam > ${sample}.sam

/home1/songlt/.local/bin/samtools import ${ref_fa}  ${sample}.sam  ${sample}.bam 

/home1/songlt/.local/bin/samtools sort ${sample}.bam  -o ${sample}.sort.bam

/home1/songlt/.local/bin/samtools index ${sample}.sort.bam

# 5 dedup
$softw/umi_tools dedup -I ${sample}.sort.bam  -S ${sample}.dedup.bam

# 6 sam to fastq
/home1/songlt/.local/bin/samtools view ${sample}.dedup.bam |  grep -v ^@ | awk 'NR%2==1 {print "@"$1"_1\n"$10"\n+\n"$11}' > ${sample}.dedup.fastq


done
