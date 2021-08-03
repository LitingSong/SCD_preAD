#!/bin/bash

# rsem 

rsem=/home1/songlt/Software/rsem/usr/local/bin
ref_gff=/home1/songlt/ref_genome/genome/hg19/gencode.v33lift37.annotation.gtf
ref_fa=/home1/songlt/ref_genome/genome/hg19/hg19.fa
star_2pass=/home1/GENE_proc/SONGLITING/mRNA/star_2pass


mkdir -p /home1/GENE_proc/SONGLITING/mRNA/rsem
rsem_out=/home1/GENE_proc/SONGLITING/mRNA/rsem

# 1. prepare-reference
##$rsme/rsem-prepare-reference --gtf ${ref_gff} ${ref_fa} ~/ref_genome/index/RSEM/hg19



# 2. calculate-expression
# For Illumina TruSeq Stranded protocols, please use  'reverse'. (Default: 'none')
#cat /home1/GENE_proc/SONGLITING/mRNA/samples_info/samples_info.txt |sed -n '193,200p' | cut -f1 | while read sample

#do

#$rsem/rsem-calculate-expression -p 20 --strandedness reverse --alignments --paired-end  ${star_2pass}/${sample}Aligned.toTranscriptome.out.bam  /home1/songlt/ref_genome/index/RSEM/hg19 ${rsem_out}/${sample}

#done

# 3. combine 
gene_res=`ls ${rsem_out}/*genes.results`
iso_res=`ls ${rsem_out}/*isoforms.results`


$rsem/rsem-generate-data-matrix ${gene_res} > ${rsem_out}/gene_output_matrix.txt
$rsem/rsem-generate-data-matrix ${iso_res} > ${rsem_out}/iso_output_matrix.txt
