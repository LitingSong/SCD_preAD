#!/bin/bash

source /home1/songlt/miniconda3/bin/activate
mature_mir=/home1/songlt/ref_genome/genome/miRNA/matureU2T_ns.fa
hairpin_mir=/home1/songlt/ref_genome/genome/miRNA/hairpinU2T_ns.fa
#Gorilla gorilla 
ggo_mir=/home1/songlt/ref_genome/genome/miRNA/ggo_mir_ns.fa
#Pan troglodytes (
ptr_mir=/home1/songlt/ref_genome/genome/miRNA/ptr_mir_ns.fa
#Pongo pygmaeus 
ppy_mir=/home1/songlt/ref_genome/genome/miRNA/ppy_mir_ns.fa
#other
other_mir=/home1/songlt/ref_genome/genome/miRNA/other.fa
softw=/home1/songlt/miniconda3/bin
ref_fa=/home1/songlt/ref_genome/genome/hg19/hg19.fa

dedup_fq=/home1/GENE_proc/SONGLITING/miRNA/dedup_fq
raw_fq=/home2/Cohorts/Genes/NeuroDegenerative/RNA/Project_s905s09001_244Samples_20200430

#cat /home1/GENE_proc/SONGLITING/miRNA/samples_info/samples_info.txt | sed -n '1,244p' |cut -f1 | while read sample
#
#do
#
#cd /home1/GENE_proc/SONGLITING/miRNA/dedup_fq
#
## fastq to  fasta
#nohup $softw/mapper.pl  ${dedup_fq}/${sample}.R1.fastq  -e -h  -j  -l 18 -m -o 8 -p /home1/songlt/ref_genome/index/bowtie/hg19 -s ${sample}_reads_collapsed.fa -t ${sample}_reads_vs_refdb.arf &
#
#done

cd /home1/GENE_proc/SONGLITING/miRNA/mirdeep2
#$softw/mapper.pl config.txt -d -c -i -j -l 18 -m -p /home1/songlt/ref_genome/index/bowtie/hg19  -s reads.fa -t reads_vs_genome.arf

# Identification of known and novel miRNAs
#$softw/miRDeep2.pl ${sample}_reads_collapsed.fa ${ref_fa} ${sample}_reads_vs_refdb.arf ${mature_mir} ${other_mir}  ${hairpin_mir} -t hsa 2 > ${sample}.log  
$softw/miRDeep2.pl reads.fa ${ref_fa}  reads_vs_genome.arf  ${mature_mir} ${other_mir}  ${hairpin_mir} -t hsa -d 2 > ${sample}.log  

