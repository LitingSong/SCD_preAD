#!/bin/bash


#Step 1. Converting bams to juncs
mkdir -p /home1/GENE_proc/SONGLITING/mRNA/leafcutter

cd /home1/GENE_proc/SONGLITING/mRNA/leafcutter


#for bamfile in `ls /home1/GENE_proc/SONGLITING/mRNA/star_2pass/*Aligned.out.bam|sed -n '2,80p'`
#do
#    echo Converting $bamfile to $bamfile.junc
#    nohup sh  /home1/songlt/Software/leafcutter/scripts/bam2junc.sh $bamfile $bamfile.junc &
#done

#mv /home1/GENE_proc/SONGLITING/mRNA/star_2pass/*Aligned.out.bam.junc ./
#ls *.junc > juncfiles.txt
# add diagnosis to juncfiles.txt -> group file
# The command line interface currently only supports two groups. 

##Step 2. Intron clustering
source /home1/songlt/miniconda3/bin/activate
python2 /home1/songlt/Software/leafcutter/clustering/leafcutter_cluster.py -j ./juncfile_nc_mci.txt -m 50 -o NC_MCI  -l 500000
#
###Step 3. Differential intron excision analysis
/home1/songlt/Software/leafcutter/scripts/leafcutter_ds.R --num_threads 20 --exon_file=/home1/songlt/Software/leafcutter/leafcutter/data/gencode19_exons.txt.gz NC_MCI_perind_numers.counts.gz group_nc_mci.txt
#
mv leafcutter_ds_cluster_significance.txt NC_MCI_leafcutter_ds_cluster_significance.txt
mv leafcutter_ds_effect_sizes.txt NC_MCI_leafcutter_ds_effect_sizes.txt
###Step 4. Plotting splice junctions
##/home1/songlt/Software/leafcutter/scripts/ds_plots.R -e /home1/songlt/Software/leafcutter/leafcutter/data/gencode19_exons.txt.gz NC_MCI_perind_numers.counts.gz group_nc_mci.txt  leafcutter_ds_cluster_significance.txt -f 0.05
#
###Generate the annotation files
##/home1/songlt/Software/leafcutter/leafviz/gtf2leafcutter.pl -o /home1/songlt/ref_genome/genome/hg19 /home1/songlt/ref_genome/genome/hg19/gencode.v33lift37.annotation.gtf
##
#Prepare the LeafCutter differential splicing results for visualisation
/home1/songlt/Software/leafcutter/leafviz/prepare_results2.R -f 0.05 -m group_nc_mci.txt NC_MCI_perind_numers.counts.gz NC_MCI_leafcutter_ds_cluster_significance.txt NC_MCI_leafcutter_ds_effect_sizes.txt /home1/songlt/ref_genome/genome/hg19/
mv leafviz.RData NC_MCI_leafviz_pvalue.RData
##
###Visualise the result
##
#~/Software/leafcutter/leafviz/run_leafviz.R ./leafviz.RData
#
##source deactivate
#
#
