## Differential expression analysis using Deseq2
## 2021.04.12
## Song Liting

library(edgeR)
library(DESeq2)
library(stringr)
library(VennDiagram)

# s1 expressed genes and transcripts
setwd('/Users/songlt/Desktop/AD2/')
load('./data/genes_samples.RData')

# functions
# f1 不同标化方法得到的结果
coldata$Age_group <- factor(coldata$Age_group, levels = c('youth','middle_aged','yonger_aged','aged'))
#coldata <- coldata[!rownames(coldata)%in%outliers1,]
exp_data <- function(level){
  
  # tpm
  tpm <- ifelse (level=='gene', './data/gene_tpm_matrix.txt', './data/isoform_tpm_matrix.txt')
  tpm <- as.matrix(read.table(tpm,header = T, stringsAsFactors = F, row.names = 1))
  colnames(tpm) <- gsub('X.home1.GENE_proc.SONGLITING.mRNA.rsem.|.genes.results|.isoforms.results','',colnames(tpm))
  tpm <- tpm[rownames(tpm)%in%c(expd_genes,expd_trans),rownames(coldata)]
  #tpm <- tpm[rownames(cts),colnames(cts)]
  
  # count
  cts <- ifelse (level=='gene', './data/gene_output_matrix.txt', './data/iso_output_matrix.txt')
  cts <- as.matrix(read.table(cts,header = T, stringsAsFactors = F, row.names = 1))
  colnames(cts) <- gsub('X.home1.GENE_proc.SONGLITING.mRNA.rsem.|.genes.results|.isoforms.results','',colnames(cts))
  cts <- cts[rownames(cts)%in%c(expd_genes,expd_trans),rownames(coldata)]
  cts <- round(cts)
  
  # log2(count+1)
  log2_cts <- log2(cts+1)
  
  # vst
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~Gender + Age_group + RIN + Diagnosis)
  
  dds <- estimateSizeFactors( dds )
  
  #dds <- dds[rowSums(counts(dds))>=20,]
  vsd <- assay(vst(dds))
  
  vsd_1 <- (vst(dds))
  
  norm_exp <- assay( vst(dds)) # log2(gene_exp$normalized_counts[1:3,1:3]+1)； normalized_counts <- counts(dds,normalized=T) 
  
  return(list(count=cts, lg_count=log2_cts, tpm= tpm, vsd = vsd, vsd_1=vsd_1,norm_exp=norm_exp))
  
}

# DEseq2
DE_seq2 <- function(ref,cond,level,coldata){
  # ref : NC
  # cond: SCD,MCI and AD
  # level: gene or transcrpt
  
  coldata <- subset(coldata, Diagnosis%in%c(ref,cond))
  coldata$Diagnosis <- factor(coldata$Diagnosis, levels=c(ref,cond))

  # count matrix
  cts <- gene_exp$count
  if(level!='gene'){cts <- iso_exp$count }
  cts <- round(cts[,rownames(coldata)])
  
  # design= ~ batch + condition
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~Gender + Age_group + RIN + Diagnosis)
  dds <- estimateSizeFactors( dds )
  
  ## DEseq
  dds <- DESeq(dds)
  
  norm_exp <- assay( normTransform(dds)) # log2(gene_exp$normalized_counts[1:3,1:3]+1)； normalized_counts <- counts(dds,normalized=T) 
  norm_exp <- assay( vst(dds)) # log2(gene_exp$normalized_counts[1:3,1:3]+1)； normalized_counts <- counts(dds,normalized=T) 
  
  
  ## DEseq
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  res$FC_med <- apply(norm_exp, 1, function(x) {median(x[which(coldata$Diagnosis==cond)])-median(x[which(coldata$Diagnosis==ref)] )})
  #res$FC_mean <- apply(norm_exp, 1, function(x) {mean(x[which(coldata$Diagnosis==cond)])/mean(x[which(coldata$Diagnosis==ref)]) )})
  
  #res$group <- paste(ref, cond, sep = '_')
  res$group <- cond
  res$level <- level
  res$ensg <- rownames(res)
  res$gene_symbol <- ensg_sym_autosome[rownames(res), 'gene_symbol']
  res$tx_symbol <- ensg_sym_autosome[rownames(res), 'tx_symbol']
  
  res$biotype <- ensg_sym_autosome[rownames(res), 'biotype']
  res$tx_biotype <- ensg_sym_autosome[rownames(res), 'tx_biotype']
  
  #res <- subset(res, biotype%in%c('lncRNA','protein_coding')  )
  #res$lFDR <- lfdr(res$pvalue)
  #res$BH <- p.adjust(res$pvalue,'BH')
  
  res_p <- res
  #res_p <- subset(res, abs(FC_med) > log2(1.3) & pvalue<0.05)
  return(res_p)
  
}

############################################################
#################### main programes ####################
############################################################


# s1 Normalization of expression
gene_exp <- exp_data('gene')
iso_exp <- exp_data('iso')

#  deseq2: Differential expression analysis
scd_gene_deseq2 <- DE_seq2('NC','SCD','gene',coldata=coldata)
mci_gene_deseq2 <- DE_seq2('NC','MCI','gene',coldata=coldata)
ad_gene_deseq2 <- DE_seq2('NC','AD','gene',coldata=coldata)

scd_iso_deseq2 <- DE_seq2('NC','SCD','iso',coldata=coldata)
mci_iso_deseq2 <- DE_seq2('NC','MCI','iso',coldata=coldata)
ad_iso_deseq2 <- DE_seq2('NC','AD','iso',coldata=coldata)

deseq2_scd <- rbind(scd_gene_deseq2, scd_iso_deseq2)
deseq2_mci <- rbind(mci_gene_deseq2, mci_iso_deseq2)
deseq2_ad <- rbind(ad_gene_deseq2, ad_iso_deseq2)

save(scd_gene_deseq2, mci_gene_deseq2,ad_gene_deseq2, scd_iso_deseq2, mci_iso_deseq2, ad_iso_deseq2, deseq2_scd,deseq2_mci,  deseq2_ad, 
     file='./data/deseq2_mrna_rmotl_0412.RData')
save.image(file='./data/DEs_0412.RData')

deseq2_scd <- subset(deseq2_scd, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )
deseq2_mci <- subset(deseq2_mci, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )
deseq2_ad <- subset(deseq2_ad, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )

deseq2_DEs_sig <- rbind(deseq2_ad,deseq2_mci,deseq2_scd)
xtabs(~group+level,deseq2_DEs_sig)


write.csv(deseq2_scd, file='./figure/deseq2_mrna_scd.csv')
write.csv(deseq2_mci, file='./figure/deseq2_mrna_mci.csv')
write.csv(deseq2_ad, file='./figure/deseq2_mrna_ad.csv')


