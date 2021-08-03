### mirna analysis
### 2020.10.20
### Song Liting


library(DESeq2)
library(ggsignif)
library(corrplot)
library(WGCNA)
library(ggplot2)
library(pheatmap)
library(UpSetR)
library(multiMiR)
library(clusterProfiler)
library(stringr)
library(reshape2)
library(dplyr)
library(FactoMineR)
library(FactoMineR)
library("factoextra")
library(mclust)
library(Rmisc)
library(ggsci)
library(ggpubr)

setwd('/Users/songlt/Desktop/AD2/')
load('./data/genes_samples.RData')

config <- read.table('./data/config.txt', sep='\t', comment.char = '' ,header = F, row.names = 2)
config$sample <- str_split_fixed(config$V1, '_',2)[,1]
coldata$Age_group <- factor(coldata$Age_group, levels = c('youth','middle_aged','yonger_aged','aged'))
coldata$Diagnosis <- factor(coldata$Diagnosis, levels = c('NC','SCD','MCI','AD'))

mirna <- read.table('./data/miRNAs_expressed_all_samples_22_06_2020_t_12_10_32.csv', sep='\t', comment.char = '', header = T, check.names = F )
rownames(mirna) <- paste(mirna$precursor, mirna$miRNA, sep='_')
mirna <- mirna[,c(5:248) ]
colnames(mirna) <- config[colnames(mirna),'sample']
mirna <- mirna[!duplicated(str_split_fixed(rownames(mirna),'_',3)[,3]),]
rownames(mirna) <- str_split_fixed(rownames(mirna),'_',3)[,3]

mirna <- mirna[, rownames(coldata)]

# step2: miRNAs: reads count > 3 in 50% sample
mirna <- mirna[rowSums((mirna) >  3) >=  nrow(coldata)*0.5, ]

mirna_norm <- DESeqDataSetFromMatrix(countData = mirna,
                                     colData = coldata,
                                     design = ~Gender + Age_group + Diagnosis)

mirna_vst <- assay( normTransform(mirna_norm)) # log2(gene_exp$normalized_counts[1:3,1:3]+1)； normalized_counts <- counts(dds,normalized=T) 


# deseq2
mirna_DE <- function(ref, cond, cts,coldata){


  #coldata <- read.table('/home1/GENE_proc/SONGLITING/mRNA/DEseq2/coldata.txt',header = T, stringsAsFactors = F, row.names = 1)
  coldata <- subset(coldata, Diagnosis%in%c(ref,cond))
  coldata$Diagnosis <- factor(coldata$Diagnosis, levels=c(ref,cond))
  
  # design= ~ batch + condition
  dds <- DESeqDataSetFromMatrix(countData = cts[,rownames(coldata)],
                                colData = coldata,
                                design = ~Gender + Age_group + Diagnosis)
  
  ## DEseq
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  
  res$group <- cond
  
  res$BH <- p.adjust(res$pvalue,'BH')
  
  #res <- subset(res,  abs(log2FoldChange) > log2(1.2) & pvalue< 0.5)
  res$mirna <- gsub('hsa-', '', str_split_fixed(rownames(res), '_',2)[,1])
  res$ma <-  rownames(res)
  return(res)
  
}

nc_ad_mir <- mirna_DE('NC','AD', mirna, coldata)
nc_scd_mir <- mirna_DE('NC','SCD', mirna, coldata)
nc_mci_mir <- mirna_DE('NC','MCI', mirna, coldata)


DEmir_all <- rbind(nc_ad_mir, nc_scd_mir,  nc_mci_mir)

DEmir <- subset(DEmir_all, abs(log2FoldChange) > log2(1.3) & pvalue<0.05)
#DEmir <- subset(DEmir_all,  pvalue<0.05)

DEmir$dir <- 'up'
DEmir$dir[DEmir$log2FoldChange<0] <- 'down'

DEmir$mir_p <- DEmir$ma
xtabs(~group+dir,DEmir)
# tables7 <- DEmir[,c(1:7,9:10)]
# tables7[,1:6] <- signif(tables7[,1:6],2)
# write.table(tables7, file='/Users/songlt/Desktop/AD2/figure/tables/TableS7.txt',row.names = F,sep='\t',quote = F)
save(mirna_vst, nc_ad_mir, nc_scd_mir, nc_mci_mir,DEmir, file='./data/deseq2_mirna.RData' )


mir_sum <- melt(xtabs(~group+dir,DEmir))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 4E The number of differentially expressed miRNAs
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

f4_de_mir_bar <- ggplot(mir_sum,aes(x=factor(group,levels = c('SCD','MCI','AD')),y=value,fill=dir,group=dir))+
  geom_text(aes(label=value),position=position_dodge(width=0.9), vjust=0,size=2.5)+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_jco()+
  scale_color_jco()+
  theme_bw()+xlab('')+ylab('Number of DE miRNAs')+
  theme(legend.position = c(0.2,0.8),legend.key.size = unit(0.3, 'cm'),legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'), legend.margin =  unit(0.01, 'cm'),
        legend.text = element_text(size=6),legend.title  = element_blank(),text = element_text(size = 7))


#dev.print(pdf, file='./figure/DE_mir_venn.pdf')

#inf_mir <- read.table('./data/mir_ipa_interferon.txt',sep='\t',header = T)# ipa的基因
inf_mir <- read.table("./data/BP_typeI_mir.txt",sep='\t',header = T) # go bp 基因

inf_mir$ID <- str_split_fixed(inf_mir$ID,':|\\*',2)[,1]

inf1_scd_mir <- merge(inf_mir, nc_scd_mir,by.x='ID',by.y='ma')

#sig_inf1_scd_mir <- subset(inf1_scd_mir,  `p.value` < 0.05) # ipa
sig_inf1_scd_mir <- subset(inf1_scd_mir,  `q.value.FDR.B.H` < 0.05) # GO:BP DUI

sig_inf1_scd_mir$log2FoldChange <- signif(sig_inf1_scd_mir$log2FoldChange,2)
sig_inf1_scd_mir$p.value <- signif(sig_inf1_scd_mir$p.value,2)
table_sig_inf1_scd_mir <- sig_inf1_scd_mir[order(sig_inf1_scd_mir$p.value),c(1,14,4,9,10,5,6)]
colnames(table_sig_inf1_scd_mir) <- c("ID", "log2FC", "Source","Hit Count \n in Query List", "Hit Count \n in Genome",   
                                      "P.value", "q.value \n Bonferroni")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 4H correlation between PC1 of mirnas and PC1 of IFN genes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

interferon_genes <- read.table("./figure/tables/typeI_IFN_gene.txt")[,1] # 以GO:BP的结果为主
#interferon_genes <- read.table("/Users/songlt/Desktop/AD2/data/IPA/interferon.genes.txt",sep='\t')[,'V1'] # 以IPA的结果为主

mir <- mirna_vst[unique(sig_inf1_scd_mir$ID),rownames(coldata)[coldata$Diagnosis%in%c('NC','SCD','MCI','AD')]]

rna <- gene_exp$norm_exp[intersect(rownames(gene_exp$norm_exp), rownames(ensg_sym_autosome)[ensg_sym_autosome$gene_symbol%in%interferon_genes]), 
                         rownames(coldata)[coldata$Diagnosis%in%c('NC','SCD','MCI','AD')]]

pca_MIR <- prcomp((mir), rank. = 2)
pca_RNA <- prcomp((rna), rank. = 2)

mir_rna <- as.data.frame(cbind(pca_MIR$rotation,pca_RNA$rotation))
mir_rna$dia <- coldata[rownames(mir_rna),'Diagnosis']
colnames(mir_rna) <- c("PC1", "PC2", "PC1.1", "PC2.1" ,"dia")
mir_rna$dia <- factor(mir_rna$dia,levels = c('NC','SCD','MCI','AD'))


ggscatter(mir_rna, x = 'PC1', y = 'PC1.1',color = "dia",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "steelblue"), 
          conf.int = TRUE,size=0.1) + 
  labs(x='PC1 of miRNAs',y="PC1 of interferon mRNAs")+theme_bw()+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  stat_cor(method = "pearson",size=2)+
  theme(legend.position = c(0.8,0.8),legend.margin =  unit(0.03, 'cm'),legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank(), legend.text = element_text(size=6),
        text = element_text(size = 7))

f4_pc1_pc2_scater <- ggscatter(mir_rna, x = 'PC1', y = 'PC1.1',color = "dia",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "steelblue"), 
          conf.int = TRUE,size=0.1) + 
  labs(x='PC1 of miRNAs',y="PC1 of type I IFN signaling genes")+theme_bw()+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+
  scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = '',legend.margin =  unit(0.03, 'cm'),legend.key.size = unit(0.3, 'cm'), 
        legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank(), legend.text = element_text(size=6),
        text = element_text(size = 7))+
  annotate("text",label=expression(paste(italic('R'),' = ',-0.37,', ',italic('P'),' = 5.8e-08')), color=1,size=2,x=-0.075,y=-0.08) 
  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 4G expression of hsa-miR-146a-5p
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

de_mirs <- unique(sig_inf1_scd_mir$ID)

ggplot_exp_df <- function(de_mirs){
  df <- t(mirna_vst)[,c(de_mirs)]
  df <- melt(df)
  colnames(df) <- c('sample','gene','exp')
  df$sample <- as.character(df$sample)
  
  df$diagnosis <- coldata[df$sample,'Diagnosis']
  df$diagnosis <- factor(df$diagnosis, levels=c('NC','SCD','MCI','AD'))
  return(df)
}

mir146_isg <- subset(ggplot_exp_df(de_mirs),gene=='hsa-miR-146a-5p')
mir146_isg$ISG_score <- ISG_score$ISG_score
mir146_isg$STAT1 <- gene_exp$vsd[rownames(subset(ensg_sym_autosome,gene_symbol=='STAT1'&tx_symbol=='STAT1')),]


mir146_exp_s <- summarySE(mir146_isg, measurevar="exp", groupvars=c("diagnosis"))

f4_mir146_box <- ggplot(mir146_isg,aes(diagnosis,exp,fill=diagnosis, color=diagnosis))+
  geom_violin(alpha=0.3,width=0.5)+facet_wrap(.~gene,nrow = 1, scales = 'free_y')+
  #geom_bar(stat = 'identity',width = 0.5)+
  geom_point(data = mir146_exp_s,aes(x=diagnosis, y=exp),pch=19,size=0.5)+ #绘制均值为点图
  geom_errorbar(data=mir146_exp_s, aes(ymin=exp-ci, ymax=exp+ci),width=.1)+
  theme_bw() + ylab('Normalized expression')+xlab('') +theme(legend.position = '', text = element_text(size = 7),strip.text=element_text(face  = 'italic') )+ ylab('Normalized expression')+xlab('')+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)]) # + #scale_fill_jco()+scale_color_jco()+
# geom_signif(test='t.test', map_signif_level=F,step_increase = 0.1, margin_top=0.5,size = 0.4, textsize = 2,
#             comparisons = list(c("NC", "SCD"),c("NC", "MCI"),c("NC", "AD") )) #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 4H correlation between hsa-miR-145a-5p and STAT1 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

ggscatter(mir146_isg, x = 'exp', y = 'STAT1',color = "diagnosis",
                                  add = "reg.line",  # Add regressin line
                                  add.params = list(color = "steelblue"), 
                                  conf.int = TRUE,size=0.1) + 
  labs(x=expression(paste("Expression of ", italic("hsa-miR-145a-5p"))),
       y=expression(paste("Expression of ", italic("STAT1")))) +
  stat_cor(method = "pearson",size = 2) +
  theme_bw()+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = c(0.8,0.2),legend.margin =  unit(0.03, 'cm'),legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank(), legend.text = element_text(size=6),
        text = element_text(size = 8))

f4_mir146_isg_scater <- ggscatter(mir146_isg, x = 'exp', y = 'STAT1',color = "diagnosis",
                                  add = "reg.line",  # Add regressin line
                                  add.params = list(color = "steelblue"), 
                                  conf.int = TRUE,size=0.1) + 
  labs(x=expression(paste("Expression of ", italic("hsa-miR-145a-5p"))),
       y=expression(paste("Expression of ", italic("STAT1")))) +
  theme_bw()+
  annotate("text",label=expression(paste(italic('R'),' = ',-0.26,', ',italic('P'),' = ',0.021)), color=1,size=2,x=13,y=14.6) +theme_bw()+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = '',legend.margin =  unit(0.03, 'cm'),legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank(), legend.text = element_text(size=6),
        text = element_text(size = 7))




