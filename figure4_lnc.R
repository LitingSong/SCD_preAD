# lncRNA analysis
# figure 4
# Song Liting
# 2021.04.10


options(stringsAsFactors = F)
library(deseq2)
library(clusterProfiler)
library(stringr)
library(gProfileR)
library(gprofiler2)
library(ggplot2)
library(cowplot)
library(reshape2)
library(UpSetR)
library(Rmisc)
library(ggsci)
library(dplyr)
library(VennDiagram)
library(FactoMineR)
library(FactoMineR)
library(factoextra)
library(mclust)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(enrichplot)

setwd('/Users/songlt/Desktop/AD2/')

load('./data/DEs_0412.RData') # output of mRNA_deseq2.R

coldata$Diagnosis <- factor(coldata$Diagnosis, levels=c('NC','SCD','MCI','AD'))

lnc_deseq2_scd <- subset(deseq2_scd, biotype%in%c('lncRNA'))
lnc_deseq2_mci <- subset(deseq2_mci, biotype%in%c('lncRNA'))
lnc_deseq2_ad <- subset(deseq2_ad, biotype%in%c('lncRNA'))

all_des <- rbind(lnc_deseq2_scd,lnc_deseq2_mci,lnc_deseq2_ad)
lnc_deseq2_scd <- subset(lnc_deseq2_scd, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05   )
lnc_deseq2_mci <- subset(lnc_deseq2_mci, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )
lnc_deseq2_ad <- subset(lnc_deseq2_ad, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )
lnc_deseq2_DEs_sig <- rbind(lnc_deseq2_ad,lnc_deseq2_mci,lnc_deseq2_scd)
lnc_deseq2_DEs_sig$gene_symbol <- ensg_sym_autosome[lnc_deseq2_DEs_sig$ensg,'gene_symbol']
lnc_deseq2_DEs_sig$tx_symbol <- ensg_sym_autosome[lnc_deseq2_DEs_sig$ensg,'tx_symbol']

lnc_deseq2_DEs_sig$group <- factor(lnc_deseq2_DEs_sig$group, levels = c('SCD','MCI','AD'))

Table_s6 <- lnc_deseq2_DEs_sig
Table_s6 <- lnc_deseq2_DEs_sig[,c(1:6,8:13)]
Table_s6[,1:6] <- signif(Table_s6[,1:6],2)
write.table(Table_s6, file='/Users/songlt/Desktop/AD2/figure/tables/TableS6.txt',row.names = F,sep='\t',quote = F)

xtabs(~group+level,lnc_deseq2_DEs_sig)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 4A. The number of differentially expressed long non-coding genes 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
lnc_sum <- melt(xtabs(~group+level,lnc_deseq2_DEs_sig))
lnc_sum$level <- as.character(lnc_sum$level)
lnc_sum$level[lnc_sum$level=='iso'] <- 'isoform'
f3_de_lunc_box <- ggplot(lnc_sum,aes(x=group,y=value,fill=level,group=level))+
  geom_text(aes(label=value),position=position_dodge(width=0.9), vjust=0,size=2.5)+
  geom_bar(stat = 'identity',position = 'dodge')+
  scale_fill_jco()+
  scale_color_jco()+
  theme_bw()+xlab('')+ylab('Number of DE lncRNAs')+
  theme(legend.position = c(0.2,0.8),legend.key.size = unit(0.3, 'cm'),legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'), legend.margin =  unit(0.01, 'cm'),
        legend.text = element_text(size=6),legend.title  = element_blank(),text = element_text(size = 7))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Figure 4B Volcano plot of lncRNA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
load('./data/DEs_0412.RData') # 是 差异表达的结果
lnc_deseq2_scd <- subset(deseq2_scd, biotype%in%c('lncRNA'))
lnc_deseq2_mci <- subset(deseq2_mci, biotype%in%c('lncRNA'))
lnc_deseq2_ad <- subset(deseq2_ad, biotype%in%c('lncRNA'))

library(dplyr)
filter_data <- lnc_deseq2_scd %>%na.omit() %>%dplyr::filter(!is.infinite(log2FoldChange)) %>%
  mutate(Gene_type = case_when(
    log2FoldChange >= log2(1.5) & pvalue < 0.05 ~ pal_jco("default")(4)[4],
    log2FoldChange <= -log2(1.5) & pvalue < 0.05 ~ pal_jco("default")(4)[1],
    TRUE ~ 'grey'))

for_label <- subset(filter_data, abs(log2FoldChange) >= log2(1.5) & pvalue < 0.05)

f3_volcano_scd <- ggplot(data = filter_data, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(shape = 21,fill = filter_data$Gene_type, color=filter_data$Gene_type,na.rm = T,position = "jitter",size=0.1) +
  scale_fill_manual(values =c('grey', pal_jco("default")(4)[4],pal_jco("default")(4)[1] )) +
  scale_color_manual(values =c('grey', pal_jco("default")(4)[4],pal_jco("default")(4)[1])) +
  xlab(expression(paste(log[2],'(FoldChange)'))) +
  ylab(expression(paste(-log[10],'(',italic('P'),'-value)')))+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed",size = 0.5,color = "grey50") +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)),linetype = "dashed",size = 0.5,color = "grey50") +
  theme_bw() +
  theme(panel.grid = element_blank(),legend.position = '', text = element_text(size = 7))+
  ggrepel::geom_text_repel(
    aes(label = tx_symbol),fontface  = 'italic',
    color =for_label$Gene_type,data = for_label,show.legend = F,segment.color = 'grey20', segment.size=0.3,size = 2)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## Figure 4C. Expression of NRIR 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

library(ggpubr)
ggplot_exp_df <- function(gene_norm,genes){
  
  df <- t(gene_norm[intersect(rownames(gene_norm),genes),])
  #colnames(df) <- ensg_sym_autosome[colnames(df), 'gene_symbol']
  mdf <- melt(df)
  if (nrow(df)==1){  colnames(mdf) <- c('gene','sample','exp'); mdf$gene <- genes;  mdf$diagnosis <- coldata[mdf$sample,'Diagnosis']; mdf$diagnosis <- factor(mdf$diagnosis, levels=c('NC','SCD','MCI','AD'))}
  if (nrow(df)>=1){   colnames(mdf) <- c('sample','gene','exp');mdf$diagnosis <- coldata[mdf$sample,'Diagnosis'];mdf$diagnosis <- factor(mdf$diagnosis, levels=c('NC','SCD','MCI','AD')) }
  
  return(mdf)
}
pg <- c('NRIR')
pgt <- rownames(ensg_sym_autosome)[ensg_sym_autosome$gene_symbol%in%pg]
pg_exp <- subset(ggplot_exp_df(rbind(gene_exp$vsd,iso_exp$norm_exp), genes = pgt))

#pg_exp$gene <- str_split_fixed(as.character(pg_exp$gene),pattern = '\\.',2)[,1]
pg_exp$gene <- ensg_sym_autosome[as.character(pg_exp$gene),'tx_symbol']
NRIR_exp <- subset(pg_exp,gene=='NRIR')

NRIR_exp_s <- summarySE(NRIR_exp, measurevar="exp", groupvars=c("diagnosis"))

f3_NRIR_vio <- ggplot(NRIR_exp,aes(diagnosis,exp,fill=diagnosis, color=diagnosis))+
  geom_violin(alpha=0.3,width=0.5)+facet_wrap(.~gene,nrow = 1, scales = 'free_y')+
  #geom_bar(stat = 'identity',width = 0.5)+
  geom_point(data = NRIR_exp_s,aes(x=diagnosis, y=exp),pch=19,size=0.5)+ #绘制均值为点图
  geom_errorbar(data=NRIR_exp_s, aes(ymin=exp-ci, ymax=exp+ci),width=.1)+
  theme_bw() + ylab('Normalized expression')+xlab('')+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+
  scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = '', text = element_text(size = 7),strip.text=element_text(face  = 'italic') ) +
  geom_signif(test='t.test',annotations = c("**", "ns.",'ns.'), map_signif_level=T,step_increase = 0.2, margin_top=0.5,size = 0.4, textsize = 2,
              comparisons = list(c("NC", "SCD"),c("NC", "MCI"),c("NC", "AD") )) 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 4H correlation between NRIR and ISG_score 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

ISG_gene <- c('IFI44L', 'IFI27', 'RSAD2', 'SIGLEC1', 'IFIT1', 'ISG15')
ISG_trans <- rownames(ensg_sym_autosome)[ensg_sym_autosome$gene_symbol%in%ISG_gene]
ISG_exp <- gene_exp$norm_exp[intersect(ISG_trans,rownames(gene_exp$norm_exp)),]
rownames(ISG_exp) <- ensg_sym_autosome[rownames(ISG_exp),'gene_symbol']
m_ISG_exp <- melt(ISG_exp)
m_ISG_exp$dia <- coldata[m_ISG_exp$Var2,'Diagnosis']
m_ISG_exp$dia <- factor(m_ISG_exp$dia, levels = c('NC','SCD','MCI','AD'))


ISG_score <- as.data.frame(colMeans(ISG_exp))
colnames(ISG_score) <- 'ISG_score'
ISG_score$dia <- coldata[,'Diagnosis']
ISG_score$dia <- factor(ISG_score$dia, levels = c('NC','SCD','MCI','AD'))

ISG_score$NRIR <- gene_exp$norm_exp['ENSG00000225964.6_7',]

library(ggpubr)
ggscatter(ISG_score, x = 'NRIR', y = 'ISG_score',color = "dia",
          add = "reg.line",  # Add regressin line
          add.params = list(color = "steelblue"), 
          conf.int = TRUE,size=0.1) + 
  labs(x=expression(paste("Expression of ", italic("NRIR"))), y="ISG score")+
  stat_cor(method = "pearson",size = 2) +theme_bw()+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = c(0.8,0.2),legend.margin =  unit(0.03, 'cm'),legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank(), legend.text = element_text(size=6),
        text = element_text(size = 8))

f3_ISG_NRIR_scater <- ggscatter(ISG_score, x = 'NRIR', y = 'ISG_score',color = "dia",
                                add = "reg.line",  # Add regressin line
                                add.params = list(color = "steelblue"), 
                                conf.int = TRUE,size=0.1) + 
  labs(x=expression(paste("Expression of ", italic("NRIR"))), y="ISG score")+
  annotate("text",label=expression(paste(italic('R'),' = ',0.85,', ',italic('P'),' < ',2.2e-16)),x=7.5,y=13.7, color=1,size=2) +theme_bw()+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = c(0.8,0.2),legend.margin =  unit(0.03, 'cm'),legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.25, 'cm'),
        legend.key.width = unit(0.25, 'cm'),legend.title = element_blank(), legend.text = element_text(size=6),
        text = element_text(size = 7))



# fig 4 lncrna
f4_lncrna <- plot_grid(f3_de_lunc_box,f3_volcano_scd,f3_NRIR_vio,f3_ISG_NRIR_scater, nrow = 1,labels=c('A','B','C','D'),label_size = 9,align = "h",axis = "b")

# fig 4 mirna
f4_mirna <- plot_grid(f4_de_mir_bar, f4_mir146_isg_scater, f4_mir146_box,f4_pc1_pc2_scater,labels = c('E','F','G','H'),nrow = 1,label_size = 9,axis = "b",align = "h")

plot_grid(f4_lncrna,f4,nrow=2) 

dev.print(pdf,file='./figure/fig4_nc.pdf', width = 7, height = 4)

