## Perform functional analysis of the differentially expressed genes obtained by DEseq2
## Figure 2
## 2021.05.10
## Song Liting

options(stringsAsFactors = F)
library(deseq2)
library(clusterProfiler)
library(stringr)
library(gProfileR)
library(gprofiler2)
library(ggpubr)
library(Hmisc)
library(ggplot2)
library(cowplot)
library(reshape2)
library(UpSetR)
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

deseq2_scd$biotype[deseq2_scd$biotype%in%c('IG_gene','protein_coding','TR_gene')] <- 'protein_coding'
deseq2_mci$biotype[deseq2_mci$biotype%in%c('IG_gene','protein_coding','TR_gene')] <- 'protein_coding'
deseq2_ad$biotype[deseq2_ad$biotype%in%c('IG_gene','protein_coding','TR_gene')] <- 'protein_coding'

deseq2_scd <- subset(deseq2_scd, biotype%in%c('protein_coding'))
deseq2_mci <- subset(deseq2_mci, biotype%in%c('protein_coding'))
deseq2_ad <- subset(deseq2_ad, biotype%in%c('protein_coding'))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 2C: Volcano plot
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

filter_data <- deseq2_scd %>%na.omit() %>%dplyr::filter(!is.infinite(log2FoldChange)) %>%
  mutate(Gene_type = case_when(
    log2FoldChange >= log2(1.3) & pvalue < 0.05 ~ pal_jco("default")(10)[4],
    log2FoldChange <= -log2(1.3) & pvalue < 0.05 ~ pal_jco("default")(10)[1],
    TRUE ~ "grey60"))

for_label <- subset(filter_data, abs(log2FoldChange) >= 1 & pvalue < 0.05)

f1_volcano <- ggplot(data = filter_data, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(shape = 21,aes(fill = Gene_type, color=Gene_type),na.rm = T,position = "jitter",size=0.001) +
  scale_color_manual(values =c(pal_jco("default")(4)[1], pal_jco("default")(10)[4],pal_jco("default")(4)[3])) +
  xlab(expression(paste(log[2],'(FoldChange)'))) +
  ylab(expression(paste(-log[10],'(',italic('P'),'-value)')))+
  geom_hline(yintercept = -log10(0.05),linetype = "dashed",size = 0.2,color = "grey50") +
  geom_vline(xintercept = c(-log2(1.3), log2(1.3)),linetype = "dashed",size = 0.2,color = "grey50") +
  theme_bw() +
  theme(panel.grid = element_blank(),legend.position = '', text = element_text(size = 7))+
  ggrepel::geom_text_repel(
    aes(label = tx_symbol),fontface  = 'italic',
    color =for_label$Gene_type,data = for_label,show.legend = F,segment.color = 'grey20', segment.size = 0.1,size = 1.5)


all_des <- rbind(deseq2_scd,deseq2_mci,deseq2_ad)
deseq2_scd <- subset(deseq2_scd, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05   )%>%arrange(pvalue)
deseq2_mci <- subset(deseq2_mci, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )%>%arrange(pvalue)
deseq2_ad <- subset(deseq2_ad, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )%>%arrange(pvalue)
deseq2_DEs_sig <- rbind(deseq2_ad,deseq2_mci,deseq2_scd)
xtabs(~group+level,deseq2_DEs_sig)

#tables2 <- deseq2_DEs_sig[,c(1:6,8:13)]
#tables2[,1:6] <- signif(tables2[,1:6],2)
#write.table(tables2, file='/Users/songlt/Desktop/AD2/figure/tables/TableS2.txt',row.names = F,sep='\t',quote = F)
deseq2_DEs_sig$group <- factor(deseq2_DEs_sig$group, levels = c('SCD','MCI','AD'))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 2A, Figure S2A,  Fold change histograms for up or down-regulated  protein-coding genes and isoforms in SCD MCI and AD
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

deseq2_DEs_sig_scd <- deseq2_DEs_sig[deseq2_DEs_sig$group=='SCD',]
deseq2_DEs_sig_scd$level[deseq2_DEs_sig_scd$level=='iso'] <- 'isoform'

f1_hist_fc_scd <- ggplot(deseq2_DEs_sig_scd[deseq2_DEs_sig_scd$log2FoldChange>0,], aes( color=level, fill=level)) +
  geom_histogram(aes(y=abs(log2FoldChange), x= ..count..), position="identity",bins = 20, alpha=0.5)+
  geom_histogram(data=deseq2_DEs_sig_scd[deseq2_DEs_sig_scd$log2FoldChange<0,], aes(y=abs(log2FoldChange),x = -..count..), position="identity",bins = 20, alpha=0.5)+
  ylim(log2(1.2),2)+theme_bw() + xlab('Number of DE genes/isoforms')+ 
  ylab(expression(paste('abs[',log[2],'(FoldChange)',']')))+
  scale_fill_manual(values=pal_jco("default")(4)[c(1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = c(0.7,0.2),legend.margin =  unit(0.01, 'cm'),legend.key.size = unit(0.2, 'cm'), legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.2, 'cm'),legend.title = element_blank(),
        text = element_text(size = 7))+coord_flip()

f1_hist_fc_mci_ad <- ggplot(deseq2_DEs_sig[deseq2_DEs_sig$log2FoldChange>0 &deseq2_DEs_sig$group!='SCD' ,], aes( color=level, fill=level)) +
  geom_histogram(aes(x=abs(log2FoldChange), y= ..count..), position="identity",bins = 20, alpha=0.5)+
  facet_grid(~group)+
  geom_histogram(data=deseq2_DEs_sig[deseq2_DEs_sig$log2FoldChange<0&deseq2_DEs_sig$group!='SCD',], aes(x=abs(log2FoldChange),y = -..count..), position="identity",bins = 20, alpha=0.5)+
  xlim(log2(1.2),3)+theme_bw() + ylab('Number of DE genes') +
  scale_fill_manual(values=pal_jco("default")(4)[c(1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = c(0.8,0.1),legend.margin =  unit(0.01, 'cm'),legend.key.size = unit(0.2, 'cm'), legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.2, 'cm'),legend.title = element_blank(), 
        text = element_text(size = 7))


f1_hist_fc_sma <- ggplot(deseq2_DEs_sig[deseq2_DEs_sig$log2FoldChange>0,], aes( color=level, fill=level)) +
  geom_histogram(aes(y=abs(log2FoldChange), x= ..count..), position="identity",bins = 20, alpha=0.5)+
  facet_grid(~group)+
  geom_histogram(data=deseq2_DEs_sig[deseq2_DEs_sig$log2FoldChange<0,], aes(y=abs(log2FoldChange),x = -..count..), position="identity",bins = 20, alpha=0.5)+
  ylim(log2(1.2),2)+theme_bw() + xlab('Number of DE genes/isoforms')+ 
  scale_fill_manual(values=pal_jco("default")(4)[c(1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = c(0.9,0.2),legend.margin =  unit(0.01, 'cm'),legend.key.size = unit(0.2, 'cm'), legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.2, 'cm'),legend.title = element_blank(), 
        text = element_text(size = 7))+coord_flip()



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 2B, Figure S2B-C,  PCA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

deseq2_gene_iso <- rbind(gene_exp$norm_exp,iso_exp$norm_exp)

f1_ad_pca <-  fviz_pca_ind(prcomp(t(deseq2_gene_iso[intersect(rownames(deseq2_gene_iso),deseq2_ad$ensg),rownames(coldata)[coldata$Diagnosis%in%c('NC','AD')] ])) ,  axes = c(1, 2), #habillage只能是分类变量
                           addEllipses = T, ellipse.type="norm",ellipse.level=0.5,label="var",pointsize=0.5,
                           habillage = coldata[rownames(coldata)[coldata$Diagnosis%in%c('NC','AD')],'Diagnosis'],
                           mean.point=T,title="NC vs. AD")+scale_color_manual(values=pal_jco("default")(4)[c(3,4)])+scale_fill_manual(values=pal_jco("default")(4)[c(3,4)])+
  theme(legend.position = c(0.1,0.1),legend.key.size = unit(0.3, 'cm'),legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'), legend.title  = element_blank(), 
        text = element_text(size = 7))


f1_mci_pca <- fviz_pca_ind(prcomp(t(deseq2_gene_iso[intersect(rownames(deseq2_gene_iso), deseq2_mci$ensg),rownames(coldata)[coldata$Diagnosis%in%c('NC','MCI')] ])) ,  axes = c(1, 2), #habillage只能是分类变量
                           addEllipses = T, ellipse.type="norm",ellipse.level=0.5,label="var",pointsize=0.5,
                           habillage = coldata[rownames(coldata)[coldata$Diagnosis%in%c('NC','MCI')],'Diagnosis'],
                           mean.point=T,title="NC vs. MCI")+scale_color_manual(values=pal_jco("default")(4)[c(3,2)])+scale_fill_manual(values=pal_jco("default")(4)[c(3,2)])+
  theme(legend.position = c(0.1,0.1),legend.key.size = unit(0.3, 'cm'),legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'), legend.title  = element_blank(), 
        text = element_text(size = 7))

f1_scd_pca <- fviz_pca_ind(prcomp(t(deseq2_gene_iso[intersect(rownames(deseq2_gene_iso), deseq2_scd$ensg),rownames(coldata)[coldata$Diagnosis%in%c('NC','SCD')] ])) , #habillage只能是分类变量
                           axes = c(1,2),addEllipses = T, ellipse.type="norm",ellipse.level=0.5,pointsize=0.5,
                           habillage = coldata[rownames(coldata)[coldata$Diagnosis%in%c('NC','SCD')],'Diagnosis'],
                           mean.point=T,title="NC vs. SCD", label="none")+scale_color_manual(values=pal_jco("default")(4)[c(3,1)])+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1)])+
  theme(legend.position = c(0.8,0.2),legend.key.size = unit(0.3, 'cm'),legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'), legend.title  = element_blank(), 
        text = element_text(size = 7))


scd_pca <- prcomp(t(deseq2_gene_iso[intersect(rownames(deseq2_gene_iso), deseq2_scd$ensg),rownames(coldata)[coldata$Diagnosis%in%c('NC','SCD')] ]))
scd_pca <- as.data.frame(scd_pca$x)
scd_pca$Diagnosis <- coldata[  rownames(scd_pca),'Diagnosis']
f1_scd_pca <- ggplot(data=as.data.frame(scd_pca),aes(x=PC1,y=PC2,color=Diagnosis,shape=Diagnosis))+
  geom_point(size=0.3,aes(color=Diagnosis,shape=Diagnosis))+
  theme_bw()+ stat_ellipse(aes(fill = Diagnosis), level = 0.3,alpha=0.2,geom="polygon")+
  scale_color_manual(values=pal_jco("default")(4)[c(3,1)])+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1)])+
  theme(legend.position = c(0.7,0.2),legend.margin =  unit(0.01, 'cm'),legend.key.size = unit(0.2, 'cm'), legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.2, 'cm'),legend.title = element_blank(), 
        text = element_text(size = 7))


plot_grid(f1_scd_pca,f1_mci_pca,f1_ad_pca)
plot_grid(f1_scd_gene_pca,f1_scd_iso_pca)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 2D  pathways    
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

scd_pat_u <- gprofiler(deseq2_scd$gene_symbol[deseq2_scd$log2FoldChange>0],src_filter=c('GO:BP','KEGG','REAC'),significant = F,min_set_size=15,max_set_size=500,hier_filtering = 'moderate')#custom_bg = ensg_sym_autosome[expd_genes,'gene_symbol']
mci_pat_u <- gprofiler(deseq2_mci$gene_symbol[deseq2_mci$log2FoldChange>0] ,src_filter=c('GO:BP','KEGG','REAC'),significant = F,min_set_size=15,max_set_size=500,hier_filtering = 'moderate')
ad_pat_u <- gprofiler(deseq2_ad$gene_symbol[deseq2_ad$log2FoldChange>0] ,src_filter=c('GO:BP','KEGG','REAC'),significant = F,min_set_size=15,max_set_size=500,hier_filtering = 'moderate')

scd_pat_d <- gprofiler(deseq2_scd$gene_symbol[deseq2_scd$log2FoldChange<0],src_filter=c('GO:BP','KEGG','REAC'),significant = F,min_set_size=15,max_set_size=500,hier_filtering = 'moderate')
mci_pat_d <- gprofiler(deseq2_mci$gene_symbol[deseq2_mci$log2FoldChange<0] ,src_filter=c('GO:BP','KEGG','REAC'),significant = F,min_set_size=15,max_set_size=500,hier_filtering = 'moderate')
ad_pat_d <- gprofiler(unique(deseq2_ad$gene_symbol[deseq2_ad$log2FoldChange<0]) ,src_filter=c('GO:BP','KEGG','REAC'),significant = F,min_set_size=15,max_set_size=500,hier_filtering = 'moderate')

scd_pat_u$group <- 'SCD'; scd_pat_u$direction <- 'upregulated'
scd_pat_d$group <- 'SCD'; scd_pat_d$direction <- 'downregulated'
mci_pat_u$group <- 'MCI'; mci_pat_u$direction <- 'upregulated'
mci_pat_d$group <- 'MCI'; mci_pat_d$direction <- 'downregulated'
ad_pat_u$group <- 'AD'; ad_pat_u$direction <- 'upregulated'
ad_pat_d$group <- 'AD'; ad_pat_d$direction <- 'downregulated'

scd_pathws_top <- rbind(scd_pat_u%>%arrange(p.value)%>% distinct(.,term.name,.keep_all= TRUE) %>% head(.,5),scd_pat_d%>%arrange(p.value)%>% distinct(.,term.name,.keep_all= TRUE) %>% head(.,5))

f1_scd_pat <- ggplot(scd_pathws_top ,aes(x=reorder(term.name, p.value),y=-log10(p.value), color=direction,fill=direction))+
  geom_bar(stat='identity',width=0.5)+facet_wrap(.~direction,scales="free",nrow = 1) +
  #facet_grid(vars(level),vars(direction),scales="free")+
  theme_bw()+ 
  theme(legend.position = '',axis.text.x = element_text(lineheight = 0.8,angle=90,hjust = 1,vjust =0.5 ,size=6),
        plot.margin=unit(c(0.5,0.5,0.5,0.5),'lines'), 
        axis.text.y = element_text(size=7),text = element_text(size = 7)) + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))+
  geom_hline(yintercept = -log10(0.05),color=2,linetype = "dashed" ) +
  xlab('') + scale_color_manual(values=pal_jco("default")(3))+
  ylab(expression(paste(-log[10],'(FDR)')))+
  scale_fill_manual(values=pal_jco("default")(3))

# IPA

DE_ad_IPA <- read.table('./data/IPA/DE_AD_pat.txt',sep='\t',header = T)
DE_mci_IPA <- read.table('./data/IPA/de_mci_pat.txt',sep='\t',header = T)
DE_scd_IPA <- read.table('./data/IPA/de_scd_pat.txt',sep='\t',header = T)
colnames(DE_scd_IPA)[2] <- c('-log10(P-value)')
DS_scd_IPA <- read.table('./data/IPA/ds_scd_pat.txt',sep='\t',header = T)


fig_scd_ipa <- ggplot(DE_scd_IPA[1:10,], aes(x=reorder(`Ingenuity.Canonical.Pathways`,-`-log10(P-value)`),
                                             y=`-log10(P-value)`, fill=z.score))+
  geom_bar(stat = 'identity',width = 0.5,color='grey')+
  theme_bw() +
  theme(legend.position = c(0.7,0.8),legend.direction = "horizontal",
        legend.margin =  unit(0.01, 'cm'),legend.key.height = unit(0.3, 'cm'),legend.key.width = unit(0.3, 'cm'),
        panel.grid.major=element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1,vjust =0.5 ,size=7),plot.margin=unit(c(0.5,0.5,0.5,0.2),'lines'),
        text=element_text(size=7),axis.title.y = element_text(vjust=-1)) + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=35))+scale_color_jco()+
  geom_hline(yintercept = -log10(0.05),color=2,linetype = "dashed" ) +xlab('') + 
  ylab(expression(paste(-log[10],'(',italic('P'),'-value)')))+
  scale_fill_gradient2(low = "#0073C2FF",mid = "white",high = "#EFC000FF",midpoint = 0)


fig_scd_ipa <- ggplot(DE_scd_IPA[1:10,], aes(x=reorder(`Ingenuity.Canonical.Pathways`,-Ratio),y=`Ratio` ,
                                             size=`-log10(P-value)`, fill=z.score ))+
  geom_point( color='grey',shape = 21,stroke = 0.5) +
  theme_bw() +
  theme(legend.position = 'top', axis.text.x = element_text(angle = 90, hjust = 1,vjust =0.5 ,size=6),
        legend.margin = unit(10,"mm"),legend.key.height = unit(0.25, 'cm'),legend.key.width = unit(0.3, 'cm'),
        #legend.key.height = unit(0.4, 'cm'),legend.key.width = unit(0.5, 'cm'),
        text=element_text(size=7),axis.title.y = element_text(vjust=-1)) + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))+ylim(c(0,0.3))+
  ylab('Gene Ratio')+xlab('')+labs( size=expression(paste(-log[10],'(',italic('P'),'-value)')))+
  scale_fill_gradient2(low = "#0073C2FF",mid = "white",high = "#EFC000FF",midpoint = 0)


ggplot(DE_scd_IPA[c(2,4,5,11),], aes(x=reorder(`Ingenuity.Canonical.Pathways`,-`-log10(P-value)`),
                                     y=`-log10(P-value)`, fill=z.score))+
  geom_bar(stat = 'identity',width = 0.5)+
  theme_bw() +
  theme(legend.position = c(0.8,0.8),legend.title = element_blank(),
        legend.margin =  unit(0.01, 'cm'),legend.key.height = unit(0.2, 'cm'),legend.key.width = unit(0.2, 'cm'),
        panel.grid.major=element_blank(),
        axis.text.y = element_text(size=8),plot.margin=unit(c(0.5,0.5,0.5,0.2),'lines'),
        text=element_text(size=8),axis.title.y = element_text(vjust=-1)) + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=30))+scale_color_jco()+
  geom_hline(yintercept = -log10(0.05),color=2,linetype = "dashed" ) +xlab('') + 
  ylab(expression(paste(-log[10],'(',italic('P'),'-value)')))+coord_flip()+
  scale_fill_gradient2(low = "#0073C2FF",mid = "white",high = "#EFC000FF",midpoint = 0)
dev.print(pdf, file='./figure/fig1_pat_scd.pdf')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 2G: expression of ISG genes (IFI44L, IFI27, RSAD2, SIGLEC1, IFIT1, and IS15) defining the IFN signature
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

load('./data/pathway_genes.RData')
#load('./data/DEs.RData')
coldata$Diagnosis <- factor(coldata$Diagnosis, levels=c('NC','SCD','MCI','AD'))

DE_pathway_genes <- Reduce(intersect, list(deseq2_scd$gene_symbol, pathway_genes, deseq2_ad$gene_symbol ))

ISG_genes_scd <- deseq2_scd[deseq2_scd$gene_symbol%in%DE_pathway_genes,'ensg']
ISG_genes_ad <- deseq2_ad[deseq2_ad$gene_symbol%in%DE_pathway_genes,'ensg']
ISG_genes_mci <- deseq2_mci[deseq2_mci$gene_symbol%in%DE_pathway_genes,'ensg']


fc_int <- subset(all_des, ensg%in%Reduce(intersect, 
                                         list(rownames(ensg_sym_autosome)[ensg_sym_autosome$gene_symbol%in% pathway_genes], 
                                              deseq2_scd$ensg, deseq2_ad$ensg)))
fc_int <- subset(all_des, ensg%in%union(ISG_genes_scd,ISG_genes_ad ))

fc_int$group <- factor(fc_int$group, levels = c('SCD','MCI','AD'))
fc_int$level[fc_int$level=='iso'] <- 'isoform'

f1_fc_int <- ggplot(fc_int, aes(x=group, y=log2FoldChange, group=ensg, color=level, fill=level))+
  geom_line(alpha=0.5)+
  geom_point(size=0.5)+
  theme_bw() +  theme(legend.position = c(0.1,0.1),legend.title = element_blank())+xlab('') +scale_fill_jco()+scale_color_jco()+
  ggrepel::geom_text_repel(aes(label = tx_symbol),data = fc_int[abs(fc_int$log2FoldChange)>=1 ,],
                           show.legend = F,segment.color = 'grey50',color='grey50' ,size=1.5,segment.size = 0.1,fontface = 'italic')+
  theme(legend.position = c(0.8,0.2),legend.margin =  unit(0.01, 'cm'),
        text = element_text(size = 7),
        legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank())+ 
  ylab(expression(paste(log[2],'(FoldChange)'))) +
  geom_hline(yintercept = 0, linetype = "dashed",size = 0.5,color = 2) 

ISG_gene <- c('IFI44L', 'IFI27', 'RSAD2', 'SIGLEC1', 'IFIT1', 'ISG15')

fc_int <- subset(all_des, gene_symbol%in%ISG_gene)

fc_int$group <- factor(fc_int$group, levels = c('SCD','MCI','AD'))
fc_int$level[fc_int$level=='iso'] <- 'isoform'
f2_fc_int <- ggplot(fc_int, aes(x=group, y=log2FoldChange, group=ensg, color=level, fill=level))+
  geom_line(alpha=0.5)+
  geom_point(size=1)+
  theme_bw() +  theme(legend.position = c(0.1,0.1),legend.title = element_blank())+xlab('') +scale_fill_jco()+scale_color_jco()+
  ggrepel::geom_text_repel(aes(label = tx_symbol),data = fc_int[abs(fc_int$log2FoldChange)>=log2(1.3)&fc_int$pvalue<0.05 ,],
                           show.legend = F,segment.color = 'grey50',color='grey50' ,size=1.5,segment.size = 0.1,fontface = 'italic')+
  theme(legend.position = c(0.8,0.2),legend.margin =  unit(0.01, 'cm'),
        text = element_text(size = 7),
        legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank())+ 
  ylab(expression(paste(log[2],'(FoldChange)'))) +
  geom_hline(yintercept = 0, linetype = "dashed",size = 0.5,color = 2) 


dev.print(pdf, file='./figure/fs_fc_int.pdf',width=6,height=3)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 2H: ISG score based on expression of six genes (IFI44L, IFI27, RSAD2, SIGLEC1, IFIT1, and ISG15) 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

ISG_gene <- c('IFI44L', 'IFI27', 'RSAD2', 'SIGLEC1', 'IFIT1', 'ISG15')
#ISG_gene <- read.table('./data/IPA/interferon.genes.txt',sep='\t')[,1]
ISG_trans <- rownames(ensg_sym_autosome)[ensg_sym_autosome$gene_symbol%in%ISG_gene]
ISG_exp <- gene_exp$norm_exp[intersect(ISG_trans,rownames(gene_exp$norm_exp)),]
rownames(ISG_exp) <- ensg_sym_autosome[rownames(ISG_exp),'gene_symbol']
m_ISG_exp <- melt(ISG_exp)
m_ISG_exp$dia <- coldata[m_ISG_exp$Var2,'Diagnosis']
m_ISG_exp$dia <- factor(m_ISG_exp$dia, levels = c('NC','SCD','MCI','AD'))


ISG_score <- as.data.frame(rowMeans(apply(ISG_exp,1,scale)))
ISG_score <- as.data.frame(colMeans(ISG_exp))

colnames(ISG_score) <- 'ISG_score'
ISG_score$dia <- coldata[,'Diagnosis']
ISG_score$dia <- factor(ISG_score$dia, levels = c('NC','SCD','MCI','AD'))
save(ISG_score, file='~/Desktop/AD2/data/ISG_score.RData')


library(Rmisc)
ISG_score_s <- summarySE(ISG_score, measurevar="ISG_score", groupvars=c("dia"))

f1_ISG_score_box <- ggplot(ISG_score,aes(dia,ISG_score,group=dia,fill=dia,color=dia))+
  geom_boxplot(alpha=0.2,width=0.5)+
  #geom_violin(alpha=0,width=0.8)+
  stat_summary(fun = mean, geom = "point", size=1)+
  stat_summary(fun = mean, geom = "line",aes(group = 1),color=2)+
  theme_bw() + ylab('ISG score')+xlab('')+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = '', text = element_text(size = 7))+
  geom_signif(test='t.test', map_signif_level=F,step_increase = 0.1, margin_top=0.5,size = 0.4, textsize = 2,
              comparisons = list(c("NC", "SCD"),c("NC", "MCI"),c("NC", "AD"), c("SCD", "MCI"),c("SCD", "AD") )) 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# GSEA
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

load('./data/DEs_0412.RData')
gene_deseq2_scd <- deseq2_scd[deseq2_scd$level=='gene',]
gene_deseq2_scd <- dplyr::distinct(gene_deseq2_scd,gene_symbol,.keep_all=TRUE)

gene_deseq2_mci <- deseq2_mci[deseq2_mci$level=='gene',]
gene_deseq2_mci <- dplyr::distinct(gene_deseq2_mci,gene_symbol,.keep_all=TRUE)

gene_deseq2_ad <- deseq2_ad[deseq2_ad$level=='gene',]
gene_deseq2_ad <- dplyr::distinct(gene_deseq2_ad,gene_symbol,.keep_all=TRUE)

geneList_scd<-gene_deseq2_scd$log2FoldChange 
names(geneList_scd)=gene_deseq2_scd$gene_symbol 
geneList_scd=sort(geneList_scd,decreasing = T) 

geneList_mci<-gene_deseq2_mci$log2FoldChange 
names(geneList_mci)=gene_deseq2_mci$gene_symbol
geneList_mci=sort(geneList_mci,decreasing = T) 

geneList_ad<-gene_deseq2_ad$log2FoldChange 
names(geneList_ad)=gene_deseq2_ad$gene_symbol 
geneList_ad=sort(geneList_ad,decreasing = T) 

#kegmt<-read.gmt("/Users/songlt/Desktop/AD2/data/gsea/h.all.v7.2.symbols.gmt") #读gmt文件
kegmt<-read.gmt("/Users/songlt/Desktop/AD2/data/gsea/c5.go.bp.v7.2.symbols.gmt") 

interferon_genes <- subset(kegmt,term==c('GO_TYPE_I_INTERFERON_SIGNALING_PATHWAY'))[,2]
save(interferon_genes,file='/Users/songlt/Desktop/AD2/data/typeI_interferon_genes.RData')
write.table(interferon_genes, file='./figure/tables/typeI_IFN_gene.txt',row.names = F, quote = F,col.names = F)

KEGG_scd<-GSEA(geneList_scd,TERM2GENE = kegmt,minGSSize=15,maxGSSize=500) #GSEA分析
#f1_gsea_scd <- gseaplot2(KEGG_scd,1,pvalue_table = F, base_size = 8,color = 'red',rel_heights = c(1.5, 0.3, 1)) 

KEGG_mci<-GSEA(geneList_mci,TERM2GENE = kegmt,minGSSize=15,maxGSSize=500) #GSEA分析
#f1_gsea_mci <- gseaplot2(KEGG_mci,1,pvalue_table = F, base_size = 5,color = 'red',rel_heights = c(1.5, 0.2, 1)) 

KEGG_ad<-GSEA(geneList_ad,TERM2GENE = kegmt,minGSSize=15,maxGSSize=500) #GSEA分析
#f1_gsea_ad <- gseaplot2(KEGG_ad,10,pvalue_table = F,base_size =5,color = 'red',rel_heights = c(1.5, 0.2, 1)) 

aaa <- KEGG_ad@result
bbb <- KEGG_scd@result
ccc <- KEGG_mci@result

aaa$group <- 'AD'
bbb$group <- 'SCD'
ccc$group <- 'MCI'

gsea_res <- rbind(aaa,bbb,ccc)

gsea_res <- subset(gsea_res, ID%in%c('GO_TYPE_I_INTERFERON_SIGNALING_PATHWAY'))

gsea_res$group <- factor(gsea_res$group, levels = c('SCD','MCI','AD'))
gsea_res$ID <- gsub('type i','type I', tolower(gsub('_',' ',gsub(pattern = '^GO_','',gsea_res$ID))))
gsea_res$ID <- gsub(' pathway','', gsea_res$ID)
gsea_res$ID <- factor(gsea_res$ID, levels = c('type I interferon signaling','innate immune response','adaptive immune response'))
gsea_res$star <- '***'
gsea_res$star[gsea_res$p.adjust>0.01] <- '**'
gsea_res$star[gsea_res$p.adjust>0.05] <- '*'
gsea_res$star[gsea_res$p.adjust>0.1] <- '+'
gsea_res$star_pos  <- gsea_res$NES+0.3
gsea_res$star_pos[gsea_res$NES<0] <- gsea_res$NES[gsea_res$NES<0] - 0.4

f1_gsea_plot <- ggplot( gsea_res, aes(x=group,y=NES, color=group, fill=group))+
  geom_bar(stat = 'identity',width = 0.5)+facet_wrap(.~ID,nrow=1) +
  geom_text(data=gsea_res,aes(x=group,y=star_pos),size = 2,label=gsea_res$star,color='black',)+
  theme_bw()+ 
  ylab('Normalized enrichment score')+ xlab('')+
  scale_y_continuous( 
    limits = c(-3,3))+  theme(legend.position = '',
                              text = element_text(size = 7),strip.text = element_text(size = 5),axis.title.y = element_text(size=6)) + 
  scale_fill_manual(values=pal_jco("default")(4)[c(1,2,4)])+
  scale_color_manual(values=pal_jco("default")(4)[c(1,2,4)])


# plot_grid

plot_grid(plot_grid(f1_hist_fc_scd,f1_scd_pca, f1_volcano,align ='h',nrow=1,labels=c('A','B','C'),label_size = 9),
          plot_grid( f1_scd_pat,NA,align = "v",axis = "b",rel_widths = c(1,1), labels=c('D','E'),label_size = 9, nrow=1),  
          plot_grid( f1_gsea_plot,f2_fc_int, f1_ISG_score_vio, nrow=1,rel_widths = c(1,1,1),labels=c('F','G','H'),label_size = 9), 
          nrow=3, rel_heights  = c(1,1.4,1))

dev.print(pdf,file='./figure/fig2_color_bg4.pdf', width = 7, height = 7.5)


