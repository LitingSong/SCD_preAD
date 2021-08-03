## local splicing 
## Figure 3
## 2021.05.10
## Song Liting

library(grid)
library(DESeq2)
library(clusterProfiler)
library(stringr)
library(gProfileR)
library(gprofiler2)
library(ggplot2)
library(cowplot)
library(reshape2)
library(UpSetR)
library(VennDiagram)
library(ggsci)
library(ensembldb)
library(ggsignif)
library(PFAM.db)
library(Gviz)
library(stringr)
library(dplyr)
library(ggsci)
library(EnsDb.Hsapiens.v75)


setwd('/Users/songlt/Desktop/AD2/')

load('./data/DEs_0412.RData')

all_des <- rbind(deseq2_scd,deseq2_mci,deseq2_ad )
deseq2_scd <- subset(deseq2_scd, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )
deseq2_mci <- subset(deseq2_mci, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )
deseq2_ad <- subset(deseq2_ad, abs(log2FoldChange) > log2(1.3) &pvalue < 0.05 )

load('~/Desktop/AD2/data/leafcutter_pvalue/202/NC_SCD_leafviz_pvalue.RData') # output of 4_leafcutter_nc_scd_pvalue.sh
scd_count <- counts
nc_scd_c <- clusters
nc_scd_c$padj <- p.adjust(nc_scd_c$FDR, method = 'fdr')
nc_scd_c <- nc_scd_c[nc_scd_c$FDR< 0.001,]
nc_scd_c$group <- 'nc_scd'
nc_scd_c$gene <- gsub(pattern = '<i>|</i>','',nc_scd_c$gene)
nc_scd_c <- nc_scd_c[nc_scd_c$gene!='.',]
nc_scd_c <- subset(nc_scd_c,gene%in%ensg_sym_autosome[expd_genes,'gene_symbol'])

load('~/Desktop/AD2/data/leafcutter_pvalue/202/NC_AD_leafviz_pvalue.RData') # output of 4_leafcutter_nc_scd_pvalue.sh
ad_count <- counts
nc_ad_c <- clusters
nc_ad_c$padj <- p.adjust(nc_ad_c$FDR,method = 'fdr')
nc_ad_c <- nc_ad_c[nc_ad_c$FDR< 0.001,]
nc_ad_c$group <- 'nc_ad'
nc_ad_c$gene <- gsub(pattern = '<i>|</i>','',nc_ad_c$gene)
nc_ad_c <- nc_ad_c[nc_ad_c$gene!='.',]
nc_ad_c <- subset(nc_ad_c,gene%in%ensg_sym_autosome[expd_genes,'gene_symbol'])

load('~/Desktop/AD2/data/leafcutter_pvalue/202/NC_MCI_leafviz_pvalue.RData') # output of 4_leafcutter_nc_mci_pvalue.sh
mci_count <- counts
nc_mci_c <- clusters
nc_mci_c$padj <- p.adjust(nc_mci_c$FDR, method = 'fdr')
nc_mci_c <- nc_mci_c[nc_mci_c$FDR< 0.001,]
nc_mci_c$group <- 'nc_mci'
nc_mci_c$gene <- gsub(pattern = '<i>|</i>','',nc_mci_c$gene)
nc_mci_c <- nc_mci_c[nc_mci_c$gene!='.',]
nc_mci_c <- subset(nc_mci_c,gene%in%ensg_sym_autosome[expd_genes,'gene_symbol'])

combine_splicing <- rbind(nc_scd_c,nc_mci_c, nc_ad_c)

#Table_s4 <-  combine_splicing[,c(1:6,8)]
#Table_s4$group <- toupper( gsub('nc_','',Table_s4$group))
#write.table(Table_s4, file='./figure/tables/TableS4.txt',sep='\t',quote = F,row.names = F)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 3A overlap between DS and DEG/DET
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

f2_DS_venn <- venn.diagram(list(SCD=(nc_scd_c$gene),MCI =(nc_mci_c$gene),AD = (nc_ad_c$gene)),
                   fill = pal_jco("default")(4)[c(1,2,4)] , alpha = 0.7, filename = NULL,col='white',
                   cex=0.5,cat.cex = 0.5);


f2_scd_mrna_splice_venn <- plot_grid( venn.diagram(list(DS=nc_scd_c$gene,
                                                        DEG = deseq2_scd$gene[deseq2_scd$level=='gene'],
                                                        DET = deseq2_scd$gene[deseq2_scd$level=='iso']),
                                fill =  pal_jco("default")(4)[c(1,2,4)], alpha = 0.7, filename = NULL,col='white',
                                cex=0.5,cat.cex = 0.5,
                                cat.just=list(c(0.5,1),c(0.5,1),c(0.5,-15))))

load('/Users/songlt/Desktop/AD2/data/typeI_interferon_genes.RData')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 3B overlap between DS and DEG/DET of IFN genes
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

interferon_genes_venn <-  plot_grid( venn.diagram(list( DS=intersect(nc_scd_c$gene,interferon_genes),
                                                        DEG = intersect(deseq2_scd$gene[deseq2_scd$level=='gene'],interferon_genes),
                                                        DET = intersect(deseq2_scd$gene[deseq2_scd$level=='iso'],interferon_genes)),
                             fill =  pal_jco("default")(4)[c(1,2,4)], alpha = 0.7, filename = NULL,col='white',cex=0.5,cat.cex = 0.5))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 3C 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
interferon_DS_bar <- ggplot(subset(nc_scd_c, gene%in%interferon_genes),aes(x=reorder(gene,FDR),y=-log10(FDR)))+
  geom_bar(stat = 'identity',width = 0.5,fill = pal_jco("default")(4)[c(1)])+
  theme_bw()+xlab('')+  ylab(expression(paste(-log[10],'(',italic('P'),'-value)')))+
  theme(axis.text.x = element_text( face="italic"),text=element_text(size = 7))


AD_genes <- read.table('~/Desktop/AD2/data/AD.combined2.txt')$V1


# Hypergeometric test, differential gene transcripts, differential splicing
phyper_test <- function(group){
  bg <- unique(ensg_sym_autosome[expd_genes,'gene_symbol'])
  
  if (group == 'SCD'){DE <- deseq2_scd$gene[deseq2_scd$level=='gene'];
  DT <- deseq2_scd$gene[deseq2_scd$level=='iso'];
  DS <- unique(nc_scd_c$gene) }
  
  if (group == 'MCI'){DE <- deseq2_mci$gene[deseq2_mci$level=='gene'];
  DT <- deseq2_mci$gene[deseq2_mci$level=='iso'];
  DS <- unique(nc_mci_c$gene) }
  
  if (group == 'AD'){DE <- deseq2_ad$gene[deseq2_ad$level=='gene'];
  DT <- deseq2_ad$gene[deseq2_ad$level=='iso'];
  DS <- unique(nc_ad_c$gene) }
  
  p1 <- phyper(length(intersect(DT,DE))-1,length(DT),length(setdiff(bg,DT)),length(DE),lower.tail = F)
  p2 <- phyper(length(intersect(DE,DS))-1,length(DE),length(setdiff(bg,DE)),length(DS),lower.tail = F)
  p3 <- phyper(length(intersect(DT,DS))-1,length(DT),length(setdiff(bg,DT)),length(DS),lower.tail = F)
  
  return(c(p1,p2,p3))
}

phyper_test('SCD')
phyper_test('MCI')
phyper_test('AD')

phyper_test1 <- function(group1, group2){
  
  bg <- unique(ensg_sym_autosome[expd_genes,'gene_symbol'])
  
  DT <- unique(subset(combine_splicing,group==group1 )[,'gene'])
  DS <- unique(subset(combine_splicing,group==group2 )[,'gene'])
  
  p3 <- phyper(length(intersect(DT,DS))-1,length(DT),length(setdiff(bg,DT)),length(DS),lower.tail = F)
  
  return(c(p3))
}

phyper_test1('nc_scd','nc_mci')
phyper_test1('nc_scd','nc_ad')
phyper_test1('nc_mci','nc_ad')

phyper_test_interferon <- function(){
  
  bg <- interferon_genes
  
  DE <- intersect(deseq2_scd$gene[deseq2_scd$level=='gene'],interferon_genes)
  DT <- intersect(deseq2_scd$gene[deseq2_scd$level=='iso'],interferon_genes)
  DS <- intersect(unique(nc_scd_c$gene) ,interferon_genes)

  p1 <- phyper(length(intersect(DT,DE))-1,length(DT),length(setdiff(bg,DT)),length(DE),lower.tail = F)
  p2 <- phyper(length(intersect(DE,DS))-1,length(DE),length(setdiff(bg,DE)),length(DS),lower.tail = F)
  p3 <- phyper(length(intersect(DT,DS))-1,length(DT),length(setdiff(bg,DT)),length(DS),lower.tail = F)
  
  return(c(p1,p2,p3))
  
}
phyper_test_interferon()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 3D  STAT1  structure
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

{
  pfam <- read.table('/Users/songlt/Documents/database/pdb_pfam_mapping.txt',sep='\t',header = T,quote = '"')[,c(5,6)]
  pfam <- unique(pfam)
  rownames(pfam) <- pfam$PFAM_ACC
  
  #de_pgt <- unique(subset(rbind(deseq2_ad,deseq2_scd,deseq2_mci),gene_symbol=='STAT1')$ensg)
  de_pgt <- unique(subset(rbind(deseq2_scd),gene_symbol=='STAT1')$ensg)
  
  edb <- EnsDb.Hsapiens.v75
  
  DB <- ensDbFromGtf('/Users/songlt/Documents/database/STAT1.gtf', organism="Homo_sapiens", genomeVersion="GRCh37",version=1)
  EDB <- EnsDb(DB)
  
  #txs1 <- getGeneRegionTrackForGviz(EDB, filter = ~ genename == c("STAT1"))
  txs1 <- getGeneRegionTrackForGviz(EDB, filter = ~ tx_id%in% de_pgt)
  txs1@strand@values <- factor('-', levels= c('+', '-' , '*'))
  txs1$symbol <- ensg_sym_autosome[txs1$transcript,'tx_symbol']
  txs1$transcript <- str_split_fixed(txs1$transcript,'\\.',2)[,1]
  
  txs <-  txs1
  
  i <- which(unique(str_split_fixed( txs$transcript,pattern = '\\.',2)[,1])%in%c(  "ENST00000361099"))[1]
  pdoms <- proteins(edb, filter = ~ tx_id %in% unique(str_split_fixed( txs$transcript,pattern = '\\.',2)[,1])[i] &
                      protein_domain_source == "pfam",
                    columns = c("protein_domain_id", "prot_dom_start",
                                "prot_dom_end"))
  
  
  pdoms_rng <- IRanges(start = pdoms$prot_dom_start, end = pdoms$prot_dom_end,
                       names = pdoms$protein_id)
  
  pdoms_gnm <- proteinToGenome(pdoms_rng, edb)
  
  
  
  ## Convert the list to a GRanges with grouping information
  pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))
  pdoms_gnm_grng$id <- rep(pdoms$protein_domain_id, lengths(pdoms_gnm))
  pdoms_gnm_grng$pfam_name <- pfam[pdoms_gnm_grng$id,'PFAM_Name']
  pdoms_gnm_grng$grp <- rep(1:nrow(pdoms), lengths(pdoms_gnm))
  
  library(Gviz)
  
  ## Define the individual tracks:
  ## - Ideagram
  #ideo_track <- IdeogramTrack(genome = "hg19", chromosome = "chr2")
  ## - Genome axis
  gaxis_track <- GenomeAxisTrack()
  ## - Transcripts
  gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",stacking="full",
                                fill=pal_jco("default")(4)[1],col=NULL,
                                name = "Transcripts", geneSymbol = TRUE, size = 0.3)
  ## - Protein domains
  pdom_track <- AnnotationTrack(pdoms_gnm_grng, group = pdoms_gnm_grng$pfam_name,
                                id = pdoms_gnm_grng$pfam_name, groupAnnotation = "pfam_name",
                                just.group = "right", shape = "box",fill=pal_jco("default")(4)[2],col=NULL,
                                name = "Domains", size = 0.3,stacking="full")
  
  #gnm_pos <- GRanges("chr2", IRanges(191874730, width = 4014))
  gnm_pos <- GRanges(c("chr2","chr2"), IRanges(c(191874730, 191878379),width = c(4014,20)))# SCD
  
  
  hl_track <- HighlightTrack(list(gene_track, pdom_track), range = gnm_pos,inBackground=F, lty=0,fill=pal_jco("default")(4)[c(1,4)], alpha=c(0.2,1))
  
  plotTracks(list( gaxis_track, hl_track),transcriptAnnotation = "symbol",cex=1,
             cex.group=1,fontsize.group=6,fontsize=8,lwd.title=0.1,lwd.title=0.1,
             title.width=0.5,extend.right=1000,extend.left=600)#transcriptAnnotation = "transcript"
  
  }
library(lattice)
myPanel <- function(x, ...) plotTracks(list( gaxis_track, hl_track),transcriptAnnotation = "symbol",cex=1,
                                       cex.group=1,fontsize.group=6,fontsize=8,lwd.title=0.1,lwd.title=0.1,extend.right=1000,extend.left=600,
                                       title.width=0.5,panel.only=TRUE)

a <- seq(191833060-2000, 191882529+9000, len=2)
f2_ensgb_stat <- xyplot(b~a, data.frame(a=a, b=1), panel=myPanel, xlab='',draw=F,alternating=0,
                        par.settings = list(axis.line = list(col = 0)),#  不要四周的边框
                        ylab='',scales=list(x=list(col = 0), y=list(at=NULL), draw = FALSE),alpha=0,height=0.5,border=F) 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure S4D NRF1  structure
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

{
  pfam <- read.table('/Users/songlt/Documents/database/pdb_pfam_mapping2.txt',sep='\t',header = T,quote = '"')[,c(5,6)]
  pfam <- unique(pfam)
  rownames(pfam) <- pfam$PFAM_ACC
  
  #de_pgt <- unique(subset(rbind(deseq2_ad,deseq2_scd,deseq2_mci),gene_symbol=='STAT1')$ensg)
  edb <- EnsDb.Hsapiens.v75
  
  DB <- ensDbFromGtf('/Users/songlt/Documents/database/NRF1.gtf', organism="Homo_sapiens", genomeVersion="GRCh37",version=1)
  EDB <- EnsDb(DB)
  
  txs1 <- getGeneRegionTrackForGviz(EDB, filter = ~ genename == c("NRF1"))
  txs1@strand@values <- factor('-', levels= c('+', '-' , '*'))
  txs1$symbol <- ensg_sym_autosome[txs1$transcript,'tx_symbol']
  txs1$transcript <- str_split_fixed(txs1$transcript,'\\.',2)[,1]
  
  txs <-  txs1
  
  i <- which(unique(str_split_fixed( txs$transcript,pattern = '\\.',2)[,1])%in%c(  "ENST00000223190"))[1]
  pdoms <- proteins(edb, filter = ~ tx_id %in% unique(str_split_fixed( txs$transcript,pattern = '\\.',2)[,1])[i] &
                      protein_domain_source == "pfam",
                    columns = c("protein_domain_id", "prot_dom_start",
                                "prot_dom_end"))
  
  
  pdoms_rng <- IRanges(start = pdoms$prot_dom_start, end = pdoms$prot_dom_end,
                       names = pdoms$protein_id)
  
  pdoms_gnm <- proteinToGenome(pdoms_rng, edb)
  
  
  
  ## Convert the list to a GRanges with grouping information
  pdoms_gnm_grng <- unlist(GRangesList(pdoms_gnm))
  pdoms_gnm_grng$id <- rep(pdoms$protein_domain_id, lengths(pdoms_gnm))
  pdoms_gnm_grng$pfam_name <- pfam[pdoms_gnm_grng$id,'PFAM_Name']
  pdoms_gnm_grng$grp <- rep(1:nrow(pdoms), lengths(pdoms_gnm))
  
  library(Gviz)
  
  ## Define the individual tracks:
  ## - Ideagram
  #ideo_track <- IdeogramTrack(genome = "hg19", chromosome = "chr2")
  ## - Genome axis
  gaxis_track <- GenomeAxisTrack()
  ## - Transcripts
  gene_track <- GeneRegionTrack(txs, showId = TRUE, just.group = "right",stacking="full",
                                fill=pal_jco("default")(4)[1],col=NULL,
                                name = "Transcripts", geneSymbol = TRUE, size = 0.3)
  ## - Protein domains
  pdom_track <- AnnotationTrack(pdoms_gnm_grng, group = pdoms_gnm_grng$pfam_name,
                                id = pdoms_gnm_grng$pfam_name, groupAnnotation = "pfam_name",
                                just.group = "right", shape = "box",fill=pal_jco("default")(4)[2],col=NULL,
                                name = "Domains", size = 0.3,stacking="full")
  
  gnm_pos <- GRanges(c("chr7","chr7"), IRanges(c(129357216, 129297414),width = c(9871,20059)))# SCD

  
  hl_track <- HighlightTrack(list(gene_track, pdom_track), range = gnm_pos,inBackground=F, lty=0,fill=pal_jco("default")(4)[1], alpha=0.2)
  
  plotTracks(list( gaxis_track, hl_track),transcriptAnnotation = "symbol",cex=1,
             cex.group=1,fontsize.group=6,fontsize=8,lwd.title=0.1,lwd.title=0.1,
             title.width=0.5,extend.right=15000,extend.left=400, 
             from=129251561, to=129396922,sizes=c(0.15,1,0.2))#transcriptAnnotation = "transcript"
  dev.print(pdf, file='./figure/nrf1_tx_domain.pdf',width=7,height=3)
}

######################## sample level psi:#######################################
#scd

scd_ratios = scd_count %>% 
  mutate(clu = str_split_fixed(rownames(scd_count), ":", 4)[,4]) %>%
  group_by(clu) %>% 
  mutate_all( funs( ./sum(.) ) ) %>% 
  ungroup() %>%
  as.data.frame()

scd_ratios <- scd_ratios[,-ncol(scd_ratios)]

rownames(scd_ratios) <- rownames(scd_count)

scd_ratios = scd_ratios[rowMeans(is.na(scd_ratios)) <= 0.4,,drop=F ]
row_means = rowMeans(scd_ratios, na.rm = T)
row_means_outer = outer(row_means, rep(1,ncol(scd_ratios)))
scd_ratios[is.na(scd_ratios)] = row_means_outer[is.na(scd_ratios)]
#mci
mci_ratios = mci_count %>% 
  mutate(clu = str_split_fixed(rownames(mci_count), ":", 4)[,4]) %>%
  group_by(clu) %>% 
  mutate_all( funs( ./sum(.) ) ) %>% 
  ungroup() %>%
  as.data.frame()

mci_ratios <- mci_ratios[,-ncol(mci_ratios)]

rownames(mci_ratios) <- rownames(mci_count)

mci_ratios = mci_ratios[rowMeans(is.na(mci_ratios)) <= 0.4,,drop=F ]
row_means = rowMeans(mci_ratios, na.rm = T)
row_means_outer = outer(row_means, rep(1,ncol(mci_ratios)))
mci_ratios[is.na(mci_ratios)] = row_means_outer[is.na(mci_ratios)]

#ad
ad_ratios = ad_count %>% 
  mutate(clu = str_split_fixed(rownames(ad_count), ":", 4)[,4]) %>%
  group_by(clu) %>% 
  mutate_all( funs( ./sum(.) ) ) %>% 
  ungroup() %>%
  as.data.frame()

ad_ratios <- ad_ratios[,-ncol(ad_ratios)]

rownames(ad_ratios) <- rownames(ad_count)

ad_ratios = ad_ratios[rowMeans(is.na(ad_ratios)) <= 0.4,,drop=F ]
row_means = rowMeans(ad_ratios, na.rm = T)
row_means_outer = outer(row_means, rep(1,ncol(ad_ratios)))
ad_ratios[is.na(ad_ratios)] = row_means_outer[is.na(ad_ratios)]

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 3E. STAT1 cluster chr2:191874730:191878744  PSI  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

STAT1_clu_scd <- melt(scd_ratios['chr2:191874730:191878744:clu_3888_NA',])
STAT1_clu_mci <- melt(mci_ratios['chr2:191874730:191878744:clu_3909_NA',])
STAT1_clu_ad <- melt(ad_ratios['chr2:191874730:191878744:clu_3811_NA',])

STAT1_clu <- unique(rbind(STAT1_clu_scd,STAT1_clu_mci,STAT1_clu_ad))

rownames(STAT1_clu) <- gsub('Aligned.out.bam','',STAT1_clu$variable)
STAT1_clu$dia <- coldata[rownames(STAT1_clu),'Diagnosis']
STAT1_clu$dia <- factor(STAT1_clu$dia, levels=c('NC','SCD','MCI','AD'))
STAT1_clu$value <- (STAT1_clu$value)*100
#STAT1_clu$color<- 'white'
STAT1_clu$color[STAT1_clu$value>10] <- pal_jco("default")(4)[c(4)]
f2_STAT1_PSI_box <- ggplot(subset(STAT1_clu,dia!='AD'), aes(x=dia,y=value,fill=dia, color=dia)) +
  geom_boxplot(alpha=0.2,width=0.5, size=0.3) +theme_bw() +
  theme(legend.position = '',text = element_text(size = 7),strip.text=element_text(face  = 'italic') )+ 
  geom_point(color=subset(STAT1_clu,dia!='AD')$color, size=0.3)+
  ylab('sample level PSI (%)')+xlab('')+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2)])+ #scale_fill_jco()+scale_color_jco()+
  geom_signif( map_signif_level=T,annotations = c("***", "ns."),step_increase = 0.15, margin_top=0.5,size = 0.4, textsize = 2,
               test = "wilcox.test",comparisons = list(c("NC", "SCD"),c("NC", "MCI") ))#,y_position = c(12)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 3G: DET of STAT1   
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

all_des$group <- factor(all_des$group, levels = c('NC','SCD','MCI','AD'))
STAT1_FC <- subset(all_des,tx_symbol%in% subset(deseq2_scd,gene_symbol=='STAT1')$tx_symbol )
f2_stat1_fc <-  ggplot(STAT1_FC, aes(x=tx_symbol,y=log2FoldChange,fill=group,color=group))+
  facet_wrap(~gene_symbol)+
  geom_bar(stat='identity',width = 0.7,position='dodge')+xlab("")  + ylab(expression(paste(log[2],'(FoldChange)'))) +theme_bw()+
  theme(strip.text=element_text(face  = 'italic'),legend.position = c(0.88,0.88),legend.margin =  unit(0.001, 'cm'),legend.key.size = unit(0.2, 'cm'), legend.key.height = unit(0.2, 'cm'),
        legend.key.width = unit(0.2, 'cm'),legend.title = element_blank(), text=element_text(size = 7),axis.text.x = element_text( face="italic", angle = 30, vjust = 0.5))+
  scale_fill_manual(values=pal_jco("default")(4)[c(1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(1,2,4)]) #scale_fill_jco()+scale_color_jco()+
  

########################  sQTL #######################################

# 影响sqtl最显著的是rs118149197
anno_s1_stat1_stat2 <- read.table('./data/stat1_rs_qtl.txt',header = T)
SCD <- rownames(coldata)[coldata$Diagnosis=='SCD']
AD <- rownames(coldata)[coldata$Diagnosis=='AD']
MCI <- rownames(coldata)[coldata$Diagnosis=='MCI']
NC <- rownames(coldata)[coldata$Diagnosis=='NC']

n_SCD <- apply(anno_s1_stat1_stat2[,SCD], 1, function(x){xf<- factor(x, levels=c('0/0','0/1','1/1','./.')); table(xf)['0/1'] + table(xf)['1/1']})
n_AD <- apply(anno_s1_stat1_stat2[,AD], 1, function(x){xf<- factor(x, levels=c('0/0','0/1','1/1','./.')); table(xf)['0/1'] + table(xf)['1/1']})
n_MCI <- apply(anno_s1_stat1_stat2[,MCI], 1, function(x){xf<- factor(x, levels=c('0/0','0/1','1/1','./.')); table(xf)['0/1'] + table(xf)['1/1']})
n_NC <- apply(anno_s1_stat1_stat2[,NC], 1, function(x){xf<- factor(x, levels=c('0/0','0/1','1/1','./.')); table(xf)['0/1'] + table(xf)['1/1']})

n_mut <- as.data.frame(cbind(n_NC,n_SCD,n_MCI,n_AD))
n_mut[is.na(n_mut)] <- 0
n_mut$id <- anno_s1_stat1_stat2$ID

#load('./data/anno_s1_stat1_stat2.RData')
#load('./data/anno_stat1.RData')


stat1_cluster <- read.table('/Users/songlt/Documents/fastqtl/stat1.splice.txt')
stat1_cluster <- subset(stat1_cluster,V1=='chr2:191874730:191878744:clu_4134_NA')
stat1_cluster$FDR <- p.adjust(stat1_cluster$V4)
stat1_cluster <- subset(stat1_cluster, FDR < 0.05)
rownames(stat1_cluster) <- str_split_fixed(stat1_cluster$V2, '_', 5)[,5]

Table_s5 <- cbind(n_mut, anno_s1_stat1_stat2[, c(1:7,12,13,138,139)] )
Table_s5 <- cbind(stat1_cluster[Table_s5$id,c('V4','FDR')], Table_s5 )
colnames(Table_s5)[1] <- 'P-value'

#write.table(Table_s5[,-7], file='./figure/tables/TableS5.txt',sep='\t',quote=F)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Figure 3F: rs118149197  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

STAT1_clu$rs118149197 <- as.character(subset(anno_s1_stat1_stat2,ID=='rs118149197')[,rownames(STAT1_clu)])#anno_s1_stat1_stat2
STAT1_clu$rs118149197[STAT1_clu$rs118149197=='0/1'] <- 'C/G'
STAT1_clu$rs118149197[STAT1_clu$rs118149197=='0/0'] <- 'C/C'


STAT1_clu_sqtl_rs118149197 <- ggplot(STAT1_clu, aes(x=rs118149197,y=value,fill=rs118149197, color=rs118149197)) +
  geom_boxplot(alpha=0.2,width=0.5) +theme_bw() +theme(legend.position = '',text = element_text(size = 7),strip.text=element_text(face  = 'italic') )+ 
  ylab('sample level PSI (%)')+xlab('')+
  annotate("text",label=bquote(italic('P')~.(sprintf("= %.2g",stat1_cluster[stat1_cluster$V2=="chr2_191878389_C_G_rs118149197","V4"]))),
           x=1.5,y=21, color=1,size=2.5) +
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2)])#+ #scale_fill_jco()+scale_color_jco()+

STAT1_clu_sqtl_rs118149197 <- ggplot(STAT1_clu, aes(x=rs118149197,y=value,fill=rs118149197, color=rs118149197)) +
  geom_boxplot(alpha=0.2,width=0.5) +theme_bw() +theme(legend.position = '',text = element_text(size = 7),strip.text=element_text(face  = 'italic') )+ 
  ylab('sample level PSI (%)')+xlab('')+
  ylim(c(0,22))+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2)])+ #scale_fill_jco()+scale_color_jco()+
  geom_signif(test='wilcox.test', map_signif_level=T,size = 0.4, textsize = 2,
            comparisons = list(c("C/C", "C/G") )) 


STAT1_clu$rs12621178 <- as.character(subset(anno_s1_stat1_stat2,ID=='rs12621178')[,rownames(STAT1_clu)])#anno_s1_stat1_stat2
STAT1_clu$rs117437688 <- as.character(subset(anno_s1_stat1_stat2,ID=='rs117437688')[,rownames(STAT1_clu)])#anno_s1_stat1_stat2
STAT1_clu$rs41474144 <- as.character(subset(anno_s1_stat1_stat2,ID=='rs41474144')[,rownames(STAT1_clu)])#anno_s1_stat1_stat2
STAT1_clu$rs16833157 <- as.character(subset(anno_s1_stat1_stat2,ID=='rs16833157')[,rownames(STAT1_clu)])#anno_s1_stat1_stat2
STAT1_clu$rs11677408 <- as.character(subset(anno_s1_stat1_stat2,ID=='rs11677408')[,rownames(STAT1_clu)])#anno_s1_stat1_stat2
STAT1_clu$rs118149197 <- as.character(subset(anno_s1_stat1_stat2,ID=='rs118149197')[,rownames(STAT1_clu)])#anno_s1_stat1_stat2

mSTAT1_clu <- melt(STAT1_clu[,-1],id=c("value" ,"dia" ))
colnames(mSTAT1_clu) <- c('PSI','dia','ID','genotype')

STAT1_clu_sqt <- ggplot(mSTAT1_clu, aes(x=genotype,y=PSI,fill=genotype, color=genotype)) +
  geom_boxplot(alpha=0.2,width=0.5) +theme_bw() +facet_wrap(.~ID)+
  ylab('sample level PSI (%)')+xlab('')+
  theme(legend.position = '')+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2)])#+ #scale_fill_jco()+scale_color_jco()+


dev.print(pdf, file='./figure/figs_STAT1_clu_sqtl.pdf',width=7,height=4)

########################################################################################################
#####                             FIGURES                           #####
########################################################################################################


f2a1_text <-grid::textGrob(bquote("SCD"~intersect(MCI)~italic('P')~.(sprintf(" = %.2g",p=phyper_test1('nc_scd','nc_mci')))  ),gp=gpar(fontsize=7),y = unit(0.9, "npc"))
f2a2_text <-grid::textGrob(bquote("SCD"~intersect(AD)~italic('P')~.(sprintf(" = %.2g",p=phyper_test1('nc_scd','nc_ad')))  ),gp=gpar(fontsize=7),y = unit(0.9, "npc"))
f2a3_text <-grid::textGrob(bquote("MCI"~intersect(AD)~italic('P')~. (sprintf(" = %.2g",p=phyper_test1('nc_mci','nc_ad'))) ),gp=gpar(fontsize=7),y = unit(0.9, "npc"))

f2b <- f2_scd_mrna_splice_venn

f2b1_text <- grid::textGrob(bquote("DEG"~intersect(DET)~italic('P')~.(sprintf(" = %.2g", p=phyper_test('SCD'))[1])),gp=gpar(fontsize=7),y = unit(0.9, "npc"))
f2b2_text <-grid::textGrob(bquote("DEG"~intersect(DS)~italic('P')~.(sprintf(" = %.2g", p=phyper_test('SCD'))[2]) ),gp=gpar(fontsize=7),y = unit(0.9, "npc"))
f2b3_text <-grid::textGrob(bquote("DET"~intersect(DS)~italic('P')~.(sprintf(" = %.2g", p=phyper_test('SCD'))[3]) ),gp=gpar(fontsize=7),y = unit(0.9, "npc"))

phyper_test_interferon

f2c1_text <- grid::textGrob(bquote("DEG"~intersect(DET)~italic('P')~.(sprintf(" = %.2g", p=phyper_test_interferon()[1]))),gp=gpar(fontsize=7),y = unit(0.9, "npc"))
f2c2_text <-grid::textGrob(bquote("DEG"~intersect(DS)~italic('P')~.(sprintf(" = %.2g", p=phyper_test_interferon()[2]) )),gp=gpar(fontsize=7),y = unit(0.9, "npc"))
f2c3_text <-grid::textGrob(bquote("DET"~intersect(DS)~italic('P')~.(sprintf(" = %.2g", phyper_test_interferon()[3]) )),gp=gpar(fontsize=7),y = unit(0.9, "npc"))


plot_grid(plot_grid(plot_grid(f2_scd_mrna_splice_venn,interferon_genes_venn,f2a1_text,f2c1_text,f2a2_text,f2c2_text,f2a3_text,f2c3_text,NA,NA,labels = c('A','B'),label_size = 9,scale = c(0.8,0.8),nrow=5,rel_heights = c(1.5,0.15,0.15,0.15,0.3)) , 
                    plot_grid(interferon_DS_bar,labels = c('C'),scale = 0.9,label_size = 9 ),nrow=1, rel_widths = c(2,1),axis='b'),
          plot_grid(f2_ensgb_stat,f2_STAT1_PSI_box, STAT1_clu_sqtl_rs118149197, labels = c('D','E','F'),label_size = 9,scale=c(1.05,0.9,0.9),rel_widths = c(1,0.3,0.3), nrow=1),
          plot_grid(NA,f2_stat1_fc,labels = c('','G'),label_size = 9,rel_widths = c(1,0.6), nrow=1,scale=c(1,0.9)),
          nrow=3, rel_heights  = c(1,0.9,1))

dev.print(pdf,file='./figure/fig3_color1.pdf', width = 7, height = 6.5)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### Figure S6 differential TRANSCRIPT USAGE (isoform percent) 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

library(limma)
setwd('/Users/songlt/Desktop/AD2/')
load('./data/genes_samples.RData')

iso_perc <- as.matrix(read.table('./data/isoform_perc_matrix.txt',header = T, stringsAsFactors = F, row.names = 1))
colnames(iso_perc) <- gsub('X.home1.GENE_proc.SONGLITING.mRNA.rsem.|.genes.results|.isoforms.results','',colnames(iso_perc))
iso_perc <- iso_perc[rownames(iso_perc)%in%c(expd_genes,expd_trans),rownames(coldata)]

stat1_trans <- rownames(ensg_sym_autosome)[ensg_sym_autosome$gene_symbol%in%'STAT1']
stat1_exp <- iso_perc[intersect(stat1_trans,rownames(iso_perc)),]
rownames(stat1_exp) <- ensg_sym_autosome[rownames(stat1_exp),'tx_symbol']
m_stat1_exp <- melt(stat1_exp)
m_stat1_exp$dia <- coldata[m_stat1_exp$Var2,'Diagnosis']
m_stat1_exp$dia <- factor(m_stat1_exp$dia, levels = c('NC','SCD','MCI','AD'))

## Boxplot of transcript usage (isoform percent) 
f1_stat1_percent_Boxplot <- ggplot(m_stat1_exp,aes(dia,value,group=dia,fill=dia,color=dia))+
  geom_boxplot(alpha=0.2)+  facet_wrap(.~Var1,nrow = 2, scales = 'free_y') +
  stat_summary(fun = mean, geom = "line",aes(group = 1),color='grey')+
  stat_summary(fun = mean, geom = "point", size=2)+
  theme_bw() + ylab('Isoform percent (%)')+xlab('')+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = '', text = element_text(size = 8),strip.text=element_text(face  = 'italic'))+
  geom_signif(test='wilcox.test', map_signif_level=F,step_increase = 0.1, margin_top=0.5,size = 0.4, textsize = 2,
              comparisons = list(c("NC", "SCD"),c("NC", "MCI"),c("NC", "AD") )) 

ggsave(f1_stat1_percent_Boxplot, file='~/Desktop/AD2/figure/stat1_percent_Boxplot.pdf',width =7 ,height =4 )

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### Figure 5B overlap with Raj et al., nature genetics 2018
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
nature_gene <- read.table('/Users/songlt/Desktop/AD2/data/naturegene_splicing.txt',sep='\t',header = T)
nature_gene$gene_id <- gsub(' ','',nature_gene$gene_id)
intersect(nc_ad_c$gene, nature_gene$gene_id)
int_splice <- venn.diagram(list(blood = nc_ad_c$gene,
                                `Raj et al.` = nature_gene$gene_id),
                           fill =  pal_jco("default")(2), cex=0.5,cat.cex = 0.5,alpha = 0.7, filename = NULL,col='white');

plot_grid(int_splice,nrow = 1,scale = 0.8)
#dev.print(pdf,file='~/Desktop/int_splice.pdf')

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
### Figure S5C overlap with Li et al., alzheimer dementia 2019?
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

alz_gene <- read.table('/Users/songlt/Desktop/AD2/data/alz12254-sup-0004-tables3.tsv',sep='\t',header = T)
alz_gene$ens <- str_split_fixed(alz_gene$introns, pattern = '-',5)[,4]
ensg_sym <- ensg_sym_autosome
ensg_sym$ensg <- str_split_fixed(rownames(ensg_sym), '\\.',2)[,1]

int_splice1 <- venn.diagram(list(blood = nc_ad_c$gene,
                                 `Li et al.` = subset(ensg_sym, ensg%in%alz_gene$ens)[,'gene_symbol']),
                            fill =  pal_jco("default")(2), cex=0.5,cat.cex = 0.5,alpha = 0.7, filename = NULL,col='white');

plot_grid(f2_DS_venn, int_splice, int_splice1, nrow = 1,scale = 0.7, labels = c('A','B','C'),label_size = 9)
dev.print(pdf,file='~/Desktop/int_splice2.pdf',width=7, height=3)

intersect(nc_ad_c$gene,subset(ensg_sym, ensg%in%alz_gene$ens)[,'gene_symbol'])
plot_grid(int_splice1,nrow = 1,scale = 0.8)
#dev.print(pdf,file='~/Desktop/int_splice1.pdf')
