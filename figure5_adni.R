# adni conversion analysis
### Song Liting

options(stringsAsFactors = F)
library(stringr)
library(limma)
library(cowplot)
library(VennDiagram)
library(gProfileR)
library(ggsci)
library(ggplot2)
library(survminer)
library(survival)

load('~/Desktop/AD2/data/hub_genes.RData')
adni <- read.csv('/Users/songlt/Documents/AD/ADNI/ADNI_Gene_Expression_Profile.csv', row.names = 1,check.names = F)

# combine the expression values of multiple probes for the same gene (maximum)
label <- as.data.frame(t(adni[1:2,-1:-2]))
label$sample_vis <- paste(label$SubjectID,label$`#Visit`,sep='--')

adni <- adni[-1:-8,-747]
adni <- aggregate( matrix(as.numeric(unlist(adni[,-1:-2])),nrow=nrow(adni)), by = list(adni$`.1`),FUN = max,na.rm = TRUE )
adni1 <- adni

rownames(adni) <- adni$Group.1
adni <- adni[,-1]
adni <- adni[-1,]

colnames(adni) <- label$sample_vis[-745]
adni1 <-  adni
DXSUM_PDXCONV_ADNIALL <- read.csv('/Users/songlt/Documents/AD/ADNI/DXSUM_PDXCONV_ADNIALL-2.csv')
DXSUM_PDXCONV_ADNIALL$sample_vis <- paste(DXSUM_PDXCONV_ADNIALL$PTID,DXSUM_PDXCONV_ADNIALL$VISCODE,sep='--')
DXSUM_PDXCONV_ADNIALL$sample_vis2 <- paste(DXSUM_PDXCONV_ADNIALL$PTID,DXSUM_PDXCONV_ADNIALL$VISCODE2,sep='--')
DXSUM_PDXCONV_ADNIALL$id_v <- paste(DXSUM_PDXCONV_ADNIALL$RID,DXSUM_PDXCONV_ADNIALL$VISCODE,sep='--')
DXSUM_PDXCONV_ADNIALL$id_v2 <- paste(DXSUM_PDXCONV_ADNIALL$RID,DXSUM_PDXCONV_ADNIALL$VISCODE2,sep='--')

DXSUM_PDXCONV_ADNIALL$DX[DXSUM_PDXCONV_ADNIALL$DXCURREN==1|DXSUM_PDXCONV_ADNIALL$DXCHANGE%in%c(1,7,9)] <- 'NC'
DXSUM_PDXCONV_ADNIALL$DX[DXSUM_PDXCONV_ADNIALL$DXCURREN==2|DXSUM_PDXCONV_ADNIALL$DXCHANGE%in%c(2,4,8)] <- 'MCI'
DXSUM_PDXCONV_ADNIALL$DX[DXSUM_PDXCONV_ADNIALL$DXCURREN==3|DXSUM_PDXCONV_ADNIALL$DXCHANGE%in%c(3,5,6)] <- 'AD'


DXSUM <- DXSUM_PDXCONV_ADNIALL[,c('RID','sample_vis','sample_vis2','id_v','id_v2','DX','VISCODE','VISCODE2')]
label_DXSUM <-  merge(label,DXSUM,by='sample_vis', all.x=T)

label_DXSUM$DX[label_DXSUM$sample_vis=="006_S_4192--v04"] <- 'AD'
label_DXSUM$DX[label_DXSUM$sample_vis=="114_S_4379--v04"] <- 'AD'
label_DXSUM$DX[label_DXSUM$sample_vis=="073_S_2182--m03"] <- 'MCI'
label_DXSUM$DX[label_DXSUM$sample_vis=="009_S_4564--v02"] <- 'MCI'
label_DXSUM$DX[label_DXSUM$sample_vis=="094_S_4282--v04"] <- 'AD'
label_DXSUM$DX[label_DXSUM$sample_vis=="072_S_4102--v04"] <- 'MCI'
label_DXSUM$DX[label_DXSUM$sample_vis=="072_S_4103--v04"] <- 'NC'
label_DXSUM$DX[label_DXSUM$sample_vis=="037_S_4432--v04"] <- 'MCI'
label_DXSUM$DX[label_DXSUM$sample_vis=="013_S_4268--v04"] <- 'MCI'
label_DXSUM$DX[label_DXSUM$sample_vis=="022_S_4266--v04"] <- 'NC'

#label_DXSUM$RID <- as.numeric(str_split_fixed(label_DXSUM$SubjectID,pattern = '_',3)[,3])

adni_ec <- read.csv('/Users/songlt/Documents/AD/ADNI/Neuropsychological/ECOGPT.csv')[,c(1,3,5,6,9)]
adni_ec$id_v <- paste(adni_ec$RID,adni_ec$VISCODE,sep='--')
adni_ec$id_v2 <- paste(adni_ec$RID,adni_ec$VISCODE2,sep='--')


label_DXSUM_ec <-  merge(label_DXSUM,adni_ec,by='id_v',all.x=T)

attach(DXSUM_PDXCONV_ADNIALL)
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==0 & DXCURREN==1] =1
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==0 & DXCURREN==2] =2
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==0 & DXCURREN==3] =3
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==1 & DXCONTYP==1] =4
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==1 & DXCONTYP==3] =5
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==1 & DXCONTYP==2] =6
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==2 & DXREV==1] = 7 
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==2 & DXREV==2] = 8 
DXSUM_PDXCONV_ADNIALL$DXCHANGE[DXCONV==2 & DXREV==3] = 9 
detach(DXSUM_PDXCONV_ADNIALL)

table(DXSUM_PDXCONV_ADNIALL$sample_vis2) [1:3,c(1:7,59)]

label_DXSUM_ec_NC <- subset(label_DXSUM_ec, DX%in%c('NC'))
DXSUM_PDXCONV_ADNIALL_nc <- subset(DXSUM_PDXCONV_ADNIALL,RID%in%label_DXSUM_ec_NC$RID.x)


conv_nc <- merge(label_DXSUM_ec_NC,DXSUM_PDXCONV_ADNIALL[,c(1:15,53,55,57:59)], by.x='RID.x',by.y='RID' )
conv_nc$time_seq <- conv_nc$VISCODE2.x
conv_nc$time_seq[conv_nc$time_seq%in%c('bl','sc')] <- 0
conv_nc$time_seq[grepl('m',conv_nc$time_seq)] <- (gsub('m','',conv_nc$time_seq[grepl('m',conv_nc$time_seq)]))

conv_nc$time_follow <- conv_nc$VISCODE2
conv_nc$time_follow[conv_nc$time_follow%in%c('bl','sc')] <- 0
conv_nc$time_follow[grepl('m',conv_nc$time_follow)] <- (gsub('m','',conv_nc$time_follow[grepl('m',conv_nc$time_follow)]))

conv_nc$time <- as.numeric(conv_nc$time_follow) - as.numeric(conv_nc$time_seq)


# conversion
conv_nc_Y <- subset(conv_nc,DXCHANGE%in%c(4,6) & time >0 )
conv_nc_Y <- conv_nc_Y%>%arrange(by=time)%>% distinct(.,RID.x,.keep_all= TRUE) # 
conv_nc_Y$status <- 1

# stable
conv_nc_N <- subset(conv_nc,!RID.x%in%conv_nc_Y$RID.x ) 
conv_nc_N$time[conv_nc_N$time>72] <- 72
conv_nc_N <- conv_nc_N%>%arrange(by=-time)%>% distinct(.,RID.x,.keep_all= TRUE) #


conv_data <- rbind(conv_nc_Y,conv_nc_N )
rownames(conv_data) <- conv_data$sample_vis.x

adni1 <- read.csv('/Users/songlt/Documents/AD/ADNI/ADNI_Gene_Expression_Profile.csv', row.names = 1,check.names = F)
adni <- adni1
label <- as.data.frame(t(adni[1:2,-1:-2]))
label$sample_vis <- paste(label$SubjectID,label$`#Visit`,sep='--')

adni <- adni[adni[,2] %in% c(hub_genes$gene,ISG_gene),-747]

adni <- aggregate( matrix(as.numeric(unlist(adni[,-1:-2])),nrow=nrow(adni)), by = list(adni$`.1`),FUN = max,na.rm = TRUE )


rownames(adni) <- adni$Group.1
adni <- adni[,-1]
adni <- adni[-1,]

colnames(adni) <- label$sample_vis[-745]

conv_exp <- cbind(conv_data,t(adni[intersect(rownames(adni),c(hub_genes$gene,ISG_gene)),rownames(conv_data)]))
conv_exp$ISG <- rowMeans(conv_exp[,ISG_gene])
surv_p <- c()


for( n in which(colnames(conv_exp)%in%c(hub_genes$gene))){
  #for( n in which(colnames(conv_exp)%in%c('STAT2'))){
  conv_exp$level <- ifelse(conv_exp[,n] >= median(conv_exp[,n])[1], 'high', 'low' )
  fit <- survfit(Surv(time, status) ~level, data =conv_exp)
  p <- signif(surv_pvalue(fit)[1,2],2)
  surv_p <- as.data.frame(rbind(surv_p, c(colnames(conv_exp)[n], p, paste(xtabs(~level+status,conv_exp)[4], (xtabs(~level+status,conv_exp)[2]+xtabs(~level+status,conv_exp)[4]), sep='/'), paste(xtabs(~level+status,conv_exp)[3], (xtabs(~level+status,conv_exp)[1]+xtabs(~level+status,conv_exp)[3]), sep='/'))))
}

surv_p$fdr <- signif(p.adjust(surv_p$V2,method = 'fdr'),2)
colnames(surv_p) <- c('hub gene','P-value','#conversion / # < median','#conversion / # >= median','FDR')
write.table(surv_p[order(surv_p$`P-value`),], file='~/Desktop/AD2/figure/tables/Table S8.txt',sep='\t',quote = F,row.names = F)

#####################      STAT1      ##############################################

n <- which(colnames(conv_exp)==c('STAT1'))
conv_exp$level <- ifelse(conv_exp[,n] >= quantile(conv_exp[,n])[3], '>= median', '< median' )
fit <- survfit(Surv(time, status) ~level, data =conv_exp)
p_value <- signif(surv_pvalue(fit)[1,2],2)

STAT1_suv_plot <- ggsurvplot(fit, data = conv_exp,  conf.int = TRUE, pval = F, risk.table=F,
                             add.all = F,xlab = "Time (months)",ylab = "Stable probability",legend.labs = c('< median', ">= median"),
                             legend.title = expression(italic('STAT1')),  size=0.5,
                             break.x.by = 12,palette = "jco",ylim = c(0.5, 1) ,xlim=c(0,72),ggtheme=theme_bw())

conv_exp$status_c <- 'Stable'
conv_exp$status_c[conv_exp$status==1] <- 'Convert'

stable.p <- ggtexttable(xtabs(~level+status_c,conv_exp ),rows = NULL,theme=ttheme(
  base_size = 6))%>%
  table_cell_bg(fill = pal_jco("default")(4)[c(1)],  row =2, column=1:2 ,alpha=0.3,color='white')%>%
  table_cell_bg(fill = pal_jco("default")(4)[c(2)],  row =3, column=1:2 ,alpha=0.3,color='white')

stat1_KM <- STAT1_suv_plot$plot + #annotation_custom( ggplotGrob(stable.p), xmin = 10, ymin = 0.5,ymax = 0.7,xmax = 20 )+
  annotate("text",label=expression(paste(italic('P'), " = 0.0031")),x=12,y=0.6, color=1,size=2.5) +# bquote(italic('P') ~ " = " ~ .(p_value)),x=45,y=0.75, color=1,size=2) +
  theme(panel.grid.minor = element_blank(),text = element_text(size = 8))

#dev.print(pdf, file='~/Desktop/AD2/figure/stat1_surv.pdf',width = 4,height=4)


#####################      ISG      ##############################################

n <- which(colnames(conv_exp)==c('ISG'))
conv_exp$level <- ifelse(conv_exp[,n] >= quantile(conv_exp[,n])[4], '>= Q1', '< Q1' )
fit <- survfit(Surv(time, status) ~level, data =conv_exp)
p_value <- signif(surv_pvalue(fit)[1,2],2)

ISG_suv_plot <- ggsurvplot(fit, data = conv_exp,  conf.int = TRUE, pval = F, risk.table=F,
                             add.all = F,xlab = "Time (months)",ylab = "Stable probability",legend.labs = c('< Q1', ">= Q1"),
                             legend.title = 'ISG score', size=0.5,
                             break.x.by = 12,palette = "jco",ylim = c(0.5, 1) ,xlim=c(0,72),ggtheme=theme_bw())

# conversation table
conv_exp$status_c <- 'Stable'
conv_exp$status_c[conv_exp$status==1] <- 'Convert'

stable.p <- ggtexttable(xtabs(~level+status_c,conv_exp ),rows = NULL,theme=ttheme(
  base_size = 6))%>%
  table_cell_bg(fill = pal_jco("default")(4)[c(1)],  row =2, column=1:2 ,alpha=0.3,color='white')%>%
  table_cell_bg(fill = pal_jco("default")(4)[c(2)],  row =3, column=1:2 ,alpha=0.3,color='white')

print(p_value)

ISG_KM <- ISG_suv_plot$plot + #annotation_custom( ggplotGrob(stable.p), xmin = 10, ymin = 0.5,ymax = 0.7,xmax = 20 )+
  annotate("text",label=expression(paste(italic('P'), " = 0.062")),x=12,y=0.6, color=1,size=2.5) +#bquote(italic('P') ~ " = " ~ .(p_value)),x=45,y=0.75, color=1,size=2)+ # 
  theme(text = element_text(size = 8),panel.grid.minor = element_blank())

#dev.print(pdf, file='~/Desktop/AD2/figure/ISG_suv_plot.pdf',width = 4,height=4)


#####################      TRIM22      ##############################################

n <- which(colnames(conv_exp)==c('TRIM22'))
conv_exp$level <- ifelse(conv_exp[,n] >= quantile(conv_exp[,n])[3], '>= median', '< median' )
fit <- survfit(Surv(time, status) ~level, data =conv_exp)
p_value <- signif(surv_pvalue(fit)[1,2],2)

TRIM22_suv_plot <- ggsurvplot(fit, data = conv_exp,  conf.int = TRUE, pval = F, risk.table=F,
                             add.all = F,xlab = "Time (months)",ylab = "Stable probability",legend.labs = c('< median', ">= median"),
                             legend.title = expression(italic('TRIM22')), size=0.5,
                             break.x.by = 12,palette = "jco",ylim = c(0.5, 1) ,xlim=c(0,72),ggtheme=theme_bw())

conv_exp$status_c <- 'Stable'
conv_exp$status_c[conv_exp$status==1] <- 'Convert'

stable.p <- ggtexttable(xtabs(~level+status_c,conv_exp ),rows = NULL,theme=ttheme("default",base_size = 6))%>%
  table_cell_bg(fill = pal_jco("default")(4)[c(1)],  row =2, column=1:2 ,alpha=0.3,color='white')%>%
  table_cell_bg(fill = pal_jco("default")(4)[c(2)],  row =3, column=1:2 ,alpha=0.3,color='white')
print(p_value)

TRIM22_KM <- TRIM22_suv_plot$plot + #annotation_custom( ggplotGrob(stable.p), xmin = 10, ymin = 0.5,ymax = 0.7,xmax = 20 )+
  annotate("text",label=expression(paste(italic('P'), " = 0.0022")),x=12,y=0.6, color=1,size=2.5) +#bquote(italic('P') ~ " = " ~ .(p_value)),x=45,y=0.75, color=1,size=2)+ # 
  theme(text = element_text(size = 8),panel.grid.minor = element_blank())


plot_grid(TRIM22_KM,stat1_KM,ISG_KM,nrow = 1,labels=c('D','E','F'),label_size = 9)

dev.print(pdf,file='./figure/fig5_surv.pdf', width = 7.3, height = 3)

  
  