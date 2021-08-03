## 1. Preprocess: to remove outliers and individuals with comorbidities
## 2. Demographic characteristics
## 2021.01.14
## Song Liting

library(DESeq2)
library(VennDiagram)
library(dplyr)
library(reshape2)
library(ggplot2)
library(clusterProfiler)
library(qvalue)
library(VennDiagram)
library(stringr)
library(haven)
library(ggsci)
options(stringsAsFactors = F)


#s1: Obtain clinical information of samples, and remove samples with hematological diseases, tumor history, other brain diseases and neuropsychiatric diseases (Parkinson|Brain|Spirit|Neurology|Tumor)
d1 <- read_sav("~/Documents/AD/AD_CLINICAL/NPT data_fudan_20200708（补录）.sav")
d1 <- as.data.frame(d1)
d2 <- read.table('~/Documents/AD/AD_CLINICAL/d1.txt', sep='\t', header = T)

rem_d1 <- d1$NPT_ID[grepl(pattern = '17|18|19|20|21|22|23|26|27|27' ,as.character(d1$Disease_history))]
rem_d2 <- d2$NPT_ID[grepl(pattern = '17|18|19|20|21|22|23|26|27|27' ,d2$Disease_history)]

rem_disease_sample <- c(rem_d1,rem_d2)

# s2 Read the expression file, remove the samples with the above-mentioned disease history, remove the outlier samples, remove the low expression transcripts
setwd('/Users/songlt/Desktop/AD2/')

# ensg
# ensg_gene <- read.table('~/Documents/database/genes1.txt',sep='\t')[,c(1,2,4,3,4,3)]
# colnames(ensg_gene) <- c('chr','ID','gene_symbol','type','tx_symbol','tx_biotype')
# ensg_tx <- read.table('~/Documents/database/txs.txt',sep='\t')[,c(1,2,4,3,6,5)]
# colnames(ensg_tx) <- c('chr','ID','gene_symbol','type','tx_symbol','tx_biotype')
# ensg_sym <- rbind(ensg_gene, ensg_tx)

ensg_sym <- read.table('/Users/songlt/Documents/database/ensg_genesymbols.txt')
colnames(ensg_sym) <- c('chr','ID','type','gene_symbol','tx_symbol')

ensg_sym <- within(ensg_sym, {
  biotype <- type
  biotype[grepl('pseudogene',biotype)] <- 'pseudogene'
  biotype[grepl('IG',biotype)] <- 'protein_coding'#'IG_gene'
  biotype[biotype%in%c('antisense','lincRNA','lncRNA', 'sense_intronic','sense_overlapping')] <- 'lncRNA'
  biotype[biotype%in%c('miRNA','misc_RNA','rRNA','scRNA','snoRNA','snRNA','vaultRNA')] <- 'ncRNA'
  biotype[grepl('TR',biotype)] <- 'protein_coding' #'TR_gene'
})

ensg_sym_autosome <- subset(ensg_sym, !chr %in% c('chrM','chrY', 'chrX'))
rownames(ensg_sym_autosome) <- ensg_sym_autosome$ID
write.table(ensg_sym_autosome, file='./data/ensg_sym_autosome.txt', sep='\t', quote = F, row.names = T )



# sample data
coldata <- read.table('./data/coldata.txt',header = T, stringsAsFactors = F, row.names = 1)
MD52 <- read.table('./data/b1md5.txt')
#MD52 <- read.table('./data/b2md5.txt')
MD52$sample <- str_split_fixed(MD52$V2,'-',3)[,2]
MD52 <- subset(MD52, sample%in%rownames(coldata))
MD52$sample <- factor(MD52$sample, levels = rownames(coldata))
MD52 <- MD52[order(MD52$sample),]
MD52_R1 <- MD52[grepl('R1.fastq.gz',MD52$V2),]
MD52_R2 <- MD52[grepl('R2.fastq.gz',MD52$V2),]
cbind(MD52_R1, MD52_R2)->a
coldata$group <- 'Case set'
coldata$group[coldata$Diagnosis=='NC'] <- 'Control set'
write.table(coldata,file='~/Desktop/c1.txt',sep='\t',quote = F, row.names = F)

write.table(a, file='./data/mrna_md5.txt',sep='\t',quote = F, row.names = F)

coldata$Age_group <- factor(coldata$Age_group, levels = c('youth','middle_aged','yonger_aged','aged'))
coldata <- coldata[setdiff(rownames(coldata),rem_disease_sample),]

## Remove low expression transcripts
cts <- as.matrix(read.table('./data/gene_tpm_matrix.txt',header = T, stringsAsFactors = F, row.names = 1))
colnames(cts) <- gsub('X.home1.GENE_proc.SONGLITING.mRNA.rsem.|.genes.results|.isoforms.results','',colnames(cts))
cts <- (cts[,rownames(coldata)])
cts <- cts[rowSums((cts) >  1) >=  nrow(coldata)*0.25, ]

cts_trans <- as.matrix(read.table('./data/isoform_tpm_matrix.txt',header = T, stringsAsFactors = F, row.names = 1))
colnames(cts_trans) <- gsub('X.home1.GENE_proc.SONGLITING.mRNA.rsem.|.genes.results|.isoforms.results','',colnames(cts_trans))
cts_trans <- (cts_trans[,rownames(coldata)])
cts_trans <- cts_trans[rowSums((cts_trans) >  1) >=  nrow(coldata)*0.25, ]

expd_genes <- intersect(rownames(cts),rownames(ensg_sym_autosome))
expd_trans <- intersect(rownames(cts_trans),rownames(ensg_sym_autosome))

save(coldata,expd_trans,expd_genes, ensg_sym_autosome, file='./data/genes_samples.RData' )


load('./data/genes_samples.RData')

cts <- log2(cts+0.001)

sampleTree = hclust(dist(t(cts)), method = "average");
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


## Remove outliers
cts <- gene_exp$norm_exp
k = colSums(cor(cts))
Z.k = scale(k)
outliers <- colnames(cts)[(which(Z.k< -2))]

pca <- prcomp(t(cts), scale. = T, rank. = 2)
outliers <- apply(pca$x, 2, function(x) which( abs(x - mean(x)) > (3 * sd(x)) ))
outliers_exp <- names(c(outliers$PC1,outliers$PC2))
fviz_pca_ind(prcomp(t(cts), scale. = T, rank. = 2),  axes = c(1, 2) ,habillage = coldata$Diagnosis)


#coldata <- coldata[setdiff(rownames(coldata),outliers_exp),]

expd_genes <- intersect(rownames(cts),rownames(ensg_sym_autosome))
expd_trans <- intersect(rownames(cts_trans),rownames(ensg_sym_autosome))

save(coldata,expd_trans,expd_genes, ensg_sym_autosome, file='./data/genes_samples.RData' )

# S3 Clinical Information Statistics

# f1 How many samples for each group
coldata$Diagnosis <- factor(coldata$Diagnosis, levels=c('NC','SCD','MCI','AD'))
sample_sta <-as.data.frame(table(coldata$Diagnosis))

ggplot(sample_sta , aes(x=Var1,y=Freq))+
  geom_bar(stat = 'identity', aes(fill = Var1))+  ylab("Number of samples") +
  labs(x="Diagnosis", fill='Diagnosis') +
  geom_text(data=sample_sta,aes(x=Var1,y=Freq,label=Freq),vjust=0) + theme_bw() +
  scale_fill_jco()+ theme(legend.position='')#text=element_text( family="Luminari")

dev.print(pdf, file='./figure/s1_sample_count_bar.pdf')


t1 <- c('# samples',sample_sta[,2])
# f2 age distribution
t2.1 <- coldata %>% group_by(Diagnosis) %>% summarise(mean(Age))
t2.2 <- coldata %>% group_by(Diagnosis) %>% summarise(sd(Age))
fit <- aov(Age~Diagnosis, data = coldata)
summary(fit)$p.value
TukeyHSD(fit, conf.level = 0.95)
t2 <- c('Age (Mean (sd))',paste0(round(t2.1$`mean(Age)`,2),' (', round(t2.2$`sd(Age)`,2),')'), summary(fit)[[1]][5][1,1])


ggplot(coldata , aes(x=Diagnosis,y=Age,color=Diagnosis,fill=Diagnosis))+
  geom_boxplot(alpha=0.5) + theme_bw() +scale_color_jco()+
  scale_fill_jco()+ theme(legend.position='')

dev.print(pdf, file='./figure/s1_age_box_allsample.pdf')

#1=文盲；2=小学；3=初中；4=高中；5＝专科（中专/大专）；6=大学本科；7=研究生
coldata <- within(coldata, {
  EDU_year <- 0
  EDU_year[Edu==1] <- 0
  EDU_year[Edu==2] <- 6
  EDU_year[Edu==3] <- 9
  EDU_year[Edu==4] <- 12
  EDU_year[Edu==5] <- 10
  EDU_year[Edu==6] <- 15
  EDU_year[Edu==6] <- 18
})

t2.1 <- coldata %>% group_by(Diagnosis) %>% summarise(mean(EDU_year,na.rm = TRUE))
t2.2 <- coldata %>% group_by(Diagnosis) %>% summarise(sd(EDU_year,na.rm = TRUE))
fit <- aov(EDU_year~Diagnosis, data = coldata)
summary(fit)
TukeyHSD(fit, conf.level = 0.95)
t_edu <- c('EDU_year (Mean (sd))',paste0(round(t2.1$`mean(EDU_year, na.rm = TRUE)`,2),' (', round(t2.2$`sd(EDU_year, na.rm = TRUE)`,2),')'), summary(fit)[[1]][5][1,1])


# f3: Distribution of age group
d2 <- coldata %>%
  group_by(Diagnosis,Age_group) %>%
  summarise(N=n()) %>%
  group_by(Diagnosis) %>%
  mutate(freq = (N/sum(N)*100))

ggplot(d2, aes(x = factor(Diagnosis), y = freq, fill = factor(Age_group))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "Diagnosis", y = "Percent", fill = "Age_group") +
  theme_minimal(base_size = 14)+theme_bw()+
  scale_fill_jco()+ theme()

dev.print(pdf, file='./figure/s1_agegroup_percent.pdf')

# RIN
RIN_f <- ggplot(coldata , aes(x=Diagnosis,y=RIN,color=Diagnosis,fill=Diagnosis))+
  geom_boxplot(alpha=0.5,width=0.5) + theme_bw() +scale_color_jco()+
  scale_fill_jco()+ xlab('')+ylab('RIN (RNA integrity number)')+
  theme(legend.position = '',text = element_text(size = 8),axis.title.y = element_text(size=8))
  

t3.1 <- coldata %>% group_by(Diagnosis) %>% summarise(mean(RIN))
t3.2 <- coldata %>% group_by(Diagnosis) %>% summarise(sd(RIN))
fit <- aov(RIN~Diagnosis, data = coldata)
summary(fit)$p.value
TukeyHSD(fit, conf.level = 0.95)
t3 <- c('RIN',paste0(round(t3.1$`mean(RIN)`,2),' (', round(t3.2$`sd(RIN)`,2),')'), summary(fit)[[1]][5][1,1])


dev.print(pdf, file='./figure/s1_RIN_box_allsample.pdf')


# f4 Gender
d2 <- coldata %>%
  group_by(Diagnosis,Gender) %>%
  summarise(N=n()) %>%
  group_by(Diagnosis) %>%
  mutate(freq = (N/sum(N)*100))

t4 <- c('Gender (M %)',paste0(round(d2$freq[c(2,4,6,8)],2),'%'), chisq.test(matrix(d2$N,nrow=2))$p.value)

ggplot(d2, aes(x = factor(Diagnosis), y = freq, fill = factor(Gender))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "Diagnosis", y = "Percent", fill = "Gender") +
  theme_minimal(base_size = 14)+theme_bw()+
  scale_fill_jco()+ theme()


dev.print(pdf, file='./figure/s1_sex_percent.pdf')

# Apoe genotypes
coldata$apoe_n[coldata$apoe_n=='-1'] <- 0
d2 <- coldata %>%
  group_by(Diagnosis,apoe_n) %>%
  summarise(N=n()) %>%
  group_by(Diagnosis) %>%
  mutate(freq = (N/sum(N)*100))

t4 <- c('APOE (+ %)',paste0(round(d2$freq[c(2,4,6,8)],2),'%'), chisq.test(matrix(d2$N,nrow=2))$p.value)

ggplot(d2, aes(x = factor(Diagnosis), y = freq, fill = factor(Gender))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "Diagnosis", y = "Percent", fill = "Gender") +
  theme_minimal(base_size = 14)+theme_bw()+
  scale_fill_jco()+ theme()


dev.print(pdf, file='./figure/s1_sex_percent.pdf')


d2 <- coldata %>%
  group_by(Diagnosis,Age_group) %>%
  summarise(N=n()) %>%
  group_by(Diagnosis) %>%
  mutate(freq = (N/sum(N)*100))

ggplot(d2, aes(x = factor(Diagnosis), y = freq, fill = factor(Age_group))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "Diagnosis", y = "Percent", fill = "Age_group") +
  theme_minimal(base_size = 14)+theme_bw()+
  scale_fill_jco()+ theme()



# f5 Education
coldata$Edu <- as.character(coldata$Edu)
coldata <- coldata[!is.na(coldata$Edu),]
d3 <- coldata %>% 
  group_by(Diagnosis,Edu) %>% 
  summarise(N=n())%>%
  group_by(Diagnosis) %>%
  mutate(freq = (N/sum(N)*100))

ggplot(d3, aes(x = factor(Diagnosis), y = freq, fill = factor(Edu))) +
  geom_bar(stat="identity", width = 0.7) +
  labs(x = "Diagnosis", y = "Percent", fill = "Edu") +
  theme_minimal(base_size = 14)+theme_bw()+
  scale_fill_jco()+ theme()


dev.print(pdf, file='./figure/s1_edu_percent.pdf')

# f6 mapped-reads

read_map <- read.table( './data/readmap.Log.final.out',sep='~')
read_map <- read_map[c(6,9,24,29,31,33),]
read_map <- as.data.frame(t(gsub('\\D','',as.matrix(read_map))))
colnames(read_map) <- c('all','Uniquely\nmappped','Multi\nmapped','Unmapped\nMismatch','Unmapped\nShort','Unmapped\nOther')
read_map <- as.data.frame(apply(read_map, 2, as.numeric))
read_map$umap <- read_map[,1] - read_map[,2] - read_map[,3]


sample <- read.table('./data/sample.reads.log')
sample <- gsub('Log.final.out','',sample$V1)
rownames(read_map) <- sample

read_map$Diagnosis <- coldata[rownames(read_map),'Diagnosis']
read_map <- melt(read_map,id.vars = 'Diagnosis')
read_map$value <- as.numeric(read_map$value)/1000000
read_map <- read_map[read_map$Diagnosis%in%c('NC','SCD','MCI','AD'), ]
read_map$Diagnosis <- factor(read_map$Diagnosis, levels=c('NC','SCD','MCI','AD'))

read_map %>% group_by(Diagnosis,variable) %>% summarise(mean(value))->a
read_map %>% group_by(Diagnosis,variable) %>% summarise(sd(value))->b

fit <- aov(value~Diagnosis, data = subset(read_map, variable=='umap'))
summary(fit)$p.value
TukeyHSD(fit, conf.level = 0.95)

read_map_f <- ggplot(subset(read_map,variable!='umap'), aes(x=variable,y=value,color=Diagnosis,fill=Diagnosis ),group=Diagnosis)+
  geom_boxplot(alpha=0.5) + labs(x='',y='Reads (Millions)')  +theme_bw() +  
  theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 0.5),legend.position = c(0.8,0.8),legend.margin =  unit(0.01, 'cm'),legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.2, 'cm'),legend.title = element_blank(), legend.text = element_text(size=8),
        text = element_text(size = 8),axis.title.y = element_text(size=8))+
  scale_fill_jco()+scale_color_jco()

fit <- aov(value~Diagnosis+variable, data = read_map)
summary(fit)$p.value
TukeyHSD(fit, conf.level = 0.95)

plot_grid(RIN_f,read_map_f,rel_widths = c(1,2),labels = c('A','B'),scale = 0.9)

dev.print(pdf, file='./figure/s1_mrna_read_map.pdf',width = 7, height = 3)


all_r <- subset(read_map,variable=='all' )
all_um <- subset(read_map,variable=='Uniquely\nmappped' )

fit <- aov(value~Diagnosis, all_r)
summary(fit)
TukeyHSD(fit, conf.level = 0.95)

fit <- aov(value~Diagnosis, all_um)
summary(fit)
TukeyHSD(fit, conf.level = 0.95)

# f7 Cognitive scale
library(ggpubr)
rec_score <- melt(coldata[,c(1,4:7)],id=c('Diagnosis','Edu'))
ggplot(rec_score, aes(x=variable,y=value,color=Diagnosis,fill=Diagnosis ),group=Diagnosis)+
  geom_boxplot(alpha=0.5) + labs(x='',y='Score')  +theme_bw() + 
  geom_signif(test='wilcox.test', comparisons = list(c("NC", "SCD"),c('NC','MCI') , c("NC", "AD") )
              ,map_signif_level=TRUE,step_increase = 0.1, margin_top=0.5,size = 0.4, textsize = 2)+
  theme(axis.text.x = element_text(angle = 0))+
  scale_fill_jco()+  scale_color_jco()
dev.print(pdf, file='./figure/s1_score.pdf')


t7.1 <- coldata %>% group_by(Diagnosis) %>% summarise(mean(MMSE))
t7.2 <- coldata %>% group_by(Diagnosis) %>% summarise(sd(MMSE))
fit <- aov(MMSE~Diagnosis, data = coldata)
summary(fit)
TukeyHSD(fit, conf.level = 0.95)
t7_mmse <- c('MMSE',paste0(round(t7.1$`mean(MMSE)`,2),' (', round(t7.2$`sd(MMSE)`,2),')'), summary(fit)[[1]][5][1,1])

t7.1 <- coldata %>% group_by(Diagnosis) %>% summarise(mean(MoCA_B))
t7.2 <- coldata %>% group_by(Diagnosis) %>% summarise(sd(MoCA_B))
fit <- aov(MoCA_B~Diagnosis, data = coldata)
summary(fit)
TukeyHSD(fit, conf.level = 0.95)
t7_MoCA_B <- c('MoCA_B',paste0(round(t7.1$`mean(MoCA_B)`,2),' (', round(t7.2$`sd(MoCA_B)`,2),')'), summary(fit)[[1]][5][1,1])


t7.1 <- coldata %>% group_by(Diagnosis) %>% summarise(mean(ACEIII_score))
t7.2 <- coldata %>% group_by(Diagnosis) %>% summarise(sd(ACEIII_score))
fit <- aov(ACEIII_score~Diagnosis, data = coldata)
summary(fit)
TukeyHSD(fit, conf.level = 0.95)
t7_ACEIII_score <- c('ACEIII_score',paste0(round(t7.1$`mean(ACEIII_score)`,2),' (', round(t7.2$`sd(ACEIII_score)`,2),')'), summary(fit)[[1]][5][1,1])


# Table 1
t1
t2
t3
t4
t5
t7_mmse
t7_MoCA_B
t7_ACEIII_score


