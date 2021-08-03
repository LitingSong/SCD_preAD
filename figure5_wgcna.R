# WGCNA
# 2021.4.15

# top 5000 gene ( coding + noncoding ) + mirna  0.15 merge, 2deepsplit,
# top 5000 isoform ( coding + noncoding ) + mirna   0.15 merge, 2deepsplit,

library(WGCNA)
allowWGCNAThreads()
library(nnet)

setwd('~/Desktop/AD2/WGCNA/')
load('~/Desktop/AD2/data/DEs_0412.RData')
load('~/Desktop/AD2/data/deseq2_mirna.RData')
load('~/Desktop/AD2/data/ISG_score.RData')


ensg_sym_autosome <- read.table('~/Desktop/AD2/data/ensg_sym_autosome.txt')
coldata$ISG_score <- ISG_score

coldata <- subset(coldata, Diagnosis%in%c('NC','SCD'))
coldata$Diagnosis <- factor(coldata$Diagnosis, levels = c('NC','SCD','MCI','AD') )
Diagi <- class.ind(coldata$Diagnosis)
coldata <- cbind(coldata,Diagi)
coldata$Gender[coldata$Gender=='F'] <- 0
coldata$Gender[coldata$Gender=='M'] <- 1

coldata <- within(coldata, {
  Education <- 0
  Education[Edu==1] <- 0
  Education[Edu==2] <- 6
  Education[Edu==3] <- 9
  Education[Edu==4] <- 12
  Education[Edu==5] <- 10
  Education[Edu==6] <- 15
  Education[Edu==6] <- 18
})

coldata <- coldata[,c( 'SCD', 'Gender', 'Age','MMSE' ,'MoCA_B','ACEIII_score','Education','apoe4','ISG_score')]
colnames(coldata) <- c( 'SCD', 'Gender', 'Age','MMSE' ,'MoCA_B','ACEIII','Education','APOE4','ISG_score')
# protein_coding 和 lncRNA
cod_nc_rnas <- rownames(ensg_sym_autosome[ensg_sym_autosome$biotype%in%c('lncRNA','protein_coding'),])

expro <- t(rbind(gene_exp$norm_exp[intersect(rownames(gene_exp$norm_exp),cod_nc_rnas ), ], mirna_vst )[,rownames(coldata)])


datExpr = expro[,order(apply(expro,2,mad), decreasing = T)[1:5000]]

sampleTree = hclust(dist(datExpr), method = "average")

powers = c(c(4:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers,networkType='signed',RsquaredCut = 0.8)
sft1 = pickSoftThreshold(datExpr, powerVector = powers,networkType='signed',RsquaredCut = 0.75)

# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.8;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.8,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#dev.print(pdf, file=paste(tissue, '_sft.pdf',sep=''))


if(is.na(sft$powerEstimate) & is.na(sft1$powerEstimate)){power <- 12}
if(is.na(sft$powerEstimate) & !is.na(sft1$powerEstimate)){power <- sft1$powerEstimate}
if(!is.na(sft$powerEstimate) & !is.na(sft1$powerEstimate)){power <- sft$powerEstimate}


minModuleSize = min(20, ncol(datExpr)/2 )
cor <- WGCNA::cor


mergeCutHeight <- 0.15
deepSplit <- 2
#for (deepSplit in c(4)){ 
net <-  blockwiseModules(
  datExpr,networkType = "signed",
  power = power,deepSplit=deepSplit,
  detectCutHeight = 0.995, 
  minModuleSize = min(20, ncol(datExpr)/2 ),
  
  # Gene reassignment, module trimming, and module "significance" criteria
  
  reassignThreshold = 1e-6,
  minCoreKME = 0.5, 
  minCoreKMESize = minModuleSize/3,
  minKMEtoStay = 0.3,
  
  # Module merging options
  
  mergeCutHeight = mergeCutHeight, 
  impute = TRUE, 
  trapErrors = FALSE, 
  
  # Output options
  
  numericLabels = T,
  
  # Options controlling behaviour
  
  useInternalMatrixAlgebra = FALSE,
  useCorOptionsThroughout = TRUE,
  verbose = 0, indent = 0,
)

table(net$colors)

moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
plotDendroAndColors(net$dendrograms[[1]],moduleColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels= FALSE, hang = 0.03,
                    addGuide= TRUE, guideHang = 0.05)

# dev.print(pdf, file='~/Desktop/AD2/figure/wgcna_module_6000.pdf')



modulec <- as.data.frame(cbind(colnames(datExpr),moduleColors ))
modulec <- modulec[order(modulec$moduleColors),]
colnames(modulec) <-c('gene','module')
modulec$symbol <- ensg_sym_autosome[modulec$gene,'gene_symbol']

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

coldata1 <- coldata[rownames(datExpr),]
coldata1$Gender <- car::recode(coldata1$Gender,"'F'=0; 'M'=1")
coldata1$APOE4 <- car::recode(coldata1$APOE4,"c('e2e3','e3e3')=0; else=1")
#coldata1 <- coldata1[,c(1:7,11)]
moduleTraitCor = cor(MEs, coldata1, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


# Plot the relationships among the eigengenes and the trait

#par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(2,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)

MEs_colpheno = orderMEs(cbind(MEs, coldata1))
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap",marDendro = c(2,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)

#dev.print(pdf, file=paste0('plotEigengeneNetworks','_',ngenes,'_',cutHeights,'_',level,'.pdf',sep='' ))

############# Display the correlation values within a heatmap plot\########################################
par(mar=c(3,8,4,2)+0.1)
par(mfrow = c(1,1));
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
textMatrix[moduleTraitPvalue>0.05] <- NA

module_trait_heatmap <- labeledHeatmap(Matrix = moduleTraitCor,xLabels = names(coldata1),yLabels = names(MEs),ySymbols = names(MEs),
               colorLabels = FALSE,colors = greenWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,
               cex.text = 0.4,zlim = c(-1,1),
               main = '',xLabelsAngle = 35,xLabelsAdj = 0.6,font.lab.x=1,cex.lab=0.7)

module_trait_heatmap <- labeledHeatmap(Matrix = moduleTraitCor[,c('SCD','ISG_score','Gender','Age','Education','APOE4')],xLabels =c('SCD','ISG_score','Gender','Age','Education','APOE4') ,
                                       yLabels = paste('M',1:length(names(MEs)),sep=''),ySymbols = gsub('ME','',names(MEs)),
                                       colorLabels = FALSE,colors = blueWhiteRed(50),
                                       textMatrix = textMatrix[, c(1,9,2,3,7,8)],setStdMargins = FALSE,
                                       cex.text = 0.4,zlim = c(-1,1),
                                       main = '',xLabelsAngle = 35,xLabelsAdj = c(0.5,1),font.lab.x=1,cex.lab=0.7)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
modNames = substring(names(MEs), 3)
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");


Diag = as.data.frame(coldata$SCD)
names(Diag) = "Diag"

geneTraitSignificance = as.data.frame(cor(datExpr, Diag, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Diag), sep="")
names(GSPvalue) = paste("p.GS.", names(Diag), sep="")


# identify hub gene

for(module in c('cyan','black')){
  #for(module in c('greenyellow')){
  
  #module = "cyan"
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  ensg_sym_autosome[colnames(datExpr)[moduleGenes], 'gene_symbol']
  PA <- gProfileR::gprofiler(ensg_sym_autosome[colnames(datExpr)[moduleGenes], 'gene_symbol'],src_filter = c('GO:BP','KEGG','REAC'), min_set_size=15,max_set_size=500,hier_filtering = 'moderate')
  
  PA$module <- paste('M',which(rownames(moduleTraitCor)==paste('ME',module,sep='')),sep='')
  
  #PA <- gProfileR::gprofiler(ensg_sym_autosome[colnames(datExpr)[moduleGenes], 'gene_symbol'],max_set_size = 1000,src_filter = c('GO:BP','KEGG','REAC'),significant = F)
  
  #PA <- gProfileR::gprofiler(colnames(datExpr)[moduleGenes],max_set_size = 1000,src_filter = c('GO:BP','KEGG','REAC'),significant = F)
  
  #scd_pat_u%>%arrange(p.value)%>% distinct(.,term.name,.keep_all= TRUE) %>% head(.,8) 
  
  #PA%>%arrange(p.value) %>%head(.,6)
  module_path <- ggplot(PA%>%arrange(p.value) %>%head(.,5) ,aes(x=reorder(term.name, p.value),y=-log10(p.value)))+
    geom_bar(stat = 'identity',width = 0.5,fill=pal_jco("default")(4)[c(1)],color=pal_jco("default")(4)[c(3)])+theme_bw()+#geom_bar(stat = 'identity',width = 0.5,fill=module)+theme_bw()+
    facet_wrap(.~module)+
    theme(legend.position = '',panel.grid.minor  = element_blank(),
          axis.text.x = element_text( size=8), 
          text = element_text(size=9)) + scale_x_discrete(labels=function(x) str_wrap(x, width=40))+
    geom_hline(yintercept = -log10(0.05),color=2,linetype = "dashed" ) +
    xlab('') +    ylab(expression(paste(-log[10],'(',italic('P'),'-value)')))+coord_flip()
  
  #dev.print(pdf, file=paste0('~/Desktop/AD2/figure/',module,'_iso_pathw.pdf'))
  ggsave(module_path, file=paste0('~/Desktop/AD2/figure/',module,'_gene_pathw.pdf'),width = 3,height = 3.1 )
}


sizeGrWindow(7, 7);
par(mfrow = c(1,1));

# Scatterplot

module <- 'cyan'
column = match(module, modNames);
moduleGenes = moduleColors==module;
ensg_sym_autosome[colnames(datExpr)[moduleGenes], 'gene_symbol']

geneModuleMembership$gene <- ensg_sym_autosome[rownames(geneModuleMembership),'gene_symbol']
geneTraitSignificance$gene <- ensg_sym_autosome[rownames(geneTraitSignificance),'gene_symbol']

verboseScatterplot((geneModuleMembership[moduleGenes, column]),
                   (geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for SCD", abline = TRUE,abline.color = module, abline.lty = 2,
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1, cex.lab = 1, cex.axis = 1,pch=20, col = module,ylim=c(-0.25,0.1))

geneModuleTrait <- cbind(geneModuleMembership,geneTraitSignificance)
hub_genes <- geneModuleTrait[as.vector(abs(geneModuleTrait$GS.Diag)> 0.15 & abs(geneModuleTrait[,paste0('MM',module)])>0.8),c('GS.Diag',paste0('MM',module),'gene')]
hub_genes <- geneModuleTrait[as.vector(abs(geneModuleTrait$GS.Diag)> 0.2 & abs(geneModuleTrait[,paste0('MM',module)])>0.8),c('GS.Diag',paste0('MM',module),'gene')]

save(hub_genes, file='~/Desktop/AD2/data/hub_genes.RData')
for (i in 1:nrow(hub_gene)){
  text(y = hub_gene[i,1]-0.005, x = hub_gene[i,2], hub_gene[i,3],cex = 0.5,col = 'red')
}

cyan_module_scater <- as.data.frame(cbind(geneModuleMembership[moduleGenes, column],geneTraitSignificance[moduleGenes, 1]))
rownames(cyan_module_scater) <- rownames(geneModuleMembership)[moduleGenes]
colnames(cyan_module_scater) <- c('geneModuleMember','geneTraitSignific')
cyan_module_scater$gene_symbol <- ensg_sym_autosome[rownames(cyan_module_scater),'tx_symbol']

filter_data <- subset(cyan_module_scater,abs(geneModuleMember)> 0.8 & abs(geneTraitSignific )>0.2)

f5_gene_module_trait_scater <- ggscatter(cyan_module_scater, x = 'geneModuleMember', y = 'geneTraitSignific',
                                add = "reg.line", color =module, star.plot.lty=2, # Add regressin line
                                add.params = list(color =module, lty=2,size=0.5), 
                                conf.int = T,size=1) + stat_cor(method = "pearson",size=2.5) +
  labs(x=paste("Module Membership in", module, "module"),y="Gene significance for SCD")+theme_bw()+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = c(0.8,0.2),legend.margin =  unit(0.03, 'cm'),legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank(), legend.text = element_text(size=6),
        text = element_text(size = 7),panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  geom_hline(yintercept = -0.2,linetype = "dashed",size = 0.4,color = 'grey') +
  geom_vline(xintercept = 0.8,linetype = "dashed",size = 0.4,color = 'grey') +
  ggrepel::geom_text_repel(
    aes(label = gene_symbol),fontface  = 'italic',
    color =2,data = filter_data,show.legend = F,segment.color = 'grey',segment.size=0.3, size = 1.5)


f5_gene_module_trait_scater <- ggscatter(cyan_module_scater, x = 'geneModuleMember', y = 'geneTraitSignific',
                                         add = "reg.line", color =module, star.plot.lty=2, # Add regressin line
                                         add.params = list(color =module, lty=2,size=0.5), 
                                         conf.int = T,size=1) + 
  annotate("text",label=expression(paste(italic('R'),' = ',-0.57,', ',italic('P'),' = ',5.4e-09)),x=0.55,y=-0.04, color=1,size=2) +theme_bw()+
  labs(x=paste("Module Membership in", module, "module"),y="Gene significance for SCD")+theme_bw()+
  scale_fill_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+scale_color_manual(values=pal_jco("default")(4)[c(3,1,2,4)])+#scale_fill_jco()+scale_color_jco()+
  theme(legend.position = c(0.8,0.2),legend.margin =  unit(0.03, 'cm'),legend.key.size = unit(0.3, 'cm'), legend.key.height = unit(0.3, 'cm'),
        legend.key.width = unit(0.3, 'cm'),legend.title = element_blank(), legend.text = element_text(size=6),
        text = element_text(size = 7),panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
  geom_hline(yintercept = -0.2,linetype = "dashed",size = 0.4,color = 'grey') +
  geom_vline(xintercept = 0.8,linetype = "dashed",size = 0.4,color = 'grey') +
  ggrepel::geom_text_repel(
    aes(label = gene_symbol),fontface  = 'italic',
    color =2,data = filter_data,show.legend = F,segment.color = 'grey',segment.size=0.3, size = 1.5)



ggsave(f5_gene_module_trait_scater, file=paste0('~/Desktop/AD2/figure/',module,'verboseScatterplot.pdf'),width = 3,height = 2.7)
dev.print(pdf,file=paste0('~/Desktop/AD2/figure/',module,'verboseScatterplot1.pdf') )
#############################################################################################################
# 5. export to cytoscape。
#############################################################################################################

#load(net$TOMFiles[1], verbose=F)
# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(datExpr, power =power , networkType='signed'); 
# 提取指定模块的基因名
# Order modules by their significance for weight
# Add module membership information in the chosen order
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
ensg_sym <- ensg_sym_autosome
probes = colnames(datExpr) ## 
inModule = (moduleColors==module);
modGene = ensg_sym[probes[inModule], 'gene_symbol']; #### 
modGene = probes[inModule]

gene_Attr <- data.frame(row.names =colnames(datExpr) )
gene_Attr$name <- colnames(datExpr)
gene_Attr$color <- moduleColors
gene_Attr$type[!grepl(pattern = 'ENS',rownames(gene_Attr))] <- 'miRNA'
gene_Attr$type[grepl(pattern = 'ENST',rownames(gene_Attr))] <- 'Isoform'
gene_Attr$type[grepl(pattern = 'ENSG',rownames(gene_Attr))] <- 'Gene'
gene_Attr$type[rownames(gene_Attr)%in% rownames(ensg_sym)[ensg_sym$biotype%in%c('lncRNA','ncRNA')]] <- 'ncRNA'
gene_Attr$hub <- 'N'
gene_Attr$hub[rownames(gene_Attr)%in% rownames(hub_genes)] <- 'Y'
gene_Attr$name[gene_Attr$name%in%rownames(ensg_sym)] <- ensg_sym[gene_Attr$name[gene_Attr$name%in%rownames(ensg_sym)],'gene_symbol']

modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modGene, modGene)
dimnames(modTOM) = list(modGene, modGene)

cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), "0.1.txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), "0.1.txt", sep=""),
  weighted = TRUE,
  threshold = 0.15,
  nodeNames = modGene, 
  nodeAttr = gene_Attr[inModule,]
)

moduleColors[inModule]
moduleColors[inModule]



#hubGenes <- chooseTopHubInEachModule(datExpr,moduleColors)
#write.table (HubGenes,file = "HubGenes_of_each_module.xls",quote=F,sep='\t')

