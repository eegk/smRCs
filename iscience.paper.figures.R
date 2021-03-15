####################################################################################################################
libs<-c("corrplot","pheatmap","","","tidyr","limma","","","" ,"UpSetR",
        "","RColorBrewer","","sva","","fitdistrplus","ff","plyranges",
        "annotables","","GenomicFeatures","ggarrange",
        "dplyr","DESeq2","ggplot2","ggrepel","edgeR","doParallel","ggradar","colorspace",
        "reshape2","rmarkdown","","treemapify","",
        "variancePartition","ggpubr","factoextra","Mfuzz", "universalmotif","ggbio","GenomicRanges",
        "randomForest","","","psych","scales","tibble",
        "colorspace","rpart","tidyverse","","fmsb","colormap","")
lapply(libs, require, character.only = TRUE)
rm(libs)
# registerDoParallel(makeCluster(8)) 
####################################################################################################################
setwd("/Users/gonzae34/Documents/00.exosome/")
####################################################################################################################
### Data * smRCs Prostate Cancer
#load("/Users/gonzae34/Documents/00.exosome/exo_RData/external_smrc_datasets_list.RData")
#prc_ann_counts <- readRDS("/Users/gonzae34/Documents/00.exosome/exo_RData/PRC_ann_counts.RDS") 
#mc_counts <- readRDS("/Users/gonzae34/Documents/00.exosome/exo_RData/mc_counts_bojan.RData") 
### smRCs Liver Cancer
#phen_prc <- readRDS("/Users/gonzae34/Documents/00.exosome/exo_RData/phen2_bojan.RDS")
### setwd("/Users/gonzae34/Documents/00.exosome/iscience.figures/")
####################################################################################################################
gg_color_hue <- function(n) { hues = seq(15, 375, length = n + 1) ; hcl(h = hues, l = 65, c = 100)[1:n]}
####################################################################################################################
#dim(prc_ann_counts)
#ix <- which(mc_counts$Locus %in% prc_ann_counts$Locus)
#mc_counts <- mc_counts[ix,]
### rename & clean
#smrc_datasets <- external_smrc_datasets
#rm(external_smrc_datasets)
### integrate into single object
#smrc_datasets$mssm_pca$smrc_counts <- mc_counts
#smrc_datasets$mssm_pca$ann_smrc_counts <- prc_ann_counts
#rm(mc_counts,prc_ann_counts,ix)
#smrc_datasets$mssm_pca$ann_smrc_counts <- smrc_datasets$mssm_pca$ann_smrc_counts[match(smrc_datasets$mssm_pca$smrc_counts$Locus,smrc_datasets$mssm_pca$ann_smrc_counts$Locus),]
#rownames(smrc_datasets$mssm_pca$smrc_counts) <- smrc_datasets$mssm_pca$smrc_counts$Locus
#smrc_datasets$mssm_pca$smrc_counts <- smrc_datasets$mssm_pca$smrc_counts[,-c(1,2,3)]
#smrc_datasets$mssm_pca$ann_smrc_counts$Locus <- as.character(smrc_datasets$mssm_pca$ann_smrc_counts$Locus)
####################################################################################################################
###
### save.image("/Users/gonzae34/Documents/00.exosome/exo_RData/data_for_iscience_paper.RData")
###
####################################################################################################################
load("/Users/gonzae34/Documents/00.exosome/exo_RData/data_for_iscience_paper.RData")
####################################################################################################################



### Verify integrity
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
identical(rownames(smrc_datasets$mssm_pca$smrc_counts),as.character(smrc_datasets$mssm_pca$ann_smrc_counts$Locus))
colnames(smrc_datasets$mssm_pca$smrc_counts) <- gsub("\\.R1","",colnames(smrc_datasets$mssm_pca$smrc_counts))
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Voom + dgelist
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dim(smrc_datasets$mssm_pca$smrc_counts)
###
smrc_deglist_obj <- DGEList(counts=smrc_datasets$mssm_pca$smrc_counts) 
smrc_deglist_obj <- calcNormFactors(smrc_deglist_obj, method = "TMM") 
smrc_deglist_obj = estimateCommonDisp(smrc_deglist_obj, verbose=TRUE)
### voom
voom_mssm_pca <- voom(smrc_deglist_obj, plot = TRUE, save.plot = TRUE) 
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Principal component analysis
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### check
table(phen_prc$full_sample %in% colnames(voom_mssm_pca$E))

### PCA
mydata <- prcomp( t(voom_mssm_pca$E) , scale=TRUE, center = TRUE )

### SCREE
Scree.plot <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw()
eigen <- get_eigenvalue(mydata)
pca1 <- format(round(eigen$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(eigen$variance.percent[2], 2), nsmall = 2)
Scree.plot ; rm(eigen)

phen_prc$Fluid_Origin <- as.character(phen_prc$Fluid_Cell)
phen_prc$Fluid_Origin[ phen_prc$Fluid_Origin %in% "Serum_Exosome"] <- "Serum exoRNA"
phen_prc$Fluid_Origin[ phen_prc$Fluid_Origin %in% "Urine_Exosome"] <- "Urine exoRNA"
phen_prc$Fluid_Origin[ phen_prc$Fluid_Origin %in% "Normal_Cell"] <- "Adj.Normal Cellular"
phen_prc$Fluid_Origin[ phen_prc$Fluid_Origin %in% "Tumor_Cell"] <- "Tumor Cellular"

prcomp_obj_mssm_pca <- as.data.frame(mydata$x)
prcomp_obj_mssm_pca$full_sample <- rownames(prcomp_obj_mssm_pca)
prcomp_obj_mssm_pca$lib.size <- smrc_deglist_obj$samples$lib.size
prcomp_obj_mssm_pca$norm.factors <- smrc_deglist_obj$samples$norm.factors
prcomp_obj_mssm_pca <- merge(prcomp_obj_mssm_pca, phen_prc, by="full_sample")

pdf(file="iscience.figures/prcomp_mssm_pca.pdf",width = 8,height = 6)
ggplot(data=prcomp_obj_mssm_pca, aes(x=PC1, y=PC2, color = Fluid_Origin)) + 
  geom_point() + 
  geom_label_repel(aes(label=Patient_ID),force = 10, max.iter = 5000,show.legend = FALSE) +
  theme_bw() + 
  labs(x=paste("PC1:",pca1,"%",sep=""), y=paste("PC1:", pca2,"%", sep=""), color="Fluid / Origin") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(4.0)),
        plot.title=element_text(size=rel(4.8)) ) + 
  theme(legend.position="bottom") + 
  guides(color=guide_legend(ncol=2,override.aes = list(size = 5))) 
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Variance Partition
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
phen_prc$Tissue <- phen_prc$Fluid_Cell
phen_prc$Tissue <- gsub("_Exosome","",phen_prc$Tissue); phen_prc$Tissue <- gsub("_Cell","",phen_prc$Tissue)

dummy_phen_prc <- phen_prc
colnames(dummy_phen_prc)[which(colnames(dummy_phen_prc) %in% c("Cell_Type"))] <- "RNA Origin"
colnames(dummy_phen_prc)[which(colnames(dummy_phen_prc) %in% c("technology"))] <- "Technology"
colnames(dummy_phen_prc)[which(colnames(dummy_phen_prc) %in% c("BX_Gleason"))] <- "Gleason"
colnames(dummy_phen_prc)[which(colnames(dummy_phen_prc) %in% c("Initial_PSA"))] <- "Initial PSA"
colnames(dummy_phen_prc)[which(colnames(dummy_phen_prc) %in% c("Patient_ID"))] <- "Patient"

dummy_phen_prc <- dummy_phen_prc[ match( colnames(voom_mssm_pca$E),  as.character(dummy_phen_prc$full_sample) ) , ]
table(as.character(dummy_phen_prc$full_sample) %in% colnames(voom_mssm_pca$E))
identical(as.character(dummy_phen_prc$full_sample),colnames(voom_mssm_pca$E))

### form from phenotype
form <- ~ (1|Patient) + Age + `Initial PSA` + (1|Gleason) + (1|Tissue) + (1|`RNA Origin`)
### run method
variance_exprs_matrix <- fitExtractVarPartModel(voom_mssm_pca$E, form, dummy_phen_prc)

### figure
pdf(file="iscience.figures/varpar_mssm_pca.pdf",width = 5,height = 3.5)
plotVarPart( sortCols( variance_exprs_matrix, decreasing=FALSE) ) + ylim(0,100) + 
  theme(axis.text.x = element_text(angle = 35,hjust = 1)) +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(4.0)),
        #axis.title=element_text(size=rel(4.0)),
        plot.title=element_text(size=rel(4.8)) ) + 
  theme(legend.position="none") + labs(y="Variance %") #+ 
  #guides(color=guide_legend(ncol=2), shape = guide_legend(override.aes = list(size = rel(4.0))))
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Differential expression Cellular to exosomal
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

phen_prc <- phen_prc[ match( colnames(voom_mssm_pca$E),  as.character(phen_prc$full_sample) ) , ]
table(as.character(phen_prc$full_sample) %in% colnames(voom_mssm_pca$E))
identical(as.character(phen_prc$full_sample),colnames(voom_mssm_pca$E))

### design
design = model.matrix( ~ 0 + Cell_Type, data = phen_prc )
colnames(design) <- gsub("Cell_Type","",colnames(design))

### contr.matrix
contr.matrix <- makeContrasts( "Exosomes vs Cells" = Exosome - Cell,  levels = colnames(design))

### statistics
vfit <- lmFit(voom_mssm_pca, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

### results into table
res_de_smRCs_cells_v_exos <- topTable(efit, coef="Exosomes vs Cells", n=Inf, adjust.method="BH") 
### add FDR converted
res_de_smRCs_cells_v_exos$nLogFDR <- -log10(res_de_smRCs_cells_v_exos$adj.P.Val)

rm(vfit,efit,contr.matrix,design,Scree.plot,mydata,pca1,pca2,ix,iy,form)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Volcano plot CElls exosomes
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
volcano.plot<- ggplot(data=res_de_smRCs_cells_v_exos,aes(x=logFC,y=nLogFDR)) + 
  geom_point(data=res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$adj.P.Val>0.05,],color="slategray",alpha=0.7, size = 1) + 
  geom_point(data=res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$logFC > 0 & res_de_smRCs_cells_v_exos$adj.P.Val<0.05,],color="steelblue",size = 1,alpha=0.6) +      
  geom_point(data=res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$logFC < 0 & res_de_smRCs_cells_v_exos$adj.P.Val<0.05,],color="firebrick",size = 1,alpha=0.6) + 
  
  geom_point(data= res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$logFC > 0 & res_de_smRCs_cells_v_exos$adj.P.Val<0.05,], color = "steelblue",size = 1) +
  geom_point(data= res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$logFC < 0 & res_de_smRCs_cells_v_exos$adj.P.Val<0.05,], color = "firebrick",size = 1) +
  
  geom_vline(xintercept = 0, linetype = 2, color = 'black',size=1.5) +
  geom_hline(yintercept = 1.3, linetype = 2, color = 'black',size=1.5) +
  
  xlab("Log2 Fold Change") + ylab("-Log10 FDR") +
  annotate(geom="text", x=-2, y=21, label='atop(bold("Cell"))',color="firebrick",size=10, parse = T) +
  annotate(geom="text", x=3, y=21, label='atop(bold("Exosome"))',color="steelblue",size=10, parse = T) +
  annotate(geom="text", x=-8.5, y=1.5, label='atop(bold("Threshold FDR < 5%"))',color="black",size=6, parse = T) +
  theme_bw() + 
  theme(text = element_text(size=rel(4.8))) 

pdf(file="iscience.figures/volcano_mssm_pca.pdf",width = 8,height = 5)
volcano.plot
dev.off()

rm(volcano.plot)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### MAplot plot CElls exosomes
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

MAPlot <- ggplot(data=res_de_smRCs_cells_v_exos,aes(x=AveExpr,y=logFC)) + 
  geom_point(data=res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$adj.P.Val>0.05,],color="slategray",alpha=0.6, size = 1) + 
  
  geom_point(data=res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$logFC > 0 & res_de_smRCs_cells_v_exos$adj.P.Val<0.05,],color="steelblue",size = 1,alpha=0.6) +      
  geom_point(data=res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$logFC < 0 & res_de_smRCs_cells_v_exos$adj.P.Val<0.05,],color="firebrick",size = 1,alpha=0.6) + 
  geom_point(data= res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$logFC > 0 & res_de_smRCs_cells_v_exos$adj.P.Val<0.05,], color = "steelblue",size = 1) +
  geom_point(data= res_de_smRCs_cells_v_exos[res_de_smRCs_cells_v_exos$logFC < 0 & res_de_smRCs_cells_v_exos$adj.P.Val<0.05,], color = "firebrick",size = 1) +
  
  geom_vline(xintercept = 0, linetype = 2, color = 'black',size=1,5) +
  geom_hline(yintercept = 0, linetype = 2, color = 'black',size=1.5) +
  xlab("Average Expression Log2CPM") + ylab("Log Fold Change") +
  annotate(geom="text", x=11.5, y=-10, label='atop(bold("Cell"))',color="firebrick",size=10, parse = T) +
  annotate(geom="text", x=11.5, y=4, label='atop(bold("Exosome"))',color="steelblue",size=10, parse = T) +
  theme_bw() + 
  theme(text = element_text(size=rel(4.8))) 

pdf(file="iscience.figures/maplot_mssm_pca.pdf",width = 8,height = 5)
MAPlot
dev.off()

rm(MAPlot)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### pHeatmap plot CElls exosomes
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
annotation_col <- data.frame( RNA.Origin = phen_prc$Cell_Type, Tissue = phen_prc$Tissue) ; rownames(annotation_col) <-  phen_prc$full_sample

gg_color_hue(4)

annotation_colors <- list ( RNA.Origin = c( Exosome="steelblue", Cell="firebrick"), High = c( Exosome="steelblue", Cell="firebrick"), 
                            Tissue = c(Normal="#F8766D",Tumor="#7CAE00",Serum="#00BFC4", Urine="#C77CFF") )
ix <- which(rownames(voom_mssm_pca) %in% rownames(res_de_smRCs_cells_v_exos)[which(res_de_smRCs_cells_v_exos$adj.P.Val<0.05)])

annotation_row <- res_de_smRCs_cells_v_exos
annotation_row$High <- "ns"
annotation_row$High[which(annotation_row$adj.P.Val<0.05 & annotation_row$logFC > 0)] <- "Exosome"
annotation_row$High[which(annotation_row$adj.P.Val<0.05 & annotation_row$logFC < 0)] <- "Cell"

iy <- c( which(annotation_row$High %in% "Exosome")[1:1000], which(annotation_row$High %in% "Cell")[1:1000] ) 
#annotation_row <- annotation_row
ix <- which(rownames(voom_mssm_pca$E) %in% rownames(annotation_row)[iy] )

pdf(file="iscience.figures/heatmap_mssm_pca.pdf",width = 5,height = 5.5)
pheatmap(voom_mssm_pca$E[ix,], color = colorRampPalette(c("steelblue", "white","firebrick"))(255),  
         cluster_cols = T, 
         cluster_rows=T, 
         angle_col = 90,
         main ='',
         name='z-score',
         annotation_col = annotation_col,
         annotation_row = annotation_row[,c("High","AveExpr")],
         annotation_colors = annotation_colors,
         border_color="white", 
         show_rownames = FALSE, 
         show_colnames = FALSE, 
         clustering_distance_rows ="euclidean", scale="row")
dev.off()
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### smrc length
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ann_smrc_counts_pca <- smrc_datasets$mssm_pca$ann_smrc_counts

pdf(file="iscience.figures/smrc_length_mssm_pca.pdf",width = 6,height = 4)
ggplot( data= ann_smrc_counts_pca, aes(log10(Length)) ) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity",color="black",fill="white")+
  stat_density(color="steelblue",fill=NA,size=2) + theme_bw() + 
  geom_rug(alpha=0.005,color="black") +
  geom_vline(aes(xintercept = log10(median(Length))), colour = "firebrick", linetype="dashed", size=2) +
  annotate(geom="text", x=3, y=0.75, label="Median", color="firebrick", size=8) +
  theme(text = element_text(size=rel(4.8)))  +
  labs(x="Log10 (smRC Length)", y="Density")
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### smrc peak fraction
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ann_smrc_counts_pca <- smrc_datasets$mssm_pca$ann_smrc_counts
ann_smrc_counts_pca$peak_read_fraction <- ann_smrc_counts_pca$MajorRNAReads/ann_smrc_counts_pca$Reads

pdf(file="iscience.figures/smrc_peak_mssm_pca.pdf", width = 6, height = 4)
ggplot( data= ann_smrc_counts_pca, aes((peak_read_fraction)) ) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity",color="black",fill="white")+
  stat_density(color="steelblue",fill=NA,size=2) + theme_bw() + 
  geom_rug(alpha=0.005,color="black") +
  geom_vline(aes(xintercept = (median(peak_read_fraction))), colour = "firebrick", linetype="dashed", size=2) +
  annotate(geom="text", x=0.5, y=1.75, label="Median", color="firebrick", size=8) +
  theme(text = element_text(size=rel(4.8)))  +
  labs(x="smRC Reads (Peak/Total)", y="Density")
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### smRCs complexity
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
pdf(file="iscience.figures/smrc_complexity_mssm_pca.pdf", width = 6, height = 4)
ggplot( data= ann_smrc_counts_pca, aes((Complexity)) ) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity",color="black",fill="white")+
  stat_density(color="steelblue",fill=NA,size=2) + theme_bw() + 
  geom_rug(alpha=0.005,color="black") +
  geom_vline(aes(xintercept = (median(Complexity))), colour = "firebrick", linetype="dashed", size=2) +
  annotate(geom="text", x=0.2, y=6, label="Median", color="firebrick", size=8) +
  theme(text = element_text(size=rel(4.8)))  +
  labs(x="smRC Complexity", y="Density")
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### smRCs unique read fraction
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### calculate this
ann_smrc_counts_pca$unique_fraction <- ann_smrc_counts_pca$UniqueReads/ann_smrc_counts_pca$Reads

pdf(file="iscience.figures/smrc_uniqueness_mssm_pca.pdf", width = 6, height = 4)
ggplot( data= ann_smrc_counts_pca, aes((unique_fraction)) ) + 
  geom_histogram(aes(y=..density..), alpha=0.5, position="identity",color="black",fill="white")+
  stat_density(color="steelblue",fill=NA,size=2) + theme_bw() + 
  geom_rug(alpha=0.005,color="black") +
  geom_vline(aes(xintercept = (median(unique_fraction))), colour = "firebrick", linetype="dashed", size=2) +
  annotate(geom="text", x=0.45, y=2, label="Median", color="firebrick", size=8) +
  theme(text = element_text(size=rel(4.8)))  +
  labs(x="smRC Reads (Unique/Total)", y="Density")
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### smRCs overlap with annotation
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### annotations used by Bojan
ann_g_gr <- readRDS(file="exo_RData/ann_g_gr.RDS")
head(ann_g_gr)

### split locus into coordinates
ann_smrc_counts_pca$seqnames <- gsub("chr","",tidyr::separate(data.frame(ann_smrc_counts_pca$Locus ), 1, sep=":", c("a","b"))$a)
ann_smrc_counts_pca$start <- tidyr::separate(data.frame( tidyr::separate(data.frame(ann_smrc_counts_pca$Locus ), 1, sep=":", c("a","b"))$b ), 1, sep="-", c("a","b"))$a
ann_smrc_counts_pca$end <- tidyr::separate(data.frame( tidyr::separate(data.frame(ann_smrc_counts_pca$Locus ), 1, sep=":", c("a","b"))$b ), 1, sep="-", c("a","b"))$b

### make ranges object
smrcs_gr <- makeGRangesFromDataFrame(ann_smrc_counts_pca[,c("seqnames","start","end","Locus")],keep.extra.columns = TRUE)

### always query, subject
hits <- findOverlaps(ann_g_gr,smrcs_gr) 
cc <- countOverlaps(ann_g_gr,smrcs_gr)
### extract data from overlaps
overlaps <- pintersect(ann_g_gr[queryHits(hits)], smrcs_gr[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(smrcs_gr[subjectHits(hits)])
percentOverlap2 <- width(overlaps) / width(ann_g_gr[queryHits(hits)])
### step1
ann_mapping_hits <- ann_g_gr[queryHits(hits)]
ann_mapping_hits<-as.data.frame(ann_mapping_hits,row.names = 1:length(overlaps))
smrc_hits <- smrcs_gr[subjectHits(hits)]
smrc_hits<-as.data.frame(smrc_hits,row.names = 1:length(overlaps))
### convert into dataframe
overlaps<-as.data.frame(overlaps,row.names = 1:length(overlaps))
percentOverlap<-as.data.frame(percentOverlap)
percentOverlap2<-as.data.frame(percentOverlap2)
### populate
smrc_hits$overlap_width<-overlaps$width
smrc_hits$overlap_with_smrc<-percentOverlap$percentOverlap
smrc_hits$overlap_with_annotation<-percentOverlap2$percentOverlap2
smrc_hits$overlap_width<-overlaps$width
smrc_hits$overlap_with_random<-percentOverlap$percentOverlap
smrc_hits$overlap_with_annotation<-percentOverlap2$percentOverlap2
smrc_hits$gene_biotype <- ann_mapping_hits$gene_biotype
smrc_hits$transcript_biotype <- ann_mapping_hits$transcript_biotype
smrc_hits$gene_name <- ann_mapping_hits$gene_name
smrc_hits$source <- ann_mapping_hits$source
smrc_hits$transcript_id <- ann_mapping_hits$transcript_id
### end
dim(ann_mapping_hits)
dim(smrc_hits)
###
head(smrc_hits)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### smRCs overlap with annotation FIGURE
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

test <- merge(ann_smrc_counts_pca, smrc_hits, by="Locus", all=TRUE)
ann_smrc_counts_pca_v2 <- test ; rm(test)

table(duplicated(smrc_hits$Locus))

### non hits to cero
ann_smrc_counts_pca_v2$overlap_with_annotation[is.na(ann_smrc_counts_pca_v2$overlap_with_annotation)] <- 0

### binary
ann_smrc_counts_pca_v2$binary_overlap_with_annotation <- ann_smrc_counts_pca_v2$overlap_with_annotation > 0
### length log10
ann_smrc_counts_pca_v2$log10.length <- log10(ann_smrc_counts_pca_v2$Length)

### true false binary for annotation
my_comparisons <- list( c( "TRUE","FALSE" ) )

pdf(file="iscience.figures/smrc_overlaps_with_annotation_mssm_pca.pdf", width = 9, height = 5)
ggarrange(

ggviolin(data=ann_smrc_counts_pca_v2, x="binary_overlap_with_annotation" , y="log10.length", color="binary_overlap_with_annotation",fill="binary_overlap_with_annotation",alpha=0.5) +
  theme_bw() +  #geom_quasirandom() + 
  geom_boxplot(width=0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
  theme(legend.position = "none") + ylim(1,5)+
  labs(x ='', y='Log10 (smRC Length)', title='',color="",fill="") + theme(text=element_text(size=21))

,

ggviolin(data=ann_smrc_counts_pca_v2, x="binary_overlap_with_annotation" , y="unique_fraction", color="binary_overlap_with_annotation",fill="binary_overlap_with_annotation",alpha=0.5) +
  theme_bw() +  #geom_quasirandom() + 
  geom_boxplot(width=0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
  theme(legend.position = "none") + ylim(0,1.25)+
  labs(x ='Annotation', y='smRC Reads\n(Unique/Total)', title='',color="",fill="") + theme(text=element_text(size=21))

,

ggviolin(data=ann_smrc_counts_pca_v2, x="binary_overlap_with_annotation" , y="Complexity", color="binary_overlap_with_annotation",fill="binary_overlap_with_annotation",alpha=0.5) +
  theme_bw() +  #geom_quasirandom() + 
  geom_boxplot(width=0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_fill_manual(values = c("firebrick","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
  theme(legend.position = "none") + ylim(0,1)+
  labs(x ='', y='smRC Complexity', title='',color="",fill="") + theme(text=element_text(size=21))

, ncol=3, nrow=1, align = 'hv')

dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


### Biotype for smRCs intersection
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
smrc_datasets$mssm_pca$smrc_counts

identical(rownames(smrc_datasets$mssm_pca$smrc_counts), (smrc_datasets$mssm_pca$ann_smrc_counts$Locus)) ### TRUE !

### counts matrix
d <- smrc_datasets$mssm_pca$smrc_counts
### length in pb
l <- smrc_datasets$mssm_pca$ann_smrc_counts$Length
### total mapped reads
cS <- colSums(d) #Total mapped reads per sample
### bam !!!!
rpkm <- (10^9)*t(t(d/l)/cS)
###
rpkm <- log10(rpkm+1)
rpkm <- melt(rpkm)
colnames(rpkm) <- c("Locus","full_sample","Log10RPKM")

rpkm <- merge(rpkm,ann_smrc_counts_pca_v2,by="Locus")

table(rpkm$gene_biotype)
#table(rpkm$transcript_biotype) # hmmm

rpkm$gene_biotype[is.na(rpkm$gene_biotype)] <- "uncharacterized"
rpkm$gene_biotype[rpkm$gene_biotype %in% "translated_processed_pseudogene" ] <- "translated_proc_pseudogene"

#library(tidytext)
#class(rpkm$Log10RPKM)

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

rpkm$clarity <- reorder(rpkm$gene_biotype, rpkm$Log10RPKM, FUN = mean)

pdf(file="iscience.figures/smrc_rpkm_per_biotype_mssm_pca.pdf", width = 12, height = 8)
ggplot(data=rpkm) +  
          aes(x=clarity,y=Log10RPKM) + theme_bw() +
          theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
          geom_boxplot(show.legend = FALSE,color="black",fill="grey") +  
          labs(x='Hg38 Gene Biotypes', y='Log10 RPKM') +
          theme(text=element_text(size=21)) +
          theme(plot.margin=unit(c(t=5,r=5,b=5,l=100),"pt"))
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### clean space
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
rm(list=setdiff(ls(), c("rpkm", "smrc_datasets", "voom_mssm_pca", "res_de_smRCs_cells_v_exos", "ann_smrc_counts_pca_v2","ann_smrc_counts_pca", 
                        "phen_prc", "prcomp_obj_mssm_pca", "smrc_deglist_obj","gg_color_hue","reorder_within","scale_x_reordered")))
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### smRCs complexity by peakness cells vs exosomes
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dim(res_de_smRCs_cells_v_exos)

### add variable with names
res_de_smRCs_cells_v_exos$Locus <- rownames(res_de_smRCs_cells_v_exos)
### sort
ann_smrc_counts_pca <- ann_smrc_counts_pca[ match(res_de_smRCs_cells_v_exos$Locus, ann_smrc_counts_pca$Locus) , ] 
### compare
identical(ann_smrc_counts_pca$Locus, res_de_smRCs_cells_v_exos$Locus)

ann_smrc_counts_pca$logFC <- res_de_smRCs_cells_v_exos$logFC
ann_smrc_counts_pca$adj.P.Val <- res_de_smRCs_cells_v_exos$adj.P.Val

ix <- which(ann_smrc_counts_pca$adj.P.Val<0.01)

pdf(file="iscience.figures/smrc_complex_peak_cell_exosome_mssm_pca.pdf", width = 7, height = 4)
ggplot(ann_smrc_counts_pca[ix,], aes(x=peak_read_fraction, y=Complexity, color=logFC)) +  theme_bw() +
  geom_point(alpha=0.5,size=0.4) +
  #geom_hex(aes(color=logFC),alpha=0.1)+
  geom_rug(alpha=0.2) +
  geom_density2d(color="black",alpha=0.5) +
  #scale_colour_gradient2(midpoint=mean(ann_smrc_counts_pca$logFC),low = "red", mid = "white", high = "blue" ,space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour") +
  scale_colour_gradient2(low = "red", mid = "white", high = "blue" ,space = "Lab", na.value = "grey50", guide = "colourbar",aesthetics = "colour") +
  theme(legend.position = c(0.9, 0.7), legend.title = element_text(size = 8),legend.key = element_rect(colour = "transparent", fill = "transparent")  ) + #, legend.justification=c(1,0)
  labs(x='smRC Reads (Peak/Total)',y='smRC Complexity') +
  annotate(geom="text", x=0.25, y=0.6, label='atop(bold("Cell"))',color="firebrick",size=10, parse = T) +
  annotate(geom="text", x=0.87, y=0.19, label='atop(bold("Exosome"))',color="steelblue",size=10, parse = T) +
  theme(text=element_text(size=21)) +
  theme(legend.title = element_text(size=rel(1)))
dev.off()
### Complexity ( N-Distinct Reads ) / ( N-SMRC Reads )
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Cell exosome correlation
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
phen_prc <- phen_prc[ match( colnames(voom_mssm_pca$E), phen_prc$full_sample), ]
identical(as.character(phen_prc$full_sample),colnames(voom_mssm_pca$E))

temporal_variable_for_figure <- data.frame( Exosome = rowMeans(voom_mssm_pca$E[,which(phen_prc$Cell_Type %in% "Exosome")]) ,
                                            Cell = rowMeans(voom_mssm_pca$E[,which(phen_prc$Cell_Type %in% "Cell")]) )

cor.test(temporal_variable_for_figure$Exosome,temporal_variable_for_figure$Cell)

pdf(file="iscience.figures/smrc_correlation_cell_exosome_mssm_pca.pdf", width = 6, height = 5)
ggplot(temporal_variable_for_figure, aes(Exosome, Cell)) + geom_point(size=0.1) + stat_density2d() + 
  geom_rug(alpha=0.01,color="black")  +
  geom_abline(intercept = 0, slope = 1, color='firebrick', linetype = 2,size=2) + 
  stat_smooth(method = 'lm', color = 'purple') + 
  theme_bw() +
  labs(x="Log2CPM smRCs (Exosome)", y="Log2CPM smRCs (Cell)", title="") +
  annotate(geom="text", x=8, y=17,    label="Spearman Rho: -0.07\np.value: 2.2e-16",color="firebrick",size=7) +
  theme(text=element_text(size=21))
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Differential expression urine serum normal and tumor
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### check names
phen_prc <- phen_prc[ match( colnames(voom_mssm_pca$E),  as.character(phen_prc$full_sample) ) , ]
table(as.character(phen_prc$full_sample) %in% colnames(voom_mssm_pca$E))
identical(as.character(phen_prc$full_sample),colnames(voom_mssm_pca$E))

### design
design = model.matrix( ~ 0 + Tissue, data = phen_prc )
### names
colnames(design) <- gsub("Tissue","",colnames(design))
### contr.matrix
contr.matrix <- makeContrasts( "Normal vs Tumor" = Normal - Tumor,  
                               "Normal vs Serum" = Normal - Serum,  
                               "Normal vs Urine" = Normal - Urine,  
                               "Tumor vs Serum" = Tumor - Serum,
                               "Tumor vs Urine" = Tumor - Urine,
                               "Serum vs Urine" = Serum - Urine,
                               "Serum" = Serum - (Urine+Tumor+Normal)/3,
                               "Urine" = Urine - (Serum+Tumor+Normal)/3,
                               "Normal" = Normal - (Urine+Tumor+Serum)/3,
                               "Tumor" = Tumor - (Urine+Serum+Normal)/3,
                               levels = colnames(design))

### statistics
block2 <- as.numeric(as.factor( phen_prc$Patient_ID ))
### Correlation factors
dupcor2 <- duplicateCorrelation(voom_mssm_pca, design, block=block2)
### linear fit with correlation factors
vfit <- lmFit(voom_mssm_pca, design, block=block2, correlation=dupcor2$consensus)
### vfit <- lmFit(voom_mssm_pca, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
### bayesian statistics
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

### chain results
lmfreq_results <- list()
for ( i in 1:ncol(summary(decideTests(efit))) ) { lmfreq_results[[i]] <- topTable(efit, coef=i,n=Inf,adjust.method="fdr") 
lmfreq_results[[i]]$Comparison <- colnames(summary(decideTests(efit)))[i]
lmfreq_results[[i]]$Locus <- rownames( topTable(efit, coef=i,n=Inf,adjust.method="fdr") ) }
### get it together
res_de_smRCs_tissues <- do.call(rbind,lmfreq_results)
rownames(res_de_smRCs_tissues) <- NULL
res_de_smRCs_tissues$nLogFDR <- -log10(res_de_smRCs_tissues$adj.P.Val)

### cleanse your soul
rm(vfit,efit,contr.matrix,design,Scree.plot,mydata,pca1,pca2,ix,iy,form,dupcor2,block2,lmfreq_results,temporal_variable_for_figure)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### pHeatmap plot top markers
#### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### annotation
annotation_col <- data.frame( RNA.Origin = phen_prc$Cell_Type, Tissue = phen_prc$Tissue, Isolation = phen_prc$technology); rownames(annotation_col) <-  phen_prc$full_sample
annotation_col$Isolation <- as.character(annotation_col$Isolation)

annotation_col$Isolation[annotation_col$Isolation %in% "bulk"] <- "none"

gg_color_hue(4)
hcl.colors(6,palette="Spectral")

annotation_colors <- list ( RNA.Origin = c( Exosome="steelblue", Cell="firebrick"), 
                            High = c(Normal="#F8766D",Tumor="#7CAE00",Serum="#00BFC4", Urine="#C77CFF"), 
                            Tissue = c(Normal="#F8766D",Tumor="#7CAE00",Serum="#00BFC4", Urine="#C77CFF"),
                            Isolation = c(none="white",DLD="black",UC="grey"),
                            Biotype = c(lincRNA="#A71B4B",`Other ncRNAs`="#ED820A",`Promoter lncRNA`="#FCDE85",`Protein Coding`="#BAEEAE",snoRNA="#00B1B5",Uncharacterized="#584B9F") )

### annotation row 
annotation_row <- res_de_smRCs_tissues
annotation_row$High <- "ns"
annotation_row$High[which(annotation_row$adj.P.Val<0.05 & annotation_row$Comparison %in% "Normal" & annotation_row$logFC > 0)][1:25] <- "Normal"
annotation_row$High[which(annotation_row$adj.P.Val<0.05 & annotation_row$Comparison %in% "Tumor" & annotation_row$logFC > 0)][1:25] <- "Tumor"
annotation_row$High[which(annotation_row$adj.P.Val<0.05 & annotation_row$Comparison %in% "Serum" & annotation_row$logFC > 0)][1:25] <- "Serum"
annotation_row$High[which(annotation_row$adj.P.Val<0.05 & annotation_row$Comparison %in% "Urine" & annotation_row$logFC > 0)][1:25] <- "Urine"

iy <- which(annotation_row$High %in% c("Normal","Tumor","Serum","Urine"))

annotation_row <- annotation_row[iy,]
annotation_row <- annotation_row[(!duplicated(annotation_row$Locus)),]
rownames(annotation_row) <- annotation_row$Locus

annotation_row <- merge(annotation_row, ann_smrc_counts_pca_v2, by="Locus")
annotation_row <- annotation_row[,c("High","AveExpr","gene_biotype","Locus")]
annotation_row <- unique(annotation_row)
annotation_row <- annotation_row[(!duplicated(annotation_row$Locus)),]
annotation_row$gene_biotype[is.na(annotation_row$gene_biotype)] <- "Uncharacterized"
colnames(annotation_row)[3] <- "Biotype"
rownames(annotation_row) <- annotation_row$Locus
annotation_row <- annotation_row[,-4]

annotation_row$Biotype[ annotation_row$Biotype %in% c("antisense","processed_transcript")] <- "Other ncRNAs"   
annotation_row$Biotype[ annotation_row$Biotype %in% c("protein_coding")] <- "Protein Coding"   
annotation_row$Biotype[ annotation_row$Biotype %in% c("bidirectional_promoter_lncRNA")] <- "Promoter lncRNA"   

ix <- which(rownames(voom_mssm_pca$E) %in% rownames(annotation_row) )

pdf(file="iscience.figures/heatmap_tissues_mssm_pca.pdf",width = 20,height = 8)
pheatmap(t(voom_mssm_pca$E[ix,]), 
         color = colorRampPalette(c("steelblue", "white","firebrick"))(255),  
         cluster_cols = T, 
         cluster_rows=T, 
         angle_col = 45,
         main ='',
         name='z-score',
         annotation_row = annotation_col,
         annotation_col = annotation_row,
         annotation_colors = annotation_colors,
         border_color="white", 
         show_rownames = FALSE, 
         show_colnames = TRUE, 
         clustering_distance_rows ="euclidean", scale="row")
dev.off()
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### smRCS common to different datasets
### venn & exosome cell/uniqueness
### Note: couldnt reproduce Bojan's figure - using previous version.
###
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#data_for_venn <- list( CRC = rownames(smrc_datasets$CRC$smrc_counts),
#      PCA = rownames(smrc_datasets$mssm_pca$smrc_counts),
#      Plasma = rownames(smrc_datasets$Plasma$smrc_counts) )
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#library(Vennerable)
#Venn(data_for_venn)



### Clustering for CRC
#### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### Voom + dgelist CRC
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dim(smrc_datasets$CRC$smrc_counts)
###
smrc_deglist_CRC_obj <- DGEList(counts=smrc_datasets$CRC$smrc_counts) 
smrc_deglist_CRC_obj <- calcNormFactors(smrc_deglist_CRC_obj, method = "TMM") 
smrc_deglist_CRC_obj = estimateCommonDisp(smrc_deglist_CRC_obj, verbose=TRUE)
### voom
voom_CRC <- voom(smrc_deglist_CRC_obj, plot = TRUE, save.plot = TRUE) 
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Principal component analysis CRC
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PCA
mydata <- prcomp( t(voom_CRC$E) , scale=TRUE, center = TRUE )

### Scree
Scree.plot <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw()
eigen <- get_eigenvalue(mydata)
pca1 <- format(round(eigen$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(eigen$variance.percent[2], 2), nsmall = 2)
Scree.plot ; rm(eigen)

phen_CRC <- data.frame( Type = colnames(voom_CRC$E)  ,  RNA.Origin = colnames(voom_CRC$E) )

phen_CRC$RNA.Origin[grep("exo",phen_CRC$RNA.Origin)] <- "Exosome"
phen_CRC$RNA.Origin[grep("cell",phen_CRC$RNA.Origin)] <- "Cell"

phen_CRC$Type[grep("DLD",phen_CRC$Type)] <- "Wild"
phen_CRC$Type[grep("DKS",phen_CRC$Type)] <- "Wild"
phen_CRC$Type[grep("DKO",phen_CRC$Type)] <- "Mutant"
phen_CRC$full_sample <- colnames(voom_CRC$E)
phen_CRC$Condition <- paste(phen_CRC$Type,phen_CRC$RNA.Origin,sep=" ")

prcomp_obj_CRC <- as.data.frame(mydata$x)
prcomp_obj_CRC$full_sample <- rownames(prcomp_obj_CRC)
prcomp_obj_CRC$lib.size <- smrc_deglist_CRC_obj$samples$lib.size
prcomp_obj_CRC$norm.factors <- smrc_deglist_CRC_obj$samples$norm.factors
prcomp_obj_CRC <- merge(prcomp_obj_CRC, phen_CRC, by="full_sample")


pdf(file="iscience.figures/prcomp_CRC.pdf",width = 8,height = 6)
ggplot(data=prcomp_obj_CRC, aes(x=PC1, y=PC2, color = Condition)) + 
  geom_point() + 
  geom_label_repel(aes(label=full_sample),force = 10, max.iter = 5000,show.legend = FALSE) +
  theme_bw() + 
  labs(x=paste("PC1:",pca1,"%",sep=""), y=paste("PC1:", pca2,"%", sep=""), color="Type / Origin") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(4.0)),
        plot.title=element_text(size=rel(4.8)) ) + 
  theme(legend.position="bottom") + 
  guides(color=guide_legend(ncol=2,override.aes = list(size = 5))) 
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Voom + dgelist plasma
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dim(smrc_datasets$Plasma$smrc_counts)
###
smrc_deglist_plasma_obj <- DGEList(counts=smrc_datasets$Plasma$smrc_counts) 
smrc_deglist_plasma_obj <- calcNormFactors(smrc_deglist_plasma_obj, method = "TMM") 
smrc_deglist_plasma_obj = estimateCommonDisp(smrc_deglist_plasma_obj, verbose=TRUE)
### voom
voom_plasma <- voom(smrc_deglist_plasma_obj, plot = TRUE, save.plot = TRUE) 
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Principal component analysis CRC
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### PCA
mydata <- prcomp( t(voom_plasma$E) , scale=TRUE, center = TRUE )

### Scree
Scree.plot <- fviz_eig(mydata, addlabels = T, barfill = "lightsteelblue3", barcolor = "lightsteelblue3") + theme_bw()
eigen <- get_eigenvalue(mydata)
pca1 <- format(round(eigen$variance.percent[1], 2), nsmall = 2)
pca2 <- format(round(eigen$variance.percent[2], 2), nsmall = 2)
Scree.plot ; rm(eigen)


prcomp_obj_plasma <- as.data.frame(mydata$x)
prcomp_obj_plasma$full_sample <- rownames(prcomp_obj_plasma)
prcomp_obj_plasma$lib.size <- smrc_deglist_plasma_obj$samples$lib.size
prcomp_obj_plasma$norm.factors <- smrc_deglist_plasma_obj$samples$norm.factors
#prcomp_obj_CRC <- merge(prcomp_obj_CRC, phen_plasma, by="full_sample")

pdf(file="iscience.figures/prcomp_plasma.pdf",width = 8,height = 6)
ggplot(data=prcomp_obj_plasma, aes(x=PC1, y=PC2)) + 
  geom_point() + 
  #geom_label_repel(aes(label=full_sample),force = 10, max.iter = 5000,show.legend = FALSE) +
  theme_bw() + 
  labs(x=paste("PC1:",pca1,"%",sep=""), y=paste("PC1:", pca2,"%", sep=""), color="Type / Origin") +
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(4.0)),
        plot.title=element_text(size=rel(4.8)) ) + 
  theme(legend.position="bottom") + 
  guides(color=guide_legend(ncol=2,override.aes = list(size = 5))) 
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### smRCs overlap with annotation CRC
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### annotations used by Bojan
ann_g_gr <- readRDS(file="exo_RData/ann_g_gr.RDS")
head(ann_g_gr)

ann_smrc_counts_crc <- smrc_datasets$CRC$ann_smrc_counts

### split locus into coordinates
ann_smrc_counts_crc$seqnames <- gsub("chr","",tidyr::separate(data.frame(ann_smrc_counts_crc$Locus ), 1, sep=":", c("a","b"))$a)
ann_smrc_counts_crc$start <- tidyr::separate(data.frame( tidyr::separate(data.frame(ann_smrc_counts_crc$Locus ), 1, sep=":", c("a","b"))$b ), 1, sep="-", c("a","b"))$a
ann_smrc_counts_crc$end <- tidyr::separate(data.frame( tidyr::separate(data.frame(ann_smrc_counts_crc$Locus ), 1, sep=":", c("a","b"))$b ), 1, sep="-", c("a","b"))$b

### make ranges object
smrcs_gr <- makeGRangesFromDataFrame(ann_smrc_counts_crc[,c("seqnames","start","end","Locus")],keep.extra.columns = TRUE)

### always query, subject
hits <- findOverlaps(ann_g_gr,smrcs_gr) 
cc <- countOverlaps(ann_g_gr,smrcs_gr)
### extract data from overlaps
overlaps <- pintersect(ann_g_gr[queryHits(hits)], smrcs_gr[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(smrcs_gr[subjectHits(hits)])
percentOverlap2 <- width(overlaps) / width(ann_g_gr[queryHits(hits)])
### step1
ann_mapping_hits <- ann_g_gr[queryHits(hits)]
ann_mapping_hits<-as.data.frame(ann_mapping_hits,row.names = 1:length(overlaps))
smrc_hits <- smrcs_gr[subjectHits(hits)]
smrc_hits<-as.data.frame(smrc_hits,row.names = 1:length(overlaps))
### convert into dataframe
overlaps<-as.data.frame(overlaps,row.names = 1:length(overlaps))
percentOverlap<-as.data.frame(percentOverlap)
percentOverlap2<-as.data.frame(percentOverlap2)
### populate
smrc_hits$overlap_width<-overlaps$width
smrc_hits$overlap_with_smrc<-percentOverlap$percentOverlap
smrc_hits$overlap_with_annotation<-percentOverlap2$percentOverlap2
smrc_hits$overlap_width<-overlaps$width
smrc_hits$overlap_with_random<-percentOverlap$percentOverlap
smrc_hits$overlap_with_annotation<-percentOverlap2$percentOverlap2
smrc_hits$gene_biotype <- ann_mapping_hits$gene_biotype
smrc_hits$transcript_biotype <- ann_mapping_hits$transcript_biotype
smrc_hits$gene_name <- ann_mapping_hits$gene_name
smrc_hits$source <- ann_mapping_hits$source
smrc_hits$transcript_id <- ann_mapping_hits$transcript_id
### end
dim(ann_mapping_hits)
dim(smrc_hits)
###
head(smrc_hits)

test <- merge(ann_smrc_counts_crc, smrc_hits, by="Locus", all=TRUE)
ann_smrc_counts_crc_v2 <- test ; rm(test)

table(duplicated(smrc_hits$Locus))

### non hits to cero
ann_smrc_counts_crc_v2$overlap_with_annotation[is.na(ann_smrc_counts_crc_v2$overlap_with_annotation)] <- 0

### binary
ann_smrc_counts_crc_v2$binary_overlap_with_annotation <- ann_smrc_counts_crc_v2$overlap_with_annotation > 0
### length log10
ann_smrc_counts_crc_v2$log10.length <- log10(ann_smrc_counts_crc_v2$Length)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### smRCs overlap with annotation PLASMA
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### annotations used by Bojan
ann_g_gr <- readRDS(file="exo_RData/ann_g_gr.RDS")
head(ann_g_gr)

ann_smrc_counts_plasma <- smrc_datasets$Plasma$ann_smrc_counts

### split locus into coordinates
ann_smrc_counts_plasma$seqnames <- gsub("chr","",tidyr::separate(data.frame(ann_smrc_counts_plasma$Locus ), 1, sep=":", c("a","b"))$a)
ann_smrc_counts_plasma$start <- tidyr::separate(data.frame( tidyr::separate(data.frame(ann_smrc_counts_plasma$Locus ), 1, sep=":", c("a","b"))$b ), 1, sep="-", c("a","b"))$a
ann_smrc_counts_plasma$end <- tidyr::separate(data.frame( tidyr::separate(data.frame(ann_smrc_counts_plasma$Locus ), 1, sep=":", c("a","b"))$b ), 1, sep="-", c("a","b"))$b

### make ranges object
smrcs_gr <- makeGRangesFromDataFrame(ann_smrc_counts_plasma[,c("seqnames","start","end","Locus")],keep.extra.columns = TRUE)

### always query, subject
hits <- findOverlaps(ann_g_gr,smrcs_gr) 
cc <- countOverlaps(ann_g_gr,smrcs_gr)
### extract data from overlaps
overlaps <- pintersect(ann_g_gr[queryHits(hits)], smrcs_gr[subjectHits(hits)])
percentOverlap <- width(overlaps) / width(smrcs_gr[subjectHits(hits)])
percentOverlap2 <- width(overlaps) / width(ann_g_gr[queryHits(hits)])
### step1
ann_mapping_hits <- ann_g_gr[queryHits(hits)]
ann_mapping_hits<-as.data.frame(ann_mapping_hits,row.names = 1:length(overlaps))
smrc_hits <- smrcs_gr[subjectHits(hits)]
smrc_hits<-as.data.frame(smrc_hits,row.names = 1:length(overlaps))
### convert into dataframe
overlaps<-as.data.frame(overlaps,row.names = 1:length(overlaps))
percentOverlap<-as.data.frame(percentOverlap)
percentOverlap2<-as.data.frame(percentOverlap2)
### populate
smrc_hits$overlap_width<-overlaps$width
smrc_hits$overlap_with_smrc<-percentOverlap$percentOverlap
smrc_hits$overlap_with_annotation<-percentOverlap2$percentOverlap2
smrc_hits$overlap_width<-overlaps$width
smrc_hits$overlap_with_random<-percentOverlap$percentOverlap
smrc_hits$overlap_with_annotation<-percentOverlap2$percentOverlap2
smrc_hits$gene_biotype <- ann_mapping_hits$gene_biotype
smrc_hits$transcript_biotype <- ann_mapping_hits$transcript_biotype
smrc_hits$gene_name <- ann_mapping_hits$gene_name
smrc_hits$source <- ann_mapping_hits$source
smrc_hits$transcript_id <- ann_mapping_hits$transcript_id
### end
dim(ann_mapping_hits)
dim(smrc_hits)
###
head(smrc_hits)

test <- merge(ann_smrc_counts_plasma, smrc_hits, by="Locus", all=TRUE)
ann_smrc_counts_plasma_v2 <- test ; rm(test)

table(duplicated(smrc_hits$Locus))

### non hits to cero
ann_smrc_counts_plasma_v2$overlap_with_annotation[is.na(ann_smrc_counts_plasma_v2$overlap_with_annotation)] <- 0

### binary
ann_smrc_counts_plasma_v2$binary_overlap_with_annotation <- ann_smrc_counts_plasma_v2$overlap_with_annotation > 0
### length log10
ann_smrc_counts_plasma_v2$log10.length <- log10(ann_smrc_counts_plasma_v2$Length)

### clean
rm(hits,mydata,overlaps,percentOverlap,percentOverlap2,data_for_venn,Scree.plot,smrc_hits,smrcs_gr,ann_g_gr,pca1,pca2,PCA,ix,iy,i,CRC,cc,ann_mapping_hits)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Differential expression Cellular to exosomal
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#phen_CRC <- phen_CRC[ match( colnames(voom_CRC$E),  as.character(phen_CRC$full_sample) ) , ]
#table(as.character(phen_CRC$full_sample) %in% colnames(voom_CRC$E))
identical(as.character(phen_CRC$full_sample),colnames(voom_CRC$E))

### design
design = model.matrix( ~ 0 + RNA.Origin, data = phen_CRC )
colnames(design) <- gsub("RNA.Origin","",colnames(design))

### contr.matrix
contr.matrix <- makeContrasts( "Exosomes vs Cells" = Exosome - Cell,  levels = colnames(design))

### statistics
vfit <- lmFit(voom_CRC, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit, robust = TRUE) ; summary(decideTests(efit))

### results into table
res_de_smRCs_CRC_cells_v_exos <- topTable(efit, coef="Exosomes vs Cells", n=Inf, adjust.method="BH") 
### add FDR converted
res_de_smRCs_CRC_cells_v_exos$nLogFDR <- -log10(res_de_smRCs_CRC_cells_v_exos$adj.P.Val)

rm(vfit,efit,contr.matrix,design,Scree.plot,mydata,pca1,pca2,ix,iy,form)
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Volcano plot CElls exosomes CRC
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

### label crc results
res_de_smRCs_CRC_cells_v_exos$Locus <- rownames(res_de_smRCs_CRC_cells_v_exos)

volcano.plot<- ggplot(data=res_de_smRCs_CRC_cells_v_exos,aes(x=logFC,y=nLogFDR)) + 
  geom_point(data=res_de_smRCs_CRC_cells_v_exos[res_de_smRCs_CRC_cells_v_exos$adj.P.Val>0.05,],color="slategray",alpha=0.7, size = 1) + 
  geom_point(data=res_de_smRCs_CRC_cells_v_exos[res_de_smRCs_CRC_cells_v_exos$logFC > 0 & res_de_smRCs_CRC_cells_v_exos$adj.P.Val<0.05,],color="steelblue",size = 1,alpha=0.6) +      
  geom_point(data=res_de_smRCs_CRC_cells_v_exos[res_de_smRCs_CRC_cells_v_exos$logFC < 0 & res_de_smRCs_CRC_cells_v_exos$adj.P.Val<0.05,],color="firebrick",size = 1,alpha=0.6) + 
  
  geom_point(data= res_de_smRCs_CRC_cells_v_exos[res_de_smRCs_CRC_cells_v_exos$logFC > 0 & res_de_smRCs_CRC_cells_v_exos$adj.P.Val<0.05,], color = "steelblue",size = 1) +
  geom_point(data= res_de_smRCs_CRC_cells_v_exos[res_de_smRCs_CRC_cells_v_exos$logFC < 0 & res_de_smRCs_CRC_cells_v_exos$adj.P.Val<0.05,], color = "firebrick",size = 1) +
  
  geom_vline(xintercept = 0, linetype = 2, color = 'black',size=1.5) +
  geom_hline(yintercept = 1.3, linetype = 2, color = 'black',size=1.5) +
  
  xlab("Log2 Fold Change") + ylab("-Log10 FDR") +
  annotate(geom="text", x=-4.5, y=16.5, label='atop(bold("Cell"))',color="firebrick",size=10, parse = T) +
  annotate(geom="text", x=6.5, y=16.5, label='atop(bold("Exosome"))',color="steelblue",size=10, parse = T) +
  annotate(geom="text", x=-8.5, y=1.5, label='atop(bold("Threshold FDR < 5%"))',color="black",size=6, parse = T) +
  theme_bw() + 
  theme(text = element_text(size=rel(4.8))) 

pdf(file="iscience.figures/volcano_CRC.pdf",width = 8,height = 5)
volcano.plot
dev.off()

rm(volcano.plot)


### Peak length vs logFC
#### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
res_de_smRCs_cells_v_exos$Locus[res_de_smRCs_cells_v_exos$adj.P.Val<0.01 & res_de_smRCs_cells_v_exos$logFC>1]

test <- ann_smrc_counts_pca_v2[,which( colnames(ann_smrc_counts_pca_v2) %in% c("MajorRNAReads","Locus") )]
test <- merge(res_de_smRCs_cells_v_exos,test,by="Locus")
ix <- which(test$adj.P.Val<0.001)
ggplot(test[ix,])+ aes(x=log10(MajorRNAReads),y=logFC) + geom_point() + geom_density_2d()

rm(test)
#### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ann_smrc_counts_plasma_v2$peak_read_fraction <- ann_smrc_counts_plasma_v2$MajorRNAReads/ann_smrc_counts_plasma_v2$Reads
ann_smrc_counts_crc_v2$peak_read_fraction <- ann_smrc_counts_crc_v2$MajorRNAReads/ann_smrc_counts_crc_v2$Reads

ann_smrc_counts_plasma_v2$unique_fraction <- ann_smrc_counts_plasma_v2$UniqueReads/ann_smrc_counts_plasma_v2$Reads
ann_smrc_counts_crc_v2$unique_fraction <- ann_smrc_counts_crc_v2$UniqueReads/ann_smrc_counts_crc_v2$Reads

dim(ann_smrc_counts_pca_v2)
dim(ann_smrc_counts_crc_v2)
dim(ann_smrc_counts_plasma_v2)

pca <- ann_smrc_counts_pca_v2[,c("Complexity","Length","peak_read_fraction","unique_fraction","binary_overlap_with_annotation","MajorRNAReads","gene_biotype","Locus")]
crc <- ann_smrc_counts_crc_v2[,c("Complexity","Length","peak_read_fraction","unique_fraction","binary_overlap_with_annotation","MajorRNAReads","gene_biotype","Locus")]
plasma <- ann_smrc_counts_plasma_v2[,c("Complexity","Length","peak_read_fraction","unique_fraction","binary_overlap_with_annotation","MajorRNAReads","gene_biotype","Locus")]

plasma$study <- "Plasma"
crc$study <- "CRC"
pca$study <- "PCA"

plasma$classification <- "Exosome"

res_de_smRCs_cells_v_exos$classification <- "ns"
res_de_smRCs_cells_v_exos$classification[ res_de_smRCs_cells_v_exos$adj.P.Val<0.05 & res_de_smRCs_cells_v_exos$logFC > 3] <- "Exosome"
res_de_smRCs_cells_v_exos$classification[ res_de_smRCs_cells_v_exos$adj.P.Val<0.05 & res_de_smRCs_cells_v_exos$logFC < -2.7] <- "Cell"

res_de_smRCs_CRC_cells_v_exos$classification <- 'ns'
res_de_smRCs_CRC_cells_v_exos$classification[ res_de_smRCs_CRC_cells_v_exos$adj.P.Val<0.05 & res_de_smRCs_CRC_cells_v_exos$logFC > 4.6] <- "Exosome"
res_de_smRCs_CRC_cells_v_exos$classification[ res_de_smRCs_CRC_cells_v_exos$adj.P.Val<0.05 & res_de_smRCs_CRC_cells_v_exos$logFC < -3.5] <- "Cell"

test <- res_de_smRCs_CRC_cells_v_exos[, c("Locus","classification")]
crc <- merge(crc,test,by="Locus")

test <- res_de_smRCs_cells_v_exos[, c("Locus","classification")]
pca <- merge(pca,test,by="Locus")

pca <- pca[ , order(colnames(pca)) ]
crc <- crc[ , order(colnames(crc)) ]
plasma <- plasma[ , order(colnames(plasma)) ]

dim(pca)
dim(crc)
dim(plasma)

test <- rbind(plasma,crc,pca)
test$condition <- paste(test$study,test$classification)
test$log10length<- log10(test$Length+1)

test$condition <- factor(test$condition, levels=c("CRC Cell","PCA Cell","CRC Exosome","Plasma Exosome","PCA Exosome","PCA ns","CRC ns"))

ix <- which(!test$condition %in% c("CRC ns","PCA ns"))
#ggplot(data=test[ix,]) + aes(x=condition,y=Complexity) + geom_boxplot()

my_comparisons <- list( c("CRC Cell", "CRC Exosome"),
                        c("CRC Cell","Plasma Exosome"),
                        c("PCA Cell","Plasma Exosome"),
                        c("CRC Exosome","Plasma Exosome"))

pdf(file="iscience.figures/smRC_properties_comparison_cross_dataset.pdf",width = 10,height = 12)
ggarrange(

ggboxplot(data=test[ix,], x="condition" , y="log10length",color='classification',size=1) +
  theme_bw() +  #geom_quasirandom() + 
  #geom_boxplot(width=0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  #scale_fill_manual(values = c("firebrick","steelblue")) +
  #scale_color_manual(values = c("firebrick","firebrick","steelblue","steelblue","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
  theme(legend.position = "none") + #ylim(1,5)+
  labs(x ='', y='Log10 (smRC Length)', title='',color="",fill="") + theme(text=element_text(size=21))
,
ggboxplot(data=test[ix,], x="condition" , y="Complexity", color='classification',size=1) +
  theme_bw() +  #geom_quasirandom() + 
  #geom_boxplot(width=0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  #scale_fill_manual(values = c("firebrick","steelblue")) +
  #scale_color_manual(values = c("firebrick","firebrick","steelblue","steelblue","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
  theme(legend.position = "none") + #ylim(1,5)+
  labs(x ='', y='smRC Complexity', title='',color="",fill="") + theme(text=element_text(size=21))
,
ggboxplot(data=test[ix,], x="condition" , y="unique_fraction", color='classification',size=1) +
  theme_bw() +  #geom_quasirandom() + 
  #geom_boxplot(width=0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  #scale_fill_manual(values = c("firebrick","steelblue")) +
  #scale_color_manual(values = c("firebrick","firebrick","steelblue","steelblue","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE,
                     oalette = c("firebrick","firebrick","steelblue","steelblue","steelblue")) +
  theme(legend.position = "none") + #ylim(1,5)+
  labs(x ='', y='smRC Reads (Unique/Total)', title='',color="",fill="") + theme(text=element_text(size=21))
,
ggboxplot(data=test[ix,], x="condition" , y="peak_read_fraction", color='classification',size=1) +
  theme_bw() +  #geom_quasirandom() + 
  #geom_boxplot(width=0.1) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  #scale_fill_manual(values = c("firebrick","steelblue")) +
  #scale_color_manual(values = c("firebrick","firebrick","steelblue","steelblue","steelblue")) +
  scale_color_manual(values = c("firebrick","steelblue")) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", paired=FALSE, label="p.signif", hide.ns = FALSE ) +
  theme(legend.position = "none") + #ylim(1,5)+
  labs(x ='', y='smRC Reads (Peak/Total)', title='',color="",fill="") + theme(text=element_text(size=21)) 
, align = 'hv', ncol=2,nrow=2
)
dev.off()
#### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Venn figures, biotype composition
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
test$gene_biotype[is.na(test$gene_biotype)] <- "Uncharacterized"

testv2 <- as.matrix(table(test$gene_biotype,test$condition))
testv2 <- rbind(testv2,colSums(testv2))
rownames(testv2)[40] <- "Total"

testv2 <- testv2[,-c(6:7)]

testv2[,1] <- (testv2[1:40,1]/testv2[40,1]) *100
testv2[,2] <- (testv2[1:40,2]/testv2[40,2]) *100
testv2[,3] <- (testv2[1:40,3]/testv2[40,3]) *100
testv2[,4] <- (testv2[1:40,4]/testv2[40,4]) *100
testv2[,5] <- (testv2[1:40,5]/testv2[40,5]) *100

testv2 <- testv2[-40,]
testv2 <- as.data.frame(testv2)
testv2$biotype <- rownames(testv2)
testv2 <- melt(testv2,id.vars = 'biotype')

testv2$biotype2 <- testv2$biotype
testv2$biotype2[grep("pseudogene",testv2$biotype2)] <- "Pseudogenes"

testv2$biotype2[which(testv2$biotype2 %in% c("IG_C_gene","IG_V_gene","TR_C_gene","TR_J_gene","TR_V_gene"))] <- "VDJ"
testv2$biotype2[which(testv2$biotype2 %in% c("antisense","sense_overlapping","TEC","sense_intronic","3prime_overlapping_ncRNA","processed_transcript","non_coding","misc_RNA"))] <- "Other ncRNAs"

testv2$biotype2[which(testv2$biotype2 %in% c("bidirectional_promoter_lncRNA","lincRNA","macro_lncRNA"))] <- "lncRNAs"
testv2$biotype2[which(testv2$biotype2 %in% c("protein_coding"))] <- "Protein Coding"

#ggplot(testv2) + aes(x=Var2,y=Freq) + geom_bar(stat='identity')

pdf(file="iscience.figures/smRCs_biotype_across_datasets.pdf",width = 16,height = 4)
ggplot(testv2, aes(fill=biotype2, y=value, x=variable)) +  
  geom_bar(position="stack", stat="identity") +
  theme_bw() + coord_flip() + theme(text = element_text(size=rel(4.8)))  +
  theme(legend.position="bottom") + 
  theme(text = element_text(size=rel(4.8)),
        legend.text=element_text(size=rel(4.8)),
        plot.title=element_text(size=rel(4.8)),
        axis.text.x=element_text(size=rel(5.5)),
        axis.text.y=element_text(size=rel(5.5))
        ) + 
  theme(legend.position="bottom") +
  labs(x="smRCs (%)",y="",fill="Biotype") +
  guides(color=guide_legend(ncol=3,override.aes = list(size = 5))) 
dev.off()
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



### Alignment Summary stats
### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

mapping_stats <- data.frame(
  Stat = c("Unmapped","Unique", "Multimapped", "QC-Filtered"),
  Value = c(40, 17.5, 28.9, 13.6))

sum(mapping_stats$Value)

mapping_stats <- mapping_stats %>% arrange(desc(Stat)) %>% mutate(lab.ypos = cumsum(Value) - 0.5*Value)

pdf(file="iscience.figures/smRCs_alignment_stats_pie.pdf",width = 5,height = 4)
ggplot(mapping_stats, aes(x = "", y = Value, fill = Stat)) + 
  geom_bar(width = 1, stat = "identity", color = "black") +
  coord_polar("y", start = 0)+ ggtitle("") +
  theme_void() + theme(axis.text = element_blank()) + theme(legend.position="bottom") +
  guides(color=guide_legend(ncol=2,override.aes = list(size = 5))) 
dev.off()

### ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





summary(smrc_datasets$mssm_pca$ann_smrc_counts$Length)
summary(smrc_datasets$mssm_pca$ann_smrc_counts$MajorRNAReads)
ix <- log10(smrc_datasets$mssm_pca$ann_smrc_counts$MajorRNAReads)<2
summary((smrc_datasets$mssm_pca$ann_smrc_counts$MajorRNAReads)[ix])






####################################################################################################################
### save.image("/Users/gonzae34/Documents/00.exosome/exo_RData/data_for_iscience_paper_analysis.RData")
### load("/Users/gonzae34/Documents/00.exosome/exo_RData/data_for_iscience_paper_analysis.RData")
####################################################################################################################