
####Salmon/Tximport####

library(Leonardi4.0)

library(GenomicFeatures)
library(tximport)
library(readr)
library(rjson)
library(dplyr)
library(AnnotationDbi)

library(tximportData)
dir <- system.file("extdata", package = "Leonardi4.0", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))
csvfile <- file.path(dir, "sample_table_Mouse_Leonardi.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata

coldata <- coldata[1:32,]
coldata$names <- coldata$Run
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf")
file.exists(coldata$files)


library(tximeta)
se <- tximeta(coldata)
dim(se)
head(rownames(se))
gse <- summarizeToGene(se)
dim(gse)

library(SummarizedExperiment)
data(gse)

gse
assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)

library(pathfindR)

Enrich_Reads1 <- read.csv("/home/john/Documents/CoExposure/PathfindR/CB_1_Sham1_GeneData.csv")
Enrich_Reads2 <- read.csv("/home/john/Documents/CoExposure/PathfindR/CB_4_Sham4_GeneData.csv")
Enrich_Reads3 <- read.csv("/home/john/Documents/CoExposure/PathfindR/O3_1_Sham_1_GeneData.csv")
Enrich_Reads4 <- read.csv("/home/john/Documents/CoExposure/PathfindR/O3_4_Sham_4_GeneData.csv")
Enrich_Reads5 <- read.csv("/home/john/Documents/CoExposure/PathfindR/CB_O3_1_Sham_1_GeneData.csv")
Enrich_Reads6 <- read.csv("/home/john/Documents/CoExposure/PathfindR/CB_O3_4_Sham_4_GeneData.csv")

DeSeq2_Enrich1 <- run_pathfindR(Enrich_Reads1, gene_sets = "mmu_KEGG", output_dir = "pathfindR_Results_CB_1", iterations = 10, n_processes = 32)
DeSeq2_Enrich2 <- run_pathfindR(Enrich_Reads2, gene_sets = "mmu_KEGG", output_dir = "pathfindR_Results_CB_4", iterations = 10, n_processes = 32)
DeSeq2_Enrich3 <- run_pathfindR(Enrich_Reads3, gene_sets = "mmu_KEGG", output_dir = "pathfindR_Results_O3_1", iterations = 10, n_processes = 32)
DeSeq2_Enrich4 <- run_pathfindR(Enrich_Reads4, gene_sets = "mmu_KEGG", output_dir = "pathfindR_Results_O3_4", iterations = 10, n_processes = 32)
DeSeq2_Enrich5 <- run_pathfindR(Enrich_Reads5, gene_sets = "mmu_KEGG", output_dir = "pathfindR_Results_CB-O3_1", iterations = 10, n_processes = 32)
DeSeq2_Enrich6 <- run_pathfindR(Enrich_Reads6, gene_sets = "mmu_KEGG", output_dir = "pathfindR_Results_CB-O3_4", iterations = 10, n_processes = 32)

clustered_1 <- cluster_enriched_terms(DeSeq2_Enrich1, plot_dend = FALSE, plot_clusters_graph = FALSE)
clustered_2 <- cluster_enriched_terms(DeSeq2_Enrich2, plot_dend = FALSE, plot_clusters_graph = FALSE)
clustered_3 <- cluster_enriched_terms(DeSeq2_Enrich3, plot_dend = FALSE, plot_clusters_graph = FALSE)
clustered_4 <- cluster_enriched_terms(DeSeq2_Enrich4, plot_dend = FALSE, plot_clusters_graph = FALSE)
clustered_5 <- cluster_enriched_terms(DeSeq2_Enrich5, plot_dend = FALSE, plot_clusters_graph = FALSE)
clustered_6 <- cluster_enriched_terms(DeSeq2_Enrich6, plot_dend = FALSE, plot_clusters_graph = FALSE)

selected_1 <- subset(clustered_1, Cluster %in% 1:3)
selected_2 <- subset(clustered_2, Cluster %in% 1:3)
selected_3 <- subset(clustered_3, Cluster %in% 1:3)
selected_4 <- subset(clustered_4, Cluster %in% 1:3)
selected_5 <- subset(clustered_5, Cluster %in% 1:3)
selected_6 <- subset(clustered_6, Cluster %in% 1:3)

enrichment_chart(selected_1, plot_by_cluster = TRUE)
enrichment_chart(selected_2, plot_by_cluster = TRUE)
enrichment_chart(selected_3, plot_by_cluster = TRUE)
enrichment_chart(selected_4, plot_by_cluster = TRUE)
enrichment_chart(selected_5, plot_by_cluster = TRUE)
enrichment_chart(selected_6, plot_by_cluster = TRUE)

#term_gene_heatmap(result_df = DeSeq2_Enrich1, genes_df = Enrich_Reads1)
#term_gene_heatmap(result_df = DeSeq2_Enrich2, genes_df = Enrich_Reads2)
#term_gene_heatmap(result_df = DeSeq2_Enrich3, genes_df = Enrich_Reads3)
#term_gene_heatmap(result_df = DeSeq2_Enrich4, genes_df = Enrich_Reads4)
#term_gene_heatmap(result_df = DeSeq2_Enrich5, genes_df = Enrich_Reads5)
#term_gene_heatmap(result_df = DeSeq2_Enrich6, genes_df = Enrich_Reads6)

term_gene_graph(result_df = DeSeq2_Enrich1, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich2, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich3, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich4, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich5, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich6, use_description = TRUE)

term_gene_graph(result_df = DeSeq2_Enrich1, num_terms = 5, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich2, num_terms = 5, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich3, num_terms = 5, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich4, num_terms = 5, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich5, num_terms = 5, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich6, num_terms = 5, use_description = TRUE)

UpSet_plot(result_df = DeSeq2_Enrich1, genes_df = Enrich_Reads1, use_description = TRUE, num_terms = 5, high = "#FF4242", mid = "#E8F086", low = "#6FDE6E")
UpSet_plot(result_df = DeSeq2_Enrich2, genes_df = Enrich_Reads2, use_description = TRUE, num_terms = 5, high = "#FF4242", mid = "#E8F086", low = "#6FDE6E")
UpSet_plot(result_df = DeSeq2_Enrich3, genes_df = Enrich_Reads3, use_description = TRUE, num_terms = 5, high = "#FF4242", mid = "#E8F086", low = "#6FDE6E")
UpSet_plot(result_df = DeSeq2_Enrich4, genes_df = Enrich_Reads4, use_description = TRUE, num_terms = 5, high = "#FF4242", mid = "#E8F086", low = "#6FDE6E")
UpSet_plot(result_df = DeSeq2_Enrich5, genes_df = Enrich_Reads5, use_description = TRUE, num_terms = 5, high = "#FF4242", mid = "#E8F086", low = "#6FDE6E")
UpSet_plot(result_df = DeSeq2_Enrich6, genes_df = Enrich_Reads6, use_description = TRUE, num_terms = 5, high = "#FF4242", mid = "#E8F086", low = "#6FDE6E")

combined_df <- combine_pathfindR_results(result_A = DeSeq2_Enrich1, 
                                         result_B = DeSeq2_Enrich2, 
                                         plot_common = TRUE)
combined_df2 <- combine_pathfindR_results(result_A = DeSeq2_Enrich3, 
                                          result_B = DeSeq2_Enrich4, 
                                          plot_common = TRUE)
combined_df3 <- combine_pathfindR_results(result_A = DeSeq2_Enrich5, 
                                          result_B = DeSeq2_Enrich6, 
                                          plot_common = TRUE)

write.csv(clustered_6, file = "/home/john/Documents/CoExposure/CB_O3_4_clustered.csv")
write.csv(combined_df3, file = "/home/john/Documents/CoExposure/CB-O3_combined.csv")


###DifferentialExpressionAnalysis###

library("DESeq2")
ddsTxi <- DESeqDataSet(gse, design = ~ group)
ddsTxi

keep <- rowSums(counts(ddsTxi) >= 10) >= 4 #(group has to have >5 counts)
dds <- ddsTxi[keep,]

ddsTxi$Exposure <- factor(ddsTxi$Exposure, levels = c("Sham","TiO2")) #The text, condition treated vs untreated, tells you that the estimates are of the logarithmic fold change log2(treated/untreated).
#dds$Exposure <- relevel(dds$Exposure, ref = "Sham")

dds <- DESeq(ddsTxi)
res <- results(dds, name="Exposure_TiO2_vs_Sham")
res
summary(res)


###Log fold change shrinkage for visualization and ranking###

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Exposure_TiO2_vs_Sham", type="apeglm")
resLFC
resOrdered <- res[order(res$pvalue),] #Orders results table by smallest p-value#
summary(res)
sum(res$padj < 0.05, na.rm=TRUE) #Tells us how many adj p-values were less than 0.1#
sum(res$padj < 0.1, na.rm=TRUE)


##Data Normalization Assessment##
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)


#vsddata <- assay(vst(dds, blind = FALSE))
#rlddata <- assay(rlog(dds, blind = FALSE))
#write.csv(vsddata, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/VSD.csv")
#write.csv(rlddata, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/RLD.csv")

library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

#df <- bind_rows(
#  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
#    mutate(transformation = "log2(x + 1)"),
#  as_data_frame(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
#  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)


##Sample Distances##
sampleDists <- dist(t(assay(vsd)))
sampleDists

library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$genotype, vsd$diet, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)


library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds)))

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$genotype, dds$diet, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


##PCA Plot##
#genotype and diet#
plotPCA(vsd, intgroup = c("genotype", "diet"))
pcaData <- plotPCA(vsd, intgroup = c("genotype", "diet"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = genotype, shape = diet)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

#genotype, diet, and sex#
plotPCA(vsd, intgroup = c("genotype", "diet", "sex"))
pcaData <- plotPCA(vsd, intgroup = c("genotype", "diet", "sex"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = genotype, shape = diet, size = sex)) + geom_point() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

write.csv(pcaData, file = "/home/john/Documents/Leonardi_mRNA/PCA2_Dimensions_Leonardi.csv")

library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$genotype <- dds$genotype
gpca.dat$diet <- dds$diet
gpca.dat$sex <- dds$sex
gpca.dat$SampleName <- dds$SampleName

write.csv(gpca.dat, file = "/home/john/Documents/Leonardi_mRNA/PCA_Dimensions_Leonardi.csv")


ggplot(gpca.dat, aes(x = dim1, y = dim2, color = genotype, shape = diet, size = sex)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


##MDS Plot##
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = genotype, shape = diet)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

write.csv(mds, file = "/home/john/Documents/Leonardi_mRNA/MDS_Dimensions_Leonardi.csv")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = genotype, shape = diet)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")

write.csv(mdsPois, file = "/home/john/Documents/Leonardi_mRNA/MDSPOIS_Dimensions_Leonardi.csv")

##Differential Expression Analysis##
#library(ashr)
dds <- DESeq(dds)
resultsNames(dds)

#PairwiseFeb2021
res001 <- results(dds, contrast=c("group","WTMCD","WTFCD"), alpha = 0.05)
mcols(res001, use.names = TRUE)
summary(res001)
table(res001$padj < 0.05)

res002 <- results(dds, contrast=c("group","WTMWD","WTFWD"), alpha = 0.05)
mcols(res002, use.names = TRUE)
summary(res002)
table(res002$padj < 0.05)

res003 <- results(dds, contrast=c("group","WTMWD","WTMCD"), alpha = 0.05)
mcols(res003, use.names = TRUE)
summary(res003)
table(res003$padj < 0.05)

res004 <- results(dds, contrast=c("group","KOMCD","WTMCD"), alpha = 0.05)
mcols(res004, use.names = TRUE)
summary(res004)
table(res004$padj < 0.05)

res005 <- results(dds, contrast=c("group","KOMWD","WTMWD"), alpha = 0.05)
mcols(res005, use.names = TRUE)
summary(res005)
table(res005$padj < 0.05)

res006 <- results(dds, contrast=c("group","KOMCD","KOFCD"), alpha = 0.05)
mcols(res006, use.names = TRUE)
summary(res006)
table(res006$padj < 0.05)

res007 <- results(dds, contrast=c("group","KOMWD","KOFWD"), alpha = 0.05)
mcols(res007, use.names = TRUE)
summary(res007)
table(res007$padj < 0.05)


#Middle group is "treatment", end group is "control"
##General##
res01 <- results(dds, contrast=c("genotype","KO","WT"), alpha = 0.05)
mcols(res01, use.names = TRUE)
summary(res01)
table(res01$padj < 0.05)

res02 <- results(dds, contrast=c("diet","WD","CD"), alpha = 0.05)
mcols(res02, use.names = TRUE)
summary(res02)
table(res02$padj < 0.05)

res03 <- results(dds, contrast=c("sex","Female","Male"), alpha = 0.05)
mcols(res03, use.names = TRUE)
summary(res03)
table(res03$padj < 0.05)

##More Specific##
#Interaction Sets#
#Do individually for male and female
resWT <- results(dds, contrast=c("diet", "WD", "CD"), alpha = 0.05)
resKO <- results(dds, list(c("diet_WD_vs_CD", "genotypeKO.dietWD")), alpha = 0.05)
resI <- results(dds, name="genotypeKO.dietWD", alpha = 0.05)

#resWT2 <- lfcShrink(dds, res = resWT, type="ashr")
#resKO2 <- lfcShrink(dds, res = resKO, type="ashr")
#resI2 <- lfcShrink(dds, res = resI, type="ashr")

ixWT = which.min(resWT$padj) # most significaggplot(pcaData, aes(x = PC1, y = PC2, color = genotype, shape = diet)) +
geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  nt
barplot(assay(dds)[ixWT,],las=2,names=coldata$Run,cex.names=0.9, main=rownames(dds)[ ixWT  ]  )

ixKO = which.min(resKO$padj) # most significant
barplot(assay(dds)[ixKO,],las=2,names=coldata$Run,cex.names=0.9, main=rownames(dds)[ ixKO  ]  )

ixI = which.min(resI$padj) # most significant
barplot(assay(dds)[ixI,],las=2,names=coldata$Run,cex.names=0.9, main=rownames(dds)[ ixI  ]  )
#barplot(assay(dds)[ixI,],las=2,names=coldata$Run,cex.names=0.9,beside=TRUE,col=rep(c("black","white"),each=8), main=rownames(dds)[ ixI  ]  )

##Most Specific##
res1 <- results(dds, contrast=c("group","KOFWD","WTFWD"), alpha = 0.05)
mcols(res1, use.names = TRUE)
summary(res1)
table(res1$padj < 0.05)

res2 <- results(dds, contrast=c("group","KOMWD","WTMWD"), alpha = 0.05)
mcols(res2, use.names = TRUE)
summary(res2)
table(res2$padj < 0.05)

res3 <- results(dds, contrast=c("group","KOFCD","WTFCD"), alpha = 0.05)
mcols(res3, use.names = TRUE)
summary(res3)
table(res3$padj < 0.05)

res4 <- results(dds, contrast=c("group","KOMCD","WTMCD"), alpha = 0.05)
mcols(res4, use.names = TRUE)
summary(res4)
table(res4$padj < 0.05)

#res5 <- results(dds, contrast=c("group","KOFWD","KOMWD"), alpha = 0.05)
#mcols(res5, use.names = TRUE)
#summary(res5)
#table(res5$padj < 0.05)

#res6 <- results(dds, contrast=c("group","WTFWD","WTMWD"), alpha = 0.05)
#mcols(res6, use.names = TRUE)
#summary(res6)
#table(res6$padj < 0.05)

#res7 <- results(dds, contrast=c("group","KOFCD","KOMCD"), alpha = 0.05)
#mcols(res7, use.names = TRUE)
#summary(res7)
#table(res7$padj < 0.05)

#res8 <- results(dds, contrast=c("group","WTFCD","WTMCD"), alpha = 0.05)
#mcols(res8, use.names = TRUE)
#summary(res8)
#table(res8$padj < 0.05)



##VennDiagram##
library(limma)

res1 <- results(dds, contrast=c("group","KOFWD","WTFWD"), alpha = 0.05)
res1 <- res1[which(res1$padj < 0.05),]
res1.genes <- row.names(res1)

res2 <- results(dds, contrast=c("group","KOMWD","WTMWD"), alpha = 0.05)
res2 <- res2[which(res2$padj < 0.05),]
res2.genes <- row.names(res2)

res3 <- results(dds, contrast=c("group","KOFCD","WTFCD"), alpha = 0.05)
res3 <- res3[which(res3$padj < 0.05),]
res3.genes <- row.names(res3)

res4 <- results(dds, contrast=c("group","KOMCD","WTMCD"), alpha = 0.05)
res4 <- res4[which(res4$padj < 0.05),]
res4.genes <- row.names(res4)

#res5 <- results(dds, contrast=c("group","KOFWD","KOMWD"), alpha = 0.05)
#res5 <- res5[which(res5$padj < 0.05),]
#res5.genes <- row.names(res5)

#res6 <- results(dds, contrast=c("group","WTFWD","WTMWD"), alpha = 0.05)
#res6 <- res6[which(res6$padj < 0.05),]
#res6.genes <- row.names(res6)

#res7 <- results(dds, contrast=c("group","KOFCD","KOMCD"), alpha = 0.05)
#res7 <- res7[which(res7$padj < 0.05),]
#res7.genes <- row.names(res7)

#res8 <- results(dds, contrast=c("group","WTFCD","WTMCD"), alpha = 0.05)
#res8 <- res8[which(res8$padj < 0.05),]
#res8.genes <- row.names(res8)

#Unique <- sort(unique(c(res1.genes, res2.genes)))
#Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes)))
Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes, res4.genes)))
#Unique <- sort(unique(c(res5.genes, res6.genes, res7.genes, res8.genes)))
#Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes, res4.genes, res5.genes, res6.genes, res7.genes, res8.genes)))


res1.genes.2 <- Unique %in% res1.genes
res2.genes.2 <- Unique %in% res2.genes
res3.genes.2 <- Unique %in% res3.genes
res4.genes.2 <- Unique %in% res4.genes
#res5.genes.2 <- Unique %in% res5.genes
#res6.genes.2 <- Unique %in% res6.genes
#res7.genes.2 <- Unique %in% res7.genes
#res8.genes.2 <- Unique %in% res8.genes


#counts.1 <- cbind(res1.genes.2,res2.genes.2)
#counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2)
counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2,res4.genes.2)
#counts.1 <- cbind(res5.genes.2,res6.genes.2,res7.genes.2,res8.genes.2)
#counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2,res4.genes.2,res5.genes.2,res6.genes.2,res7.genes.2,res8.genes.2)

results.1 <- vennCounts(counts.1, include="both")

#vennDiagram(results.1, include="both", cex = 1, names = c("",""), circle.col = c("blue", "red"))
#vennDiagram(results.1, include="both", cex = 1, names = c("","",""), circle.col = c("blue", "red", "green"))
vennDiagram(results.1, include="both", cex = 1, names = c("","","",""), circle.col = c("blue", "red", "green", "black"))

#library(UpSetR)
#UpSetHisto <- as.data.frame(array(as.numeric(unlist(results.1)), dim=c(256, 9)))
#upset(UpSetHisto, sets = c("V1", "V2", "V3", "V4", "V5", 
#                           "V6", "V7", "V8"), order.by = "freq", nsets = 8)
#upset(UpSetHisto)


##Other Visualizations##
library(vidger)

vsBoxPlot(dds, d.factor = "group", type = "deseq")
vsDEGMatrix(dds, d.factor = "group", padj = 0.05, type = "deseq")
#vsFourWay("trt2", "trt", "untrt2", d.factor = "group", dds, type = "deseq")
vsMAMatrix(dds, d.factor = "group", type = "deseq", y.lim = c(-10,10))
vsScatterMatrix(dds, d.factor = "group", type = "deseq")
#vsScatterPlot("untrt2", "trt2", dds, d.factor = "group", type = "deseq")
#vsVolcano("untrt2", "trt2", dds, d.factor = "group", type = "deseq", x.lim = c(-10,10), padj = 0.05)
vsVolcanoMatrix(dds, d.factor = "group", type = "deseq", lfc = 2, padj = 0.05, x.lim = c(-8,8),
                title = FALSE, legend = TRUE, grid = TRUE, counts = FALSE, facet.title.size = 10)
##Bubble Chart##
library(pathfindR)
#Need at least 100 significant genes#
Enrich_Reads <- read.csv("/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Leonardi_mRNA/Output/Final/Pairwise 2-19-2021/PathFindR/KOFWDvsKOFCD_Gene_Data.csv")
DeSeq2_Enrich <- run_pathfindR(Enrich_Reads, output_dir = "pathfindR_Results")
enrichment_chart(DeSeq2_Enrich, top_terms = 10, num_bubbles = 4, even_breaks = TRUE) #top_terms= determines the number of most enriched pathways to be included


##Volcano Plots##
library(EnhancedVolcano)

#write.csv(Unique, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_mRNA_unique_1_trt6.csv")
#write.csv(counts.1, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_mRNA_unique_2_trt6.csv")

#resUNIQUE <- read.csv("/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_mRNA_unique_3_trt6.csv", row.names = 1)

#EnhancedVolcano(resUNIQUE, lab = rownames(resUNIQUE), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
#                xlim = c(-10, 8),  ylim = c(0, -log10(10e-11)))

EnhancedVolcano(res01, lab = rownames(res01), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))

EnhancedVolcano(res02, lab = rownames(res02), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))

EnhancedVolcano(res3, lab = rownames(res3), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))pheatmap

EnhancedVolcano(res4, lab = rownames(res4), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-7.5, 7.5))




##Plotting Results##
topGene <- rownames(res01)[which.min(res01$padj)]
plotCounts(dds, gene = topGene, intgroup=c("genotype"))

#GeneofChoice <- "ENSMUSG00000092341"
#plotCounts(dds, gene = GeneofChoice, intgroup=c("genotype"))


library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("genotype","diet"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = genotype, y = count, color = diet)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = genotype, y = count, color = diet, group = diet)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


##MA-Plot##
library("apeglm")
resultsNames(dds)

#res <- lfcShrink(dds, coef="genotype_WT_vs_KO", type="apeglm")
#plotMA(res01, ylim = c(-5, 5))
#res.noshr <- results(dds, name="genotype_WT_vs_KO")
#plotMA(res.noshr, ylim = c(-5, 5))

plotMA(res01, ylim = c(-5,5))
topGene <- rownames(res01)[which.min(res01$padj)]
with(res01[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})


plotMA(res02, ylim = c(-5,5))
topGene <- rownames(res02)[which.min(res02$padj)]
with(res02[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})


plotMA(res03, ylim = c(-10,10))
topGene <- rownames(res03)[which.min(res03$padj)]
with(res03[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

##Histogram of p vaules for genes with mean normalized count larger than 2##
hist(res1$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

##Gene Clustering##

library("genefilter")pheatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
#or#
topGenes <- head(order(res$padj),decreasing = TRUE, 50)

mat  <- assay(vsd)[ topGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("dex","cell")])

my_colour = list(
  cell = c(Four = "#8B4513", One = "#708090"),
  dex = c(trt = "#FF0000", trt2 = "#F08080", trt3 = "#FFFF00", trt4 = "#FFD700",
          trt5 = "#2E8B57", trt6 = "#00FF00", untrt = "#00BFFF", untrt2 = "#1E90FF"))

pheatmap(mat, annotation_col = anno, annotation_colors = my_colour)


#Annotating and Exporting
library("AnnotationDbi")
library("org.Mm.eg.db")

columns(org.Mm.eg.db)

ens.str001 <- rownames(res001)
res001$symbol <- mapIds(org.Mm.eg.db,
                       keys=ens.str001,
                       column="GENENAME",
                       keytype="ENSEMBL",
                       multiVals="first")
res001$entrez <- mapIds(org.Mm.eg.db,
                       keys=ens.str001,
                       column="SYMBOL",
                       keytype="ENSEMBL",
                       multiVals="first")

resOrdered001 <- res001[order(res001$padj),]
head(resOrdered001)

ens.str002 <- rownames(res002)
res002$symbol <- mapIds(org.Mm.eg.db,
                        keys=ens.str002,
                        column="GENENAME",
                        keytype="ENSEMBL",
                        multiVals="first")
res002$entrez <- mapIds(org.Mm.eg.db,
                        keys=ens.str002,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

resOrdered002 <- res002[order(res002$padj),]
head(resOrdered002)

ens.str003 <- rownames(res003)
res003$symbol <- mapIds(org.Mm.eg.db,
                        keys=ens.str003,
                        column="GENENAME",
                        keytype="ENSEMBL",
                        multiVals="first")
res003$entrez <- mapIds(org.Mm.eg.db,
                        keys=ens.str003,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

resOrdered003 <- res003[order(res003$padj),]
head(resOrdered003)

ens.str004 <- rownames(res004)
res004$symbol <- mapIds(org.Mm.eg.db,
                        keys=ens.str004,
                        column="GENENAME",
                        keytype="ENSEMBL",
                        multiVals="first")
res004$entrez <- mapIds(org.Mm.eg.db,
                        keys=ens.str004,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

resOrdered004 <- res004[order(res004$padj),]
head(resOrdered004)

ens.str005 <- rownames(res005)
res005$symbol <- mapIds(org.Mm.eg.db,
                        keys=ens.str005,
                        column="GENENAME",
                        keytype="ENSEMBL",
                        multiVals="first")
res005$entrez <- mapIds(org.Mm.eg.db,
                        keys=ens.str005,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

resOrdered005 <- res005[order(res005$padj),]
head(resOrdered005)

ens.str006 <- rownames(res006)
res006$symbol <- mapIds(org.Mm.eg.db,
                        keys=ens.str006,
                        column="GENENAME",
                        keytype="ENSEMBL",
                        multiVals="first")
res006$entrez <- mapIds(org.Mm.eg.db,
                        keys=ens.str006,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

resOrdered006 <- res006[order(res006$padj),]
head(resOrdered006)

ens.str007 <- rownames(res007)
res007$symbol <- mapIds(org.Mm.eg.db,
                        keys=ens.str007,
                        column="GENENAME",
                        keytype="ENSEMBL",
                        multiVals="first")
res007$entrez <- mapIds(org.Mm.eg.db,
                        keys=ens.str007,
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")

resOrdered007 <- res007[order(res007$padj),]
head(resOrdered007)

##Exporting Results##
resOrderedDF1 <- as.data.frame(resOrdered001)[1:13580, ]
write.csv(resOrderedDF1, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/WTMCDvsWTFCD.csv")

resOrderedDF2 <- as.data.frame(resOrdered002)[1:13580, ]
write.csv(resOrderedDF2, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/WTMWDvsWTFWD.csv")

resOrderedDF3 <- as.data.frame(resOrdered003)[1:13580, ]
write.csv(resOrderedDF3, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/WTMWDvsWTMCD.csv")

resOrderedDF4 <- as.data.frame(resOrdered004)[1:13580, ]
write.csv(resOrderedDF4, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/KOMCDvsWTMCD.csv")
  
resOrderedDF5 <- as.data.frame(resOrdered005)[1:13580, ]
write.csv(resOrderedDF5, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/KOMWDvsWTMWD.csv")

resOrderedDF6 <- as.data.frame(resOrdered006)[1:13580, ]
write.csv(resOrderedDF6, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/KOMCDvsKOFCD.csv")

resOrderedDF7 <- as.data.frame(resOrdered007)[1:13580, ]
write.csv(resOrderedDF7, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/KOMWDvsKOFWD.csv")


##Genomic Space##
resGR <- lfcShrink(dds, coef="genotype_WT_vs_KO", type="apeglm", format="GRanges")
resGR

ens.str <- rownames(res1)
resGR$symbol <- mapIds(org.Mm.eg.db, ens.str, "SYMBOL", "ENSEMBL")

library("Gviz")

window <- resGR[topGene] + 1e6
strand(window) <- "*"
resGRsub <- resGR[resGR %over% window]
naOrDup <- is.na(resGRsub$symbol) | duplicated(resGRsub$symbol)
resGRsub$group <- ifelse(naOrDup, names(resGRsub), resGRsub$symbol)

status <- factor(ifelse(resGRsub$padj < 0.05 & !is.na(resGRsub$padj),
                        "sig", "notsig"))

options(ucscChromosomeNames = FALSE)
g <- GenomeAxisTrack()
a <- AnnotationTrack(resGRsub, name = "gene ranges", feature = status)
d <- DataTrack(resGRsub, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g, d, a), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")



##Gene Ontology##
library(topGO)
library(KEGG.db)

#resUNIQUE2 <- read.csv("/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_mRNA_unique_4_trt6.csv", row.names = 1)

rowsum.threshold <- 1 # user chosen
fdr.threshold <- 0.05 # user chosen
rs <- rowSums(counts(dds))
dds <- dds[ rs > rowsum.threshold ,]
dds <- DESeq(dds)
res <- results(dds, contrast=c("group","KOFWD","KOFCD"), independentFiltering=FALSE, alpha = 0.05) # use count threshold instead of IF
assayed.genes <- rownames(res)
de.genes <- rownames(res)[ which(res$padj < fdr.threshold) ]
gene.vector=as.integer(assayed.genes%in%de.genes)
names(gene.vector)=assayed.genes
head(gene.vector)

library(goseq)
pwf=nullp(gene.vector,"mm9","ensGene")
head(pwf)

GO.wall=goseq(pwf,"mm9","ensGene")
head(GO.wall)

#Random Re-sampling
#GO.samp=goseq(pwf,"mm9","ensGene",method="Sampling",repcnt=1000)
#head(GO.samp)

write.csv(GO.wall, file = "/home/john/Documents/Bioinformatics_Seq_Files/Fastq/Salik_mRNA/Salmon_1.1.0_Output/Salik_GO_Unique_CB-O3_Day4.csv")

KEGG=goseq(pwf,'mm9','ensGene',test.cats="KEGG")
head(KEGG)
