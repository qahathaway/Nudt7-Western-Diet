
####Salmon/Tximport####

library(QuantFiles)

library(GenomicFeatures)
library(tximport)
library(readr)
library(rjson)
library(dplyr)
library(AnnotationDbi)

library(tximportData)
dir <- system.file("extdata", package = "QuantFiles", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))
csvfile <- file.path(dir, "sample_table.csv")
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

Enrich_Reads1 <- read.csv("path/to/files.csv")
Enrich_Reads2 <- read.csv("path/to/files.csv")

DeSeq2_Enrich1 <- run_pathfindR(Enrich_Reads1, gene_sets = "mmu_KEGG", output_dir = "pathfindR_Results_CB_1", iterations = 10, n_processes = 32)
clustered_1 <- cluster_enriched_terms(DeSeq2_Enrich1, plot_dend = FALSE, plot_clusters_graph = FALSE)
selected_1 <- subset(clustered_1, Cluster %in% 1:3)
enrichment_chart(selected_1, plot_by_cluster = TRUE)
term_gene_heatmap(result_df = DeSeq2_Enrich1, genes_df = Enrich_Reads1)
term_gene_graph(result_df = DeSeq2_Enrich1, use_description = TRUE)
term_gene_graph(result_df = DeSeq2_Enrich1, num_terms = 5, use_description = TRUE)
UpSet_plot(result_df = DeSeq2_Enrich1, genes_df = Enrich_Reads1, use_description = TRUE, num_terms = 5, high = "#FF4242", mid = "#E8F086", low = "#6FDE6E")

combined_df <- combine_pathfindR_results(result_A = DeSeq2_Enrich1, 
                                         result_B = DeSeq2_Enrich2, 
                                         plot_common = TRUE)

write.csv(clustered_1, file = "path/to/files.csv")
write.csv(combined_df, file = "path/to/files.csv")


###DifferentialExpressionAnalysis###

library("DESeq2")
ddsTxi <- DESeqDataSet(gse, design = ~ group)
ddsTxi

keep <- rowSums(counts(ddsTxi) >= 10) >= 4 #(group has to have >5 counts)
dds <- ddsTxi[keep,]


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


vsddata <- assay(vst(dds, blind = FALSE))
rlddata <- assay(rlog(dds, blind = FALSE))
write.csv(vsddata, file = "path/to/files.csv")
write.csv(rlddata, file = "path/to/files.csv")

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

write.csv(pcaData, file = "path/to/files.csv")

library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$genotype <- dds$genotype
gpca.dat$diet <- dds$diet
gpca.dat$sex <- dds$sex
gpca.dat$SampleName <- dds$SampleName

write.csv(gpca.dat, file = "path/to/files.csv")


ggplot(gpca.dat, aes(x = dim1, y = dim2, color = genotype, shape = diet, size = sex)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


##MDS Plot##
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = genotype, shape = diet)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

write.csv(mds, file = "path/to/files.csv")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = genotype, shape = diet)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")

write.csv(mdsPois, file = "path/to/files.csv")


##Differential Expression Analysis##
dds <- DESeq(dds)
resultsNames(dds)

#Pairwise
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

res5 <- results(dds, contrast=c("group","KOFWD","KOMWD"), alpha = 0.05)
res5 <- res5[which(res5$padj < 0.05),]
res5.genes <- row.names(res5)

res6 <- results(dds, contrast=c("group","WTFWD","WTMWD"), alpha = 0.05)
res6 <- res6[which(res6$padj < 0.05),]
res6.genes <- row.names(res6)

res7 <- results(dds, contrast=c("group","KOFCD","KOMCD"), alpha = 0.05)
res7 <- res7[which(res7$padj < 0.05),]
res7.genes <- row.names(res7)

res8 <- results(dds, contrast=c("group","WTFCD","WTMCD"), alpha = 0.05)
res8 <- res8[which(res8$padj < 0.05),]
res8.genes <- row.names(res8)

Unique <- sort(unique(c(res1.genes, res2.genes)))
Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes)))
Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes, res4.genes)))
Unique <- sort(unique(c(res5.genes, res6.genes, res7.genes, res8.genes)))
Unique <- sort(unique(c(res1.genes, res2.genes, res3.genes, res4.genes, res5.genes, res6.genes, res7.genes, res8.genes)))


res1.genes.2 <- Unique %in% res1.genes
res2.genes.2 <- Unique %in% res2.genes
res3.genes.2 <- Unique %in% res3.genes
res4.genes.2 <- Unique %in% res4.genes
res5.genes.2 <- Unique %in% res5.genes
res6.genes.2 <- Unique %in% res6.genes
res7.genes.2 <- Unique %in% res7.genes
res8.genes.2 <- Unique %in% res8.genes


counts.1 <- cbind(res1.genes.2,res2.genes.2)
counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2)
counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2,res4.genes.2)
counts.1 <- cbind(res5.genes.2,res6.genes.2,res7.genes.2,res8.genes.2)
counts.1 <- cbind(res1.genes.2,res2.genes.2,res3.genes.2,res4.genes.2,res5.genes.2,res6.genes.2,res7.genes.2,res8.genes.2)

results.1 <- vennCounts(counts.1, include="both")

vennDiagram(results.1, include="both", cex = 1, names = c("",""), circle.col = c("blue", "red"))
vennDiagram(results.1, include="both", cex = 1, names = c("","",""), circle.col = c("blue", "red", "green"))
vennDiagram(results.1, include="both", cex = 1, names = c("","","",""), circle.col = c("blue", "red", "green", "black"))


##Other Visualizations##
library(vidger)

vsBoxPlot(dds, d.factor = "group", type = "deseq")
vsDEGMatrix(dds, d.factor = "group", padj = 0.05, type = "deseq")
vsMAMatrix(dds, d.factor = "group", type = "deseq", y.lim = c(-10,10))
vsScatterMatrix(dds, d.factor = "group", type = "deseq")
vsVolcanoMatrix(dds, d.factor = "group", type = "deseq", lfc = 2, padj = 0.05, x.lim = c(-8,8),
                title = FALSE, legend = TRUE, grid = TRUE, counts = FALSE, facet.title.size = 10)


##Bubble Chart##
library(pathfindR)
#Need at least 100 significant genes#
Enrich_Reads <- read.csv("path/to/files.csv")
DeSeq2_Enrich <- run_pathfindR(Enrich_Reads, output_dir = "pathfindR_Results")
enrichment_chart(DeSeq2_Enrich, top_terms = 10, num_bubbles = 4, even_breaks = TRUE) #top_terms= determines the number of most enriched pathways to be included


##Volcano Plots##
library(EnhancedVolcano)

write.csv(Unique, file = "path/to/files.csv")
write.csv(counts.1, file = "path/to/files.csv")

resUNIQUE <- read.csv("path/to/files.csv", row.names = 1)

EnhancedVolcano(resUNIQUE, lab = rownames(resUNIQUE), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-10, 8),  ylim = c(0, -log10(10e-11)))

EnhancedVolcano(res01, lab = rownames(res01), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
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


##Exporting Results##
resOrderedDF1 <- as.data.frame(resOrdered001)[1:13580, ]
write.csv(resOrderedDF1, file = "/home/john/Documents/Leonardi_mRNA/Leonardi_01272022/WTMCDvsWTFCD.csv")


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

resUNIQUE2 <- read.csv("path/to/files.csv", row.names = 1)

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

write.csv(GO.wall, file = "path/to/files.csv")

KEGG=goseq(pwf,'mm9','ensGene',test.cats="KEGG")
head(KEGG)
