# Setting path to the working directory
setwd("C:/Users/aaisu/Box/Carter Lab/Projects/RNAseq_EthiopiaNUCI")
## Written on R version 4.5.1

# Importing necessary libraries
library(rtracklayer)
library(GenomicFeatures)
library(DESeq2)
library(dplyr)
library(tibble)
library(ggplot2)
library(tximport)
library(readr)
library(AnnotationDbi)     # if you need to map IDs
library(apeglm)            # for LFC shrinkage
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(AnnotationForge)
library(txdbmaker)
library(biomaRt)
library(org.Astephensi.eg.db)
library(GO.db)

# 1. Read in sample metadata
# Expect a table with columns: sampleName, stage, quant_dir
sampleTable <- read_csv("./quant_files/All_sample_info.csv")
# e.g. sample_info.csv:
# sampleName,stage,quant_dir
# S1,L1,/path/to/S1
# S2,L1,/path/to/S2
# ...
sampleTable$stage <- factor(sampleTable$stage, levels=c("IS1","IS2","IS3","IS4", "DK"))
sampleTable$Collection <- factor(sampleTable$Collection, levels = c("Lab", "Wild"))
sampleTable$Site <- factor(sampleTable$Site, levels = c("UCI", "Adama", "Erer", "Jijiga") )
rownames(sampleTable) <- sampleTable$sampleName

# 2. Build vector of Salmon quant.sf files
files <- file.path(sampleTable$quant_dir, "quant.sf")
names(files) <- sampleTable$sampleName

# 3. Import transcript-level counts and aggregate to genes
# You need a two-column data.frame tx2gene mapping transcripts to genes
transcript_info <- read.csv("./Reference_transcriptome/Anstep_UCI_V1.0_transcripts_info.csv")
# Each transcript ID ends with ".1" extension indicating the version. However, this has been disregarded in the Salmon quantification.
# Therefore I had to remove it from the tx2gene table
# Note: Remember to add that part before any conversion of transcript IDs to another format.
# For the KEGGS analysis we need the gene ID without the perfix "LOC". Therefore it was also removed
transcript_info$accession <- sub("\\.1$", "", transcript_info$accession)
transcript_info$gene_id <- sub("^LOC", "", transcript_info$gene_id)
tx2gene <- transcript_info[,c(1,5)]
txi <- tximport(files,
                type="salmon",
                tx2gene=tx2gene,
                ignoreTxVersion=TRUE)

# 4. Create DESeqDataSet
# In this step I'm creating 4 data sets corresponding to the 4 larval stages
# This is for the purpose of using each larval stage as reference factor
dds <- DESeqDataSetFromTximport(txi,
                                   colData = sampleTable,
                                   design = ~ Site)
ddsIS1 <- DESeqDataSetFromTximport(txi,
                                colData = sampleTable,
                                design = ~ stage)
ddsIS2 <- DESeqDataSetFromTximport(txi,
                                   colData = sampleTable,
                                   design = ~ stage)
ddsIS3 <- DESeqDataSetFromTximport(txi,
                                   colData = sampleTable,
                                   design = ~ stage)
ddsIS4 <- DESeqDataSetFromTximport(txi,
                                   colData = sampleTable,
                                   design = ~ stage)
ddsSiteJijiga <- DESeqDataSetFromTximport(txi,
                                   colData = sampleTable,
                                   design = ~ Site)
ddsSiteErer <- DESeqDataSetFromTximport(txi,
                                          colData = sampleTable,
                                          design = ~ Site)
ddsSiteAdama <- DESeqDataSetFromTximport(txi,
                                          colData = sampleTable,
                                          design = ~ Site)
ddsSiteUCI <- DESeqDataSetFromTximport(txi,
                                         colData = sampleTable,
                                         design = ~ Site)
ddsCollection <- DESeqDataSetFromTximport(txi,
                                    colData = sampleTable,
                                    design = ~ Collection)

# 5. Prefilter low-count genes (optional but recommended)
# Here I'm using the same Keep input to filter the genes from dds object I created only for this purpose
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
ddsIS1 <- ddsIS1[keep,]
ddsIS2 <- ddsIS2[keep,]
ddsIS3 <- ddsIS3[keep,]
ddsIS4 <- ddsIS4[keep,]
ddsSiteJijiga <- ddsSiteJijiga[keep,]
ddsSiteErer <- ddsSiteErer[keep,]
ddsSiteAdama <- ddsSiteAdama[keep,]
ddsSiteUCI <- ddsSiteUCI[keep,]
ddsCollection <- ddsCollection[keep,]

# 5.5 Setting the reference level for each larval stage separately
ddsIS1$stage <- relevel(ddsIS1$stage, ref = "IS1")
ddsIS2$stage <- relevel(ddsIS1$stage, ref = "IS2")
ddsIS3$stage <- relevel(ddsIS1$stage, ref = "IS3")
ddsIS4$stage <- relevel(ddsIS1$stage, ref = "IS4")
ddsSiteJijiga$Site <- relevel(ddsSiteJijiga$Site, ref = "Jijiga")
ddsSiteErer$Site <- relevel(ddsSiteErer$Site, ref = "Erer")
ddsSiteAdama$Site <- relevel(ddsSiteAdama$Site, ref = "Adama")
ddsSiteUCI$Site <- relevel(ddsSiteAdama$Site, ref = "UCI")
ddsCollection$Collection <- relevel(ddsCollection$Collection, ref = "Lab")

# 6. Run the DESeq pipeline
ddsIS1 <- DESeq(ddsIS1)
ddsIS2 <- DESeq(ddsIS2)
ddsIS3 <- DESeq(ddsIS3)
ddsIS4 <- DESeq(ddsIS4)
ddsSiteJijiga <- DESeq(ddsSiteJijiga)
ddsSiteErer <- DESeq(ddsSiteErer)
ddsSiteAdama <- DESeq(ddsSiteAdama)
ddsSiteUCI <- DESeq(ddsSiteUCI)
ddsCollection <- DESeq(ddsCollection)

# 7. Examine results names (coefficients)
resultsNames(ddsIS1)
resultsNames(ddsIS2)
resultsNames(ddsIS3)
resultsNames(ddsIS4)
resultsNames(ddsSiteJijiga)
resultsNames(ddsSiteErer)
resultsNames(ddsSiteAdama)
resultsNames(ddsSiteUCI)
resultsNames(ddsCollection)
# e.g. "Intercept", "stage_L2_vs_L1", "stage_L3_vs_L1", "stage_L4_vs_L1"

# 8. Extract results of important comparisons. It is important to pick the correct model for each comparison
# The reason behind having multiple modles and picking different models for each comparison is, the set of upregulated vs. downregulated genes
# depends on the reference larval. For example if we want get the upregulated genes in IS3 compared to IS2 we have to pick the model ddsIS2
# because this model has the IS2 as the reference and in the coefficient stage_IS3_vs_IS2 willhave the upregulated and downregulated genes in IS3 against IS2.
res_IS2_vs_IS1 <- results(ddsIS1, contrast=c("stage","IS2","IS1"))
res_IS3_vs_IS2 <- results(ddsIS2, contrast=c("stage","IS3","IS2"))
res_IS4_vs_IS3 <- results(ddsIS3, contrast=c("stage","IS4","IS3"))
res_Erer_vs_Jijiga <- results(ddsSiteJijiga, contrast = c("Site", "Erer", "Jijiga"))
res_Adama_vs_Jijiga <- results(ddsSiteJijiga, contrast = c("Site", "Adama", "Jijiga"))
res_Jijiga_vs_Erer <- results(ddsSiteErer, contrast = c("Site", "Jijiga", "Erer"))
res_Adama_vs_Erer <- results(ddsSiteErer, contrast = c("Site", "Adama", "Erer"))
res_Jijiga_vs_Adama <- results(ddsSiteAdama, contrast = c("Site", "Jijiga", "Adama"))
res_Erer_vs_Adama <- results(ddsSiteAdama, contrast = c("Site", "Erer", "Adama"))
res_Jijiga_vs_UCI <- results(ddsSiteUCI, contrast = c("Site", "Jijiga", "UCI"))
res_Erer_vs_UCI <- results(ddsSiteUCI, contrast = c("Site", "Erer", "UCI"))
res_Adama_vs_UCI <- results(ddsSiteUCI, contrast = c("Site", "Adama", "UCI"))
res_Wild_vs_Lab <- results(ddsCollection, contrast = c("Collection", "Wild", "Lab"))

# 9. Shrink log2 fold-changes for more accurate effect sizes. It is very important to pick the correct model.
resLFC_IS2_vs_IS1 <- lfcShrink(ddsIS1,
                             coef = 2,
                             type="apeglm")
resLFC_IS3_vs_IS2 <- lfcShrink(ddsIS2,
                               coef = 3,
                             type="apeglm")
resLFC_IS4_vs_IS3 <- lfcShrink(ddsIS3,
                               coef = 4,
                             type="apeglm")
resLFC_Erer_vs_Jijiga <- lfcShrink(ddsSiteJijiga,
                               coef = 4,
                               type="apeglm")
resLFC_Adama_vs_Jijiga <- lfcShrink(ddsSiteJijiga,
                                   coef = 3,
                                   type="apeglm")
resLFC_Jijiga_vs_Erer <- lfcShrink(ddsSiteErer,
                                    coef = 4,
                                    type="apeglm")
resLFC_Adama_vs_Erer <- lfcShrink(ddsSiteErer,
                                   coef = 3,
                                   type="apeglm")
resLFC_Jijiga_vs_Adama <- lfcShrink(ddsSiteAdama,
                                  coef = 4,
                                  type="apeglm")
resLFC_Erer_vs_Adama <- lfcShrink(ddsSiteAdama,
                                    coef = 3,
                                    type="apeglm")
resLFC_Adama_vs_UCI <- lfcShrink(ddsSiteUCI,
                                  coef = 2,
                                  type="apeglm")
resLFC_Erer_vs_UCI <- lfcShrink(ddsSiteUCI,
                                 coef = 3,
                                 type="apeglm")
resLFC_Jijiga_vs_UCI <- lfcShrink(ddsSiteUCI,
                                coef = 4,
                                type="apeglm")
resLFC_Wild_vs_Lab <- lfcShrink(ddsCollection,
                                  coef = 2,
                                  type="apeglm")

# 10. Quick summaries
summary(resLFC_IS2_vs_IS1)
summary(resLFC_IS3_vs_IS2)
summary(resLFC_IS4_vs_IS3)
summary(resLFC_Erer_vs_Jijiga)
summary(resLFC_Adama_vs_Jijiga)
summary(resLFC_Jijiga_vs_Erer)
summary(resLFC_Adama_vs_Erer)
summary(resLFC_Jijiga_vs_Adama)
summary(resLFC_Erer_vs_Adama)
summary(resLFC_Jijiga_vs_UCI)
summary(resLFC_Erer_vs_UCI)
summary(resLFC_Adama_vs_UCI)
summary(resLFC_Wild_vs_Lab)

# 11. Export top tables (e.g. padj < 0.05 & |log2FC| > 1)
sigUp_IS2_vs_IS1 <- subset(resLFC_IS2_vs_IS1, padj < 0.05 & log2FoldChange > 1)
sigDown_IS2_vs_IS1 <- subset(resLFC_IS2_vs_IS1, padj < 0.05 & log2FoldChange < 1)
sigUp_IS3_vs_IS2 <- subset(resLFC_IS3_vs_IS2, padj < 0.05 & log2FoldChange > 1)
sigDown_IS3_vs_IS2 <- subset(resLFC_IS3_vs_IS2, padj < 0.05 & log2FoldChange < 1)
sigUp_IS4_vs_IS3 <- subset(resLFC_IS4_vs_IS3, padj < 0.05 & log2FoldChange > 1)
sigDown_IS4_vs_IS3 <- subset(resLFC_IS4_vs_IS3, padj < 0.05 & log2FoldChange < 1)
sigUp_Erer_vs_Jijiga <- subset(resLFC_Erer_vs_Jijiga, padj < 0.05 & log2FoldChange > 1)
sigDown_Erer_vs_Jijiga <- subset(resLFC_Erer_vs_Jijiga, padj < 0.05 & log2FoldChange < 1)
sigUp_Adama_vs_Jijiga <- subset(resLFC_Adama_vs_Jijiga, padj < 0.05 & log2FoldChange > 1)
sigDown_Adama_vs_Jijiga <- subset(resLFC_Adama_vs_Jijiga, padj < 0.05 & log2FoldChange < 1)
sigUp_Jijiga_vs_Erer <- subset(resLFC_Jijiga_vs_Erer, padj < 0.05 & log2FoldChange > 1)
sigDown_Jijiga_vs_Erer <- subset(resLFC_Jijiga_vs_Erer, padj < 0.05 & log2FoldChange < 1)
sigUp_Adama_vs_Erer <- subset(resLFC_Adama_vs_Erer, padj < 0.05 & log2FoldChange > 1)
sigDown_Adama_vs_Erer <- subset(resLFC_Adama_vs_Erer, padj < 0.05 & log2FoldChange < 1)
sigUp_Jijiga_vs_Adama <- subset(resLFC_Jijiga_vs_Adama, padj < 0.05 & log2FoldChange > 1)
sigDown_Jijiga_vs_Adama <- subset(resLFC_Jijiga_vs_Adama, padj < 0.05 & log2FoldChange < 1)
sigUp_Erer_vs_Adama <- subset(resLFC_Erer_vs_Adama, padj < 0.05 & log2FoldChange > 1)
sigDown_Erer_vs_Adama <- subset(resLFC_Erer_vs_Adama, padj < 0.05 & log2FoldChange < 1)
sigUp_Jijiga_vs_UCI <- subset(resLFC_Jijiga_vs_UCI, padj < 0.05 & log2FoldChange > 1)
sigDown_Jijiga_vs_UCI <- subset(resLFC_Jijiga_vs_UCI, padj < 0.05 & log2FoldChange < 1)
sigUp_Erer_vs_UCI <- subset(resLFC_Erer_vs_UCI, padj < 0.05 & log2FoldChange > 1)
sigDown_Erer_vs_UCI <- subset(resLFC_Erer_vs_UCI, padj < 0.05 & log2FoldChange < 1)
sigUp_Adama_vs_UCI <- subset(resLFC_Adama_vs_UCI, padj < 0.05 & log2FoldChange > 1)
sigDown_Adama_vs_UCI <- subset(resLFC_Adama_vs_UCI, padj < 0.05 & log2FoldChange < 1)
sigUp_Wild_vs_Lab <- subset(resLFC_Wild_vs_Lab, padj < 0.05 & log2FoldChange > 1)
sigDown_Wild_vs_Lab <- subset(resLFC_Wild_vs_Lab, padj < 0.05 & log2FoldChange < 1)

# write.csv(as.data.frame(sig_L2_vs_L1), file="DEG_L2_vs_L1.csv")

# 12. PCA of samples on VST-transformed counts
vsdIS1 <- vst(ddsIS1, blind=FALSE)
vsdIS2 <- vst(ddsIS2, blind=FALSE)
vsdIS3 <- vst(ddsIS3, blind=FALSE)
vsdSite <- vst(ddsCollection, blind = FALSE)

# Plotting PCA
stages_PCA_plot <- plotPCA(vsdIS3, intgroup="stage") +
  ggtitle("PCA of An. stephensi larval stages")

Site_PCA_plot <- plotPCA(vsdSite, intgroup=c("Site", "stage")) +
  ggtitle("PCA of An. stephensi by site")

png(file = paste("./UCI_Diff_Exp_plots/Anstep_Larval_stages_PCA_plot.png"), width = 10, height = 5, units = "in", res = 300)
print(stages_PCA_plot)
dev.off()

# 13. Sample-to-sample distance heatmap
sampleDists <- dist(t(assay(vsdSite)))
mat <- as.matrix(sampleDists)
rownames(mat) <- colnames(mat) <- sampleTable$sampleName
pheatmap(mat,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main="Sample distances")

# 14. Optional: export normalized counts
norm_counts <- counts(dds, normalized=TRUE)
# write.csv(as.data.frame(norm_counts), file="normalized_counts.csv")

# 15. Prepare your gene list
#    Here we assume `sig_IS2_vs_IS1` is the DESeq2 result
#    with rownames = Entrez gene IDs. If yours are SYMBOLs or
#    other IDs, you must map them to Entrez first (see note below).
sigUp_IS2_vs_IS1_gene_ids <- rownames(sigUp_IS2_vs_IS1[ order(sigUp_IS2_vs_IS1$log2FoldChange, decreasing = TRUE), ]) 
sigUp_IS3_vs_IS2_gene_ids <- rownames(sigUp_IS3_vs_IS2[ order(sigUp_IS3_vs_IS2$log2FoldChange, decreasing = TRUE), ])
sigUp_IS4_vs_IS3_gene_ids <- rownames(sigUp_IS4_vs_IS3[ order(sigUp_IS4_vs_IS3$log2FoldChange, decreasing = TRUE), ])
sigDown_IS2_vs_IS1_gene_ids <- rownames(sigDown_IS2_vs_IS1[ order(sigDown_IS2_vs_IS1$log2FoldChange, decreasing = TRUE), ])
sigDown_IS3_vs_IS2_gene_ids <- rownames(sigDown_IS3_vs_IS2[ order(sigDown_IS3_vs_IS2$log2FoldChange, decreasing = TRUE), ])
sigDown_IS4_vs_IS3_gene_ids <- rownames(sigDown_IS4_vs_IS3[ order(sigDown_IS4_vs_IS3$log2FoldChange, decreasing = TRUE), ])
sigUp_Erer_vs_Jijiga_gene_ids <- rownames(sigUp_Erer_vs_Jijiga[ order(sigUp_Erer_vs_Jijiga$log2FoldChange, decreasing = TRUE), ])
sigDown_Erer_vs_Jijiga_gene_ids <- rownames(sigDown_Erer_vs_Jijiga[ order(sigDown_Erer_vs_Jijiga$log2FoldChange, decreasing = TRUE), ])
sigUp_Adama_vs_Jijiga_gene_ids <- rownames(sigUp_Adama_vs_Jijiga[ order(sigUp_Adama_vs_Jijiga$log2FoldChange, decreasing = TRUE), ])
sigDown_Adama_vs_Jijiga_gene_ids <- rownames(sigDown_Adama_vs_Jijiga[ order(sigDown_Adama_vs_Jijiga$log2FoldChange, decreasing = TRUE), ])
sigUp_Jijiga_vs_Erer_gene_ids <- rownames(sigUp_Jijiga_vs_Erer[ order(sigUp_Jijiga_vs_Erer$log2FoldChange, decreasing = TRUE), ])
sigDown_Jijiga_vs_Erer_gene_ids <- rownames(sigDown_Jijiga_vs_Erer[ order(sigDown_Jijiga_vs_Erer$log2FoldChange, decreasing = TRUE), ])
sigUp_Adama_vs_Erer_gene_ids <- rownames(sigUp_Adama_vs_Erer[ order(sigUp_Adama_vs_Erer$log2FoldChange, decreasing = TRUE), ])
sigDown_Adama_vs_Erer_gene_ids <- rownames(sigDown_Adama_vs_Erer[ order(sigDown_Adama_vs_Erer$log2FoldChange, decreasing = TRUE), ])
sigUp_Jijiga_vs_Adama_gene_ids <- rownames(sigUp_Jijiga_vs_Adama[ order(sigUp_Jijiga_vs_Adama$log2FoldChange, decreasing = TRUE), ])
sigDown_Jijiga_vs_Adama_gene_ids <- rownames(sigDown_Jijiga_vs_Adama[ order(sigDown_Jijiga_vs_Adama$log2FoldChange, decreasing = TRUE), ])
sigUp_Erer_vs_Adama_gene_ids <- rownames(sigUp_Erer_vs_Adama[ order(sigUp_Erer_vs_Adama$log2FoldChange, decreasing = TRUE), ])
sigDown_Erer_vs_Adama_gene_ids <- rownames(sigDown_Erer_vs_Adama[ order(sigDown_Erer_vs_Adama$log2FoldChange, decreasing = TRUE), ])
sigUp_Jijiga_vs_UCI_gene_ids <- rownames(sigUp_Jijiga_vs_UCI[ order(sigUp_Jijiga_vs_UCI$log2FoldChange, decreasing = TRUE), ])
sigDown_Jijiga_vs_UCI_gene_ids <- rownames(sigDown_Jijiga_vs_UCI[ order(sigDown_Jijiga_vs_UCI$log2FoldChange, decreasing = TRUE), ])
sigUp_Erer_vs_UCI_gene_ids <- rownames(sigUp_Erer_vs_UCI[ order(sigUp_Erer_vs_UCI$log2FoldChange, decreasing = TRUE), ])
sigDown_Erer_vs_UCI_gene_ids <- rownames(sigDown_Erer_vs_UCI[ order(sigDown_Erer_vs_UCI$log2FoldChange, decreasing = TRUE), ])
sigUp_Adama_vs_UCI_gene_ids <- rownames(sigUp_Adama_vs_UCI[ order(sigUp_Adama_vs_UCI$log2FoldChange, decreasing = TRUE), ])
sigDown_Adama_vs_UCI_gene_ids <- rownames(sigDown_Adama_vs_UCI[ order(sigDown_Adama_vs_UCI$log2FoldChange, decreasing = TRUE), ])
sigUp_Wild_vs_Lab_gene_ids <- rownames(sigUp_Wild_vs_Lab[ order(sigUp_Wild_vs_Lab$log2FoldChange, decreasing = TRUE), ])
sigDown_Wild_vs_Lab_gene_ids <- rownames(sigDown_Wild_vs_Lab[ order(sigDown_Wild_vs_Lab$log2FoldChange, decreasing = TRUE), ])

# Here I'm creating a vector to hold the names of the vairables containing gene IDs to be used in the for loops.
sigGeneVars <- c("sigUp_IS2_vs_IS1_gene_ids", "sigUp_IS3_vs_IS2_gene_ids", "sigUp_IS4_vs_IS3_gene_ids",
                 "sigDown_IS2_vs_IS1_gene_ids", "sigDown_IS3_vs_IS2_gene_ids", "sigDown_IS4_vs_IS3_gene_ids",
                 "sigUp_Erer_vs_Jijiga_gene_ids", "sigDown_Erer_vs_Jijiga_gene_ids",
                 "sigUp_Adama_vs_Jijiga_gene_ids", "sigDown_Adama_vs_Jijiga_gene_ids",
                 "sigUp_Jijiga_vs_Erer_gene_ids", "sigDown_Jijiga_vs_Erer_gene_ids",
                 "sigUp_Adama_vs_Erer_gene_ids", "sigDown_Adama_vs_Erer_gene_ids",
                 "sigUp_Jijiga_vs_Adama_gene_ids", "sigDown_Jijiga_vs_Adama_gene_ids",
                 "sigUp_Erer_vs_Adama_gene_ids", "sigDown_Erer_vs_Adama_gene_ids",
                 "sigUp_Jijiga_vs_UCI_gene_ids", "sigDown_Jijiga_vs_UCI_gene_ids",
                 "sigUp_Erer_vs_UCI_gene_ids", "sigDown_Erer_vs_UCI_gene_ids",
                 "sigUp_Adama_vs_UCI_gene_ids", "sigDown_Adama_vs_UCI_gene_ids",
                 "sigUp_Wild_vs_Lab_gene_ids", "sigDown_Wild_vs_Lab_gene_ids")

# 16. Run KEGG over-representation analysis
# Here I will be running the KEGGS analysis and will be generating the plots in the same loop.

# Looping over the variable names
for (sgv in sigGeneVars){
  sigGeneIDs <- get(sgv)
  kegg_res <- enrichKEGG(
    gene         = sigGeneIDs,
    organism     = "aste",         # An. stephensi KEGG code :contentReference[oaicite:1]{index=1}
    keyType      = "ncbi-geneid",         # assuming your IDs are NCBI-GeneIDs
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.10
  )
  
  # 4. Inspect results
  # print(head(as.data.frame(kegg_res)))
  
  # Saving the KEGGS results for each gene set in CSV file
  # This replaces "_gene_ids" with nothing, e.g., "sigUp_IS2_vs_IS1_gene_ids" -> "sigUp_IS2_vs_IS1"
  file_identifier <- gsub("_gene_ids", "", sgv)
  # Convert the enrichKEGG result to a data frame
  kegg_results_df <- as.data.frame(kegg_res)
  
  if (nrow(kegg_results_df) > 0) {
    # Define the output file path for the CSV
    output_filename <- paste0("./KEGGS_results/", file_identifier, "_KEGG_results.csv")
    
    # Write the data frame to a CSV file
    write.csv(kegg_results_df, file = output_filename, row.names = FALSE)
  }
  
  # Getting gene list identifier
  genListIDElements <- strsplit(sgv, "_gene_ids")
  
  # 5. Visualize top pathways
  # 5a. Barplot
  keggBarPlot <- barplot(kegg_res, showCategory=20, title=paste("Top KEGG pathways", genListIDElements[1], sep = " " ))
  png(file = paste("./UCI_Diff_Exp_plots/", genListIDElements[1], "KEGGS_PW_BP.png"), width = 10, height = 10, units = "in", res = 300)
  print(keggBarPlot)
  dev.off()
  
  # 5b. Dotplot
  keggdotplot <- dotplot(kegg_res, showCategory=20) + ggtitle(paste("Top KEGG pathways", genListIDElements[1], sep = " " ))
  png(file = paste("./UCI_Diff_Exp_plots/", genListIDElements[1], "KEGGS_PW_DP.png"), width = 10, height = 10, units = "in", res = 300)
  print(keggdotplot)
  dev.off()
  
  # 5c. Pathway‐gene network
  # cnetplot(kegg_res_IS2_vs_IS1, categorySize="pvalue", foldChange=sig_IS2_vs_IS1$log2FoldChange)
}

# 17. Running GSEA analysis
# I have created a org.db package for anopheles stephensi named "org.Astephensi.eg.db" and I need to have it loaded into the environment to run the rest of the analysis.

# Before running GSEA taking a look at keytypes in org.db
keytypes(org.Astephensi.eg.db)

# Set the org.db
orgdb <- org.Astephensi.eg.db

# Pull GO annotations (direct only) and keep BP
ann <- AnnotationDbi::select(
  orgdb,
  keys = keys(orgdb, keytype = "GID"), # or whatever ID you want to use
  columns = c("GO","ONTOLOGY"),
  keytype = "GID"
)
ann <- ann[!is.na(ann$GO) & ann$ONTOLOGY == "BP", c("GO","GID")]
colnames(ann) <- c("term", "gene")  # TERM2GENE

# 2) Optional TERM2NAME (GO term names)
term_names <- AnnotationDbi::select(GO.db, keys = unique(ann$term),
                                    columns = "TERM", keytype = "GOID")
colnames(term_names) <- c("term","name")  # TERM2NAME

# First I need to get the ranked gene lists
# Ensure rownames are gene IDs (must be Entrez IDs or mapped to them)
resLFC_IS2_vs_IS1 <- as.data.frame(resLFC_IS2_vs_IS1)
resLFC_IS2_vs_IS1 <- rownames_to_column(resLFC_IS2_vs_IS1, var = "gene")

resLFC_IS3_vs_IS2 <- as.data.frame(resLFC_IS3_vs_IS2)
resLFC_IS3_vs_IS2 <- rownames_to_column(resLFC_IS3_vs_IS2, var = "gene")

resLFC_IS4_vs_IS3 <- as.data.frame(resLFC_IS4_vs_IS3)
resLFC_IS4_vs_IS3 <- rownames_to_column(resLFC_IS4_vs_IS3, var = "gene")

# Drop NA adjusted p-values and log2FCs
resLFC_IS2_vs_IS1 <- resLFC_IS2_vs_IS1 %>%
  filter(!is.na(log2FoldChange), !is.na(padj))

resLFC_IS3_vs_IS2 <- resLFC_IS3_vs_IS2 %>%
  filter(!is.na(log2FoldChange), !is.na(padj))

resLFC_IS4_vs_IS3 <- resLFC_IS4_vs_IS3 %>%
  filter(!is.na(log2FoldChange), !is.na(padj))

# Create geneList named vector
geneList_IS2_vs_IS1 <- resLFC_IS2_vs_IS1$log2FoldChange
names(geneList_IS2_vs_IS1) <- resLFC_IS2_vs_IS1$gene

geneList_IS3_vs_IS2 <- resLFC_IS3_vs_IS2$log2FoldChange
names(geneList_IS3_vs_IS2) <- resLFC_IS3_vs_IS2$gene

geneList_IS4_vs_IS3 <- resLFC_IS4_vs_IS3$log2FoldChange
names(geneList_IS4_vs_IS3) <- resLFC_IS4_vs_IS3$gene

# Sort in decreasing order
geneList_IS2_vs_IS1 <- sort(geneList_IS2_vs_IS1, decreasing = TRUE)

geneList_IS3_vs_IS2 <- sort(geneList_IS3_vs_IS2, decreasing = TRUE)

geneList_IS4_vs_IS3 <- sort(geneList_IS4_vs_IS3, decreasing = TRUE)

# Running a for loop for the three gene lists to perform GSEA

geneLists <- c("geneList_IS2_vs_IS1", "geneList_IS3_vs_IS2", "geneList_IS4_vs_IS3")

for (gl in geneLists){
  
  # Getting the gene list
  gene_list <- get(gl)

  # Running GSEA
  gsea_result <- GSEA(geneList = gene_list,
                      TERM2GENE = ann,
                      TERM2NAME = term_names,
                      pvalueCutoff = 0.05,
                      minGSSize = 10,
                      maxGSSize = 500,
                      verbose = TRUE)
  
  # printing the head of gsea results table
  print(gsea_result@result)
  print(gl)
  
  # Dotplot
  png(file = paste("./UCI_Diff_Exp_plots/", gl, "_GSEA_dotplot.png", sep = ""), width = 10, height = 10, units = "in", res = 300)
  print(dotplot(gsea_result, showCategory = 20, title = paste(gl, "GSEA dotplot", sep = " ")))
  dev.off()
  
  # Ridgeplot
  png(file = paste("./UCI_Diff_Exp_plots/", gl, "_GSEA_ridgeplot.png", sep = ""), width = 10, height = 12, units = "in", res = 300)
  print(ridgeplot(gsea_result) + ggtitle(paste(gl, "GSEA Ridge Plot", sep = " ")))
  dev.off()
  
  # Enrichment Map
  #emapplot(pairwise_termsim(gsea_result), showCategory = 20)
  
  # GSEA plot for a specific GO term
  #gseaplot2(gsea_result, geneSetID = gsea_result$ID[1])

}

# 18. Running Over representation analysis for significantly up and down reguated genes between stages

# Looping over the variable names
for (sgv in sigGeneVars){
  
  # Separating variable name elements to produce experiment name
  var_name_elements <- strsplit(sgv, "_")
  #print(var_name_elements[[1]][1])
  
  # Preparing universe input
  uni_input <- get(paste("res", var_name_elements[[1]][2], var_name_elements[[1]][3], var_name_elements[[1]][4], sep = "_"))
  #print(head(uni_input))
  
  ## Define background (universe): genes that were actually tested
  universe <- rownames(uni_input)[!is.na(uni_input$padj)]
  
  enricher_res <- enricher(
    gene          = get(sgv),
    TERM2GENE     = ann,
    TERM2NAME     = term_names,
    universe      = universe,
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.05,
    minGSSize     = 10,
    maxGSSize     = 500
  )
  
  ## Visualizing top 20 terms – ORA
  png(file = paste("./UCI_Diff_Exp_plots/", sgv, "_enricher_GO_dotplot.png", sep = ""), width = 10, height = 10, units = "in", res = 300)
  print(dotplot(enricher_res, showCategory = 20))
  dev.off()
  
  png(file = paste("./UCI_Diff_Exp_plots/", sgv, "_enricher_GO_barplot.png", sep = ""), width = 10, height = 10, units = "in", res = 300)
  print(barplot(enricher_res, showCategory = 20))
  dev.off()

}


