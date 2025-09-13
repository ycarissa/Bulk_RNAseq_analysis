# Load libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(AnnotationHub)
library(clusterProfiler)
library(org.Mm.eg.db)


# Directories
work_dir <- "C:/Users/caris/Desktop/Projects/GSE205459"
save_dir <- file.path(work_dir, "Analysis")
dir.create(save_dir)

cts <- read.table(file.path(work_dir, "GSE205459_WT_and_KO_readcounts.txt"), header = T)
cts <- aggregate(cts[ , -1], list(gene_id = cts[, 1]), FUN = sum)
rownames(cts) <- cts$gene_id
cts$gene_id <- NULL

coldata <- read.csv(file.path(work_dir, "SraRunTable.csv"))
rownames(coldata) <- coldata$Name
coldata$Name <- NULL

cts <- cts[, rownames(coldata)]

dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design    = ~ condition
)

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  ggtitle("PCA of samples")

ggsave(filename = file.path(save_dir, "PCA.png"), plot = pca, 
       width = 8, height = 6, dpi=600)

## DE
de_dir <- file.path(save_dir, "Differential_Expression")
dir.create(de_dir)

res <- results(dds, contrast = c("condition", "KO", "WT"))
write.csv(res, file = file.path(de_dir, "DE_results.csv"), row.names = TRUE)

v1 <- EnhancedVolcano(res,
                      lab = rownames(res),
                      x = "log2FoldChange",
                      y = "padj",
                      pCutoff = 0.05,
                      FCcutoff = 1.0,
                      title = "Volcano Plot",
                      subtitle = "WDR6-WKO vs WT")
ggsave(filename = file.path(de_dir, "Volcano_plot.png"), 
       plot = v1, width = 8, height = 6, dpi=600)

res_clean <- res[!is.na(res$log2FoldChange) & !is.na(res$padj), ]
sig_genes_list <- rownames(res_clean[abs(res_clean$log2FoldChange) > 1 & 
                                       res_clean$padj < 0.05, ])

mat <- counts(dds, normalized=T)
mat <- mat[sig_genes_list, , drop=F]

mat.z <- t(scale(t(as.matrix(mat))))
Heatmap(mat.z,
        cluster_rows = T,
        cluster_columns = F,
        name = "Z-Transformed Counts",
        row_labels = rownames(mat.z),
        row_names_gp = gpar(fontsize = 8))

### Pathway Enrichment

pe_dir <- file.path(save_dir, "Pathway_enrichment")
dir.create(pe_dir)

sig_genes <- res_clean[abs(res_clean$log2FoldChange) > 0.35 & 
                         res_clean$padj < 0.05, ]

sig_genes$entrez_id <- mapIds(org.Mm.eg.db,
                              keys = rownames(sig_genes),
                              column = "ENTREZID",
                              keytype = "SYMBOL",
                              multiVals = "first")

run_go <- function(gene_list) {
  
  onts <- c("BP", "CC", "MF")
  
  for (ont in onts) {
    
    ont_dir <- paste0(pe_dir, "/", ont)
    
    if (!dir.exists(ont_dir)) {
      dir.create(ont_dir, recursive = TRUE)
    }
    
    #Total
    ego <- enrichGO(gene = gene_list$entrez_id,
                    OrgDb = org.Mm.eg.db,
                    keyType = "ENTREZID",
                    ont = ont,
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.2,
                    readable = TRUE)
    
    if (!is.null(ego) && nrow(ego@result) > 0) {
      
      d3 <- dotplot(ego, showCategory=10) + ggtitle(paste("Total GO Enrichment -", ont))
      ggsave(filename = file.path(ont_dir, paste0("Total_GO_", ont, "_Enrichment_Dotplot.png")),
             plot = d3, width = 8, height = 6, dpi=600)
      write.csv(as.data.frame(ego), 
                file = file.path(ont_dir, paste0("Total_GO_", ont, "_Enrichment.csv")))
    }
    
    genes_up <- gene_list$entrez_id[gene_list$log2FoldChange > 0]
    genes_down <- gene_list$entrez_id[gene_list$log2FoldChange < 0]
    #Up
    up <- enrichGO(gene = genes_up,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = ont,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)
    
    if (!is.null(up) && nrow(up@result) > 0) {
      d3 <- dotplot(up, showCategory=10) + 
        ggtitle(paste("Upregulated GO Enrichment -", ont))
      ggsave(filename = file.path(ont_dir, paste0("Upregulated_GO_", ont, 
                                                  "_Enrichment_Dotplot.png")), 
             plot = d3, width = 8, height = 6, dpi=600)
      write.csv(as.data.frame(up), 
                file = file.path(ont_dir, paste0("Upregulated_GO_", 
                                                 ont, "_Enrichment.csv")))
    }
    
    #Down
    down <- enrichGO(gene = genes_down,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = ont,
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2,
                     readable = TRUE)
    
    if (!is.null(down) && nrow(down@result) > 0) {
      d3 <- dotplot(down, showCategory=10) + 
        ggtitle(paste("Downregulated GO Enrichment -", ont))
      ggsave(filename = file.path(ont_dir, paste0("Downregulated_GO_", 
                                                  ont, "_Enrichment_Dotplot.png")), 
             plot = d3, width = 8, height = 6, dpi = 600)
      write.csv(as.data.frame(down),
                file = file.path(ont_dir, paste0("Downregulated_GO_", 
                                                 ont, "_Enrichment.csv")))
    }
  }
}


run_kegg <- function(gene_list) {
  
  ont_dir <- file.path(pe_dir, "KEGG")
  
  if (!dir.exists(ont_dir)) {
    dir.create(ont_dir, recursive = TRUE)
  }
  
  #Total
  kegg <- enrichKEGG(gene = gene_list$entrez_id,
                     organism = 'mmu',            
                     keyType = "kegg",           
                     pAdjustMethod= "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.2)
  
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    d3 <- dotplot(kegg, showCategory = 10) + ggtitle("Total KEGG Enrichment")
    ggsave(filename = file.path(ont_dir, "Total_KEGG_Enrichment_Dotplot.png"),
           plot = d3, width = 8, height = 6, dpi = 600)
    write.csv(as.data.frame(kegg),
              file = file.path(ont_dir, "Total_KEGG_Enrichment.csv"))
  }
  
  genes_up <- gene_list$entrez_id[gene_list$log2FoldChange > 0]
  genes_down <- gene_list$entrez_id[gene_list$log2FoldChange < 0]
  
  #Up
  kegg_up <- enrichKEGG(gene = genes_up,
                        organism = 'mmu',
                        keyType = "kegg",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.2)
  if (!is.null(kegg_up) && nrow(kegg_up@result) > 0) {
    d3 <- dotplot(kegg_up, showCategory = 10) + 
      ggtitle("Upregulated KEGG Enrichment")
    ggsave(filename = file.path(ont_dir, "Upregulated_KEGG_Enrichment_Dotplot.png"),
           plot = d3, width = 8, height = 6, dpi = 600)
    write.csv(as.data.frame(kegg_up),
              file = file.path(ont_dir, "Upregulated_KEGG_Enrichment.csv"))
  }
  
  
  #Down
  kegg_down <- enrichKEGG(gene = genes_down,
                          organism = 'mmu',
                          keyType = "kegg",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)
  if (!is.null(kegg_down) && nrow(kegg_down@result) > 0) {
    d3 <- dotplot(kegg_down, showCategory = 10) + 
      ggtitle("Downregulated KEGG Enrichment")
    ggsave(filename = file.path(ont_dir, "Downregulated_KEGG_Enrichment_Dotplot.png"),
           plot = d3, width = 8, height = 6, dpi = 600)
    write.csv(as.data.frame(kegg_down),
              file = file.path(ont_dir, "Downregulated_KEGG_Enrichment.csv"))
  }
  
}

run_go(sig_genes)
run_kegg(sig_genes)
