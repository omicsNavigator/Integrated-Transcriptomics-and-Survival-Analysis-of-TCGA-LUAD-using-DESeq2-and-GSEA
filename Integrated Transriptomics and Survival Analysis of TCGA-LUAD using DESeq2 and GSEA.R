#Data acquisition
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
query <- GDCquery(project = "TCGA-LUAD",data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",workflow.type = "STAR - Counts")

GDCdownload(query)
data <- GDCprepare(query)

#Extract expression matrix
counts <- assay(data)

#Define groups
colnames(colData(data)) #just seeing inside the coldata if sample_type is available
metadata <- colData(data)

group <- ifelse(metadata$sample_type == "Primary Tumor","Tumor","Normal")
group <- factor(group)

table(group)

#Create DESeq2 dataset
colData_df <- data.frame(condition = group)
rownames(colData_df) <- colnames(counts)
dds <- DESeqDataSetFromMatrix(countData = counts,colData = colData_df,
                              design = ~ condition)
dds <- dds[rowSums(counts(dds))>10, ]

#Differential expression
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
summary(res)
head(res)

#Clean data by removing missing values
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$pvalue), ]

#Significant genes
sigGenes <- res[!is.na(res$padj) & res$padj < 0.05, ]
nrow(sigGenes)

#Volcano plot
library(ggplot2)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                             "Significant","Not Significant")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj),color = significant)) +
  geom_point(alpha = 0.5) +theme_minimal() +
  labs(title = "TCGA Lung Cancer Volcano Plot", x= "Log2 Fold Change",y="-log10(p-value)")

#Pathway Enrichment(GO + KEGG)
#install required packages
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
#BiocManager::install("org.Hs.eg.db", lib = .libPaths()[1], force = TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidydr)

#Prepare gene list & Convert ENSEMBL TO ENTREZ IDs
gene_list <- rownames(sigGenes)
head(rownames(sigGenes))
gene_list_clean <- gsub("\\..*", "", gene_list) 
head(gene_list_clean)

gene_df <- bitr(gene_list_clean,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
gene_ids <- gene_df$ENTREZID
head(gene_df)

#GO Enrichment Analysis
go_result <- enrichGO(gene=gene_ids,OrgDb=org.Hs.eg.db,keyType="ENTREZID",
                      ont="BP",pAdjustMethod="BH",qvalueCutoff=0.05)
head(go_result)

#Plot GO result
dotplot(go_result, showCategory=10)

#KEGG Pathway Analysis
kegg_result <- enrichKEGG(gene=gene_ids,organism="hsa",pvalueCutoff=0.05)
head(kegg_result)

#KEGG visualization
dotplot(kegg_result, showCategory=10)

#Survival Analysis
#Get clinical data
clinical <- colData(data)

#Build survival time & status
time <- ifelse(is.na(clinical$days_to_death),clinical$days_to_last_follow_up,
               clinical$days_to_death)

status <- ifelse(clinical$vital_status == "Dead", 1, 0)

#Extract normalized expression
head(rownames(counts))
rownames(counts) <- gsub("\\..*","", rownames(counts))
kras_id <- bitr("KRAS",fromType = "SYMBOL",toType = "ENSEMBL",OrgDb = org.Hs.eg.db)$ENSEMBL[1]
expr <- counts(dds, normalized = TRUE)[kras_id, ]
head(rownames(counts))

#Get survival data

surv_data <- data.frame(expr=as.numeric(expr),time=time,status=status)

#Clean data
surv_data <- na.omit(surv_data)

#Split groups (Stratify high vs low expression)
surv_data$group <- ifelse(surv_data$expr > median(surv_data$expr), "High","Low")

#Kaplan-Meier model
#install.packages(c("survival","survminer"))
library(survival)
library(survminer)
fit <- survfit(Surv(time, status)~group,data = surv_data)

#Plot survival curve
ggsurvplot(fit,data = surv_data,pval = TRUE,risk.table = TRUE,title = "KRAS Survival Analysis (TCGA-LUAD")

#Lets study all the genes rather than just KRAS
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$stat), ]

res_df$ENSEMBL <- rownames(res_df)
res_df$ENSEMBL <- gsub("\\..*", "", res_df$ENSEMBL)

gene_map <- bitr(res_df$ENSEMBL,
                 fromType="ENSEMBL",
                 toType="ENTREZID",
                 OrgDb=org.Hs.eg.db)

res_df <- res_df[res_df$ENSEMBL %in% gene_map$ENSEMBL, ]
gene_map <- gene_map[match(res_df$ENSEMBL, gene_map$ENSEMBL), ]

res_df$ENTREZID <- gene_map$ENTREZID

res_df <- res_df[order(res_df$stat, decreasing=TRUE), ]
res_df <- res_df[!duplicated(res_df$ENTREZID), ]

gene_list <- res_df$stat
names(gene_list) <- res_df$ENTREZID
gene_list <- sort(gene_list, decreasing=TRUE)

#GSEA-GO Analysis
gsea_go <- gseGO(geneList = gene_list,OrgDb = org.Hs.eg.db,ont = "BP",
                keyType = "ENTREZID",minGSSize = 10,maxGSSize = 500,pvalueCutoff = 0.05,
                verbose = FALSE)

#GSEA-KEGG pathways
gsea_kegg <- gseKEGG(geneList = gene_list,organism = "hsa",minGSSize = 10,
                      pvalueCutoff = 0.05,verbose = FALSE)

#Dotplots
dotplot(gsea_go, showCategory = 10, title = "GSEA GO-BP")

dotplot(gsea_kegg, showCategory = 10, title = "GSEA KEGG")

#Enrichment plot
gseaplot2(gsea_kegg,geneSetID = 1,title = gsea_kegg$Description[1])

#Ridgeplot
ridgeplot(gsea_go, showCategory = 10)

length(gene_list)