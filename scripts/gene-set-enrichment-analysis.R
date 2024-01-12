# load required packages
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(pathview)


# set working directory
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

#############################
# Load and Prepare Input Data
#############################

# reading in data from deseq2
df = read.csv("output_files/DGE/Deseq2-results-all.tsv", header=TRUE, sep = "\t")
# df <- df[0:100,]

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$GeneName

# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


#############################
# Gene Set Enrichment
#############################
# perform gene set enrichment analysis
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ALIAS",
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")



###########################################
# Gene Set Enrichment Results Visualization
###########################################
# import the dose package
require(DOSE)

# Dot plot of enriched terms
dotplot(gse, showCategory=30, split=".sign", font.size =7, label_format = 50)

# Dot Plot of activated & Suppressed GO terms
dotplot(gse, showCategory=30, split=".sign", font.size =7, label_format = 50) + facet_grid(.~.sign)

# Gene-Concept Network plot of enriched terms
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, font.size =7)

# Gene-Concept Network plot of enriched terms Circular
cnetplot(gse, foldChange=gene_list, circular = TRUE, colorEdge = TRUE, font.size =5) +opts(legend.position="bottom")

# Heatmap plot of enriched terms. default (A), foldChange=geneList (B) 
p1 <- heatplot(gse, showCategory=5)
p2 <- heatplot(gse, foldChange=gene_list, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

# Enrichment Map
x2 <- pairwise_termsim(gse)
emapplot(x2, showCategory = 20)

# Ridgeplot for gene set enrichment analysis.
ridgeplot(gse, label_format = 50, showCategory =20) + labs(x = "enrichment distribution")

# gseaplot for GSEA result(by = "runningScore"). by = "runningScore" (A), by = "preranked" (B), default (C) 
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)

# Pmcplot of enrichment analysis (pubmed trend of enriched terms)
terms <- gse$Description[1:3]
terms <- gsub("(.{30})", "\\1\n", terms)  
pmcplot(terms, 2010:2023, proportion=FALSE)

################################################################################
# KEGG Gene Set Enrichment Analysis
################################################################################
# Convert gene IDs for gseKEGG function. We will lose some genes here because not all IDs will be converted
ids<-bitr(names(original_gene_list), fromType = "ALIAS", toType = "ENTREZID", OrgDb=organism)
# remove duplicate IDS (here I use "ALIAS", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ALIAS")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$GeneName %in% dedup_ids$ALIAS,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

# create a gseKEGG object
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

# Dot plot of enriched terms kegg
dotplot(kk2, showCategory=20, split=".sign", font.size =7, label_format = 50)

# Dot Plot of activated & Suppressed GO terms kegg
dotplot(kk2, showCategory=20, split=".sign", font.size =7, label_format = 50) + facet_grid(.~.sign)

# Gene-Concept Network plot of enriched terms kegg
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list, font.size =7)

# Gene-Concept Network plot of enriched terms Circular kegg
cnetplot(kk2, foldChange=gene_list, circular = TRUE, colorEdge = TRUE, font.size =5, showCategory=3)

# Heatmap plot of enriched terms. default (A), foldChange=geneList (B) kegg
p1 <- heatplot(kk2, showCategory=3)
p2 <- heatplot(kk2, foldChange=gene_list, showCategory=3)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

# Enrichment Map kegg
x2 <- pairwise_termsim(kk2)
emapplot(x2, showCategory = 20)

# Ridgeplot for gene set enrichment analysis. kegg
ridgeplot(kk2, label_format = 50, showCategory =20) + labs(x = "enrichment distribution")

# gseaplot for GSEA result(by = "runningScore"). by = "runningScore" (A), by = "preranked" (B), default (C) kegg
gseaplot(kk2, by = "all", title = gse$Description[1], geneSetID = 1)

# Pmcplot of enrichment analysis (pubmed trend of enriched terms)
terms <- kk2$Description[1:3]
terms <- gsub("(.{25})", "\\1\n", terms)  
pmcplot(terms, 2010:2023, proportion=FALSE)

# Produce the native KEGG plot (PNG)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05168", species = kegg_organism)

# Produce a different plot (PDF) (not displayed here)
dme <- pathview(gene.data=kegg_gene_list, pathway.id="hsa05168", species = kegg_organism, kegg.native = F)

