library(AnnotationDbi)
library(cluster)
library(clusterProfiler)
library(DESeq2)
library(fastDummies)
library(ggplot2)
library(GOstats)
library(igraph)
library(netZooR)
library(org.Hs.eg.db)
library(psych)
library(r2r)
library(recount3)
library(igraph)
#library(WGCNA)

# Extracting data from TCGA
data <- recount3::create_rse_manual(
  project = "THCA",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

# Gene expression pre-processing, for more information see 
# https://github.com/netZoo/netbooks/blob/main/netbooks/netZooR/gene_expression_for_coexpression_nets.ipynb
G <- transform_counts(data, by = "mapped_reads")
G <- G[data@rowRanges@elementMetadata@listData$gene_type == "protein_coding",]
G <- G[-which(rowSums(G) <= 1),] # Filtering: remove genes with no counts
countMat=SummarizedExperiment::assay(DESeqDataSetFromMatrix(G, data.frame(row.names=seq_len(ncol(G))), ~1), 1)
expression <- vst(countMat, blind=FALSE)

# Extract metadata that we want to include in COBRA
metadata_url <- locate_url(
  "THCA",
  "data_sources/tcga")
metadata <- read_metadata(file_retrieve(url = metadata_url))

cancer <- metadata$tcga.gdc_cases.samples.sample_type
cancer <- ifelse(cancer == "Solid Tissue Normal", 0, 1)

write.csv(expression[,cancer == 0], "THCA_nat.csv")
write.csv(expression[,cancer == 1], "THCA_case.csv")
