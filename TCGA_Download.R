#---- Load libraries------------------------------------------------------------
# BiocManager::install("TCGAbiolinks")
# BiocManager::install("SummarizedExperiment")
# BiocManager::install("biomaRt")
library(TCGAbiolinks)
library(SummarizedExperiment)
library(biomaRt)
library(dplyr)
#---- Set up parameters --------------------------------------------------------
project <- 'TCGA-BLCA'
dir_data <- './data'
if (!dir.exists(dir_data)) {
  dir.create(dir_data)
}
#---- Project summary ----------------------------------------------------------
tryCatch({
  getProjectSummary(project)
}, error = function(e) {
  cat("Error in fetching project summary:", e$message, "\n")
})
#---- Count matrix -------------------------------------------------------------
# Download the count matrix
Expr_counts <- GDCquery(project = project,
                        data.category = 'Transcriptome Profiling', 
                        data.type = 'Gene Expression Quantification',
                        workflow.type = "STAR - Counts" ) 
GDCdownload(Expr_counts, method = "api")
# Prepare data and convert to data frame
df_counts <- as.data.frame(assay(GDCprepare(query = Expr_counts)))
df_counts$gene_id <- row.names(df_counts)
# Convert gene ids to gene names using biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensg_ids <- gsub("\\..*", "", df_counts$gene_id)
gene_mapping <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = ensg_ids,
  mart = ensembl
)
gene_mapping <- gene_mapping[match(ensg_ids, gene_mapping$ensembl_gene_id), ]
df_counts$gene_name <- gene_mapping$hgnc_symbol
df_counts <- df_counts %>%
  select(tail(names(.), 2), everything())
# Save the count matrix
write.csv(df_counts, file.path(dir_data, paste0(project, '_counts.csv')), 
          row.names = FALSE)
#---- Clinical data ------------------------------------------------------------
# Download the clinical data
df_clinical <- GDCquery_clinic(project = project, type = "clinical")
# Save the clinical data
write.csv(df_clinical, file.path(dir_data, paste0(project, '_clinical.csv')), 
          row.names = FALSE)
#---- Biospecimen data ---------------------------------------------------------
# Download the Biospecimen data
df_biospecimen <- GDCquery_clinic(project = project, type = "biospecimen")
# Identify list columns and convert them to character
df_biospecimen[] <- lapply(df_biospecimen, function(x) {
  if (is.list(x)) {
    sapply(x, toString)
  } else {
    x
  }
})
# Save the Biospecimen data
write.csv(df_biospecimen, file.path(dir_data, paste0(project, '_biospecimen.csv')), 
          row.names = FALSE)
