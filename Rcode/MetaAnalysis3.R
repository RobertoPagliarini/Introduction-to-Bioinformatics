# --- Install / load packages -----------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma", "annotate", "hgu133a.db"), update = FALSE, ask = FALSE)
if (!requireNamespace("metaMA", quietly = TRUE)) install.packages("metaMA")

library(GEOquery)
library(limma)
library(metaMA)
library(annotate)
library(hgu133a.db)  # annotation package for GPL96

# --- Step 1: download GEO datasets -----------------------------------------
gse1 <- getGEO("GSE2034", GSEMatrix = TRUE)[[1]]
gse2 <- getGEO("GSE2990", GSEMatrix = TRUE)[[1]]

expr1 <- exprs(gse1)
expr2 <- exprs(gse2)

# --- Step 2: define phenotype classes (0/1) --------------------------------
# GSE2034: relapse column exists
group1 <- as.numeric(pData(gse1)[,"bone relapses (1=yes, 0=no):ch1"])
table(group1)  # sanity check

# GSE2990: define 0/1 automatically from 'grade:ch1'
grade_col <- as.numeric(pData(gse2)[,"grade:ch1"])
valid_idx <- grade_col %in% c(1,3)           # keep only grade 1 or 3
expr2 <- expr2[, valid_idx]                  # subset expression matrix
group2 <- ifelse(grade_col[valid_idx] == 1, 0, 1)
table(group2)  # sanity check

# --- Step 3: prepare 'esets' with common genes --------------------------------
common_genes <- intersect(rownames(expr1), rownames(expr2))
if(length(common_genes) == 0) stop("No common genes!")

esets <- list(
  expr1[common_genes, , drop=FALSE],
  expr2[common_genes, , drop=FALSE]
)

# --- Step 4: prepare class labels as a list --------------------------------
classes <- list(group1, group2)

# --- Step 5: run meta-analysis using p-value combination -------------------
res <- pvalcombination(esets = esets, classes = classes, moderated = "limma", BHth = 0.05)

# --- Step 6: map probes to gene symbols ------------------------------------
# Using hgu133a.db for GPL96
probe2symbol <- unlist(mget(common_genes, hgu133aSYMBOL, ifnotfound=NA))
DE_genes_symbols <- probe2symbol[res$Meta]
head(DE_genes_symbols)  # top DE genes

# Compute raw and adjusted p-values
raw_pval_meta <- 2 * (1 - pnorm(abs(res$TestStatistic)))
adj_pval_meta <- p.adjust(raw_pval_meta, method="BH")

# Show top 10 DE genes by adjusted p-value
top_idx <- order(adj_pval_meta)[1:10]
data.frame(
  gene = probe2symbol[top_idx],
  rawP = raw_pval_meta[top_idx],
  adjP = adj_pval_meta[top_idx]
)

# --- Step 7: visualization -------------------------------------------------
hist(raw_pval_meta, breaks = 100, main = "Histogram of raw meta p-values", xlab = "p-value")
