# --- Install / load packages -----------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "limma"), update=FALSE, ask=FALSE)
if (!requireNamespace("metaMA", quietly = TRUE)) install.packages("metaMA")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(GEOquery)
library(limma)
library(metaMA)
library(ggplot2)

# --- Step 1: download GEO datasets -----------------------------------------
gse1 <- getGEO("GSE2034", GSEMatrix=TRUE)[[1]]
gse2 <- getGEO("GSE2990", GSEMatrix=TRUE)[[1]]

expr1 <- exprs(gse1)
expr2 <- exprs(gse2)

# --- Step 2: define phenotype classes (0/1) --------------------------------
# GSE2034: "bone relapses (1=yes, 0=no):ch1"
group1 <- as.numeric(pData(gse1)[, "bone relapses (1=yes, 0=no):ch1"])
if(length(unique(group1)) != 2) stop("GSE2034 does not have exactly 2 classes")

# GSE2990: choose "relapse:ch1" (or any column with exactly 2 classes)
# First, inspect candidate columns:
sapply(pData(gse2), function(x) length(unique(na.omit(x))))  # see which are binary

# Suppose we choose "distant rfs:ch1" and remove NAs
valid_idx <- which(pData(gse2)[, "distant rfs:ch1"] %in% c(0,1))
expr2 <- expr2[, valid_idx]
group2 <- as.numeric(pData(gse2)[valid_idx, "distant rfs:ch1"])

if(length(unique(group2)) != 2) stop("GSE2990 does not have exactly 2 classes after removing NAs")

# --- Step 3: keep common genes/probes --------------------------------------
common_genes <- intersect(rownames(expr1), rownames(expr2))
esets <- list(expr1[common_genes, ], expr2[common_genes, ])

# --- Step 4: prepare class list --------------------------------------------
classes <- list(group1, group2)

# --- Step 5: meta-analysis using p-value combination -----------------------
res <- pvalcombination(esets = esets, classes = classes, moderated="limma", BHth=0.05)

# --- Step 6: extract results ------------------------------------------------
DE_genes <- rownames(esets[[1]])[res$Meta]
length(DE_genes)
head(DE_genes)

# Compute raw and adjusted p-values
raw_pval_meta <- 2*(1 - pnorm(abs(res$TestStatistic)))
adj_pval_meta <- p.adjust(raw_pval_meta, method="BH")

# --- Step 7: volcano plot with top 10 genes labeled -------------------------
top_idx <- order(adj_pval_meta)[1:10]
volcano_data <- data.frame(
  gene = rownames(esets[[1]]),
  logFC = res$TestStatistic,
  adjP = adj_pval_meta
)

top_genes <- volcano_data$gene[top_idx]

ggplot(volcano_data, aes(x=logFC, y=-log10(adjP), label=gene)) +
  geom_point(alpha=0.6) +
  geom_text(data=subset(volcano_data, gene %in% top_genes), vjust=-1, size=3) +
  theme_minimal() +
  xlab("Meta-analysis statistic") +
  ylab("-log10(adj P-value)") +
  ggtitle("Volcano plot of meta-analysis") +
  geom_hline(yintercept=-log10(0.05), color="red", linetype="dashed")
