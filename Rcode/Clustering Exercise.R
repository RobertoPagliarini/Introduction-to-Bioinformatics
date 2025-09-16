#Ensures reproducibility (so random results, like k-means, will be the same each run).
set.seed(2)

#Installs the pasilla package (you only need to do this once).
BiocManager::install("pasilla")

#Loads the package into your session.
library (pasilla)

#Finds the path to the pasilla_gene_counts.tsv file that comes with the package.
pasCts <- system.file ("extdata", "pasilla_gene_counts.tsv", package ="pasilla", mustWork = TRUE )

#Loads the gene count table into R, with genes as row names.
cts <- as.matrix ( read.csv(pasCts, sep="\t", row.names="gene_id"))

#What is the model organism being studied? How many genes and how many samples per condition were assayed
packageDescription("pasilla")

#Print some informations
print(paste("The Model organism being studied is Drosophila_melanogaster from package info."))

print(paste("Number of genes:", nrow(cts)))

print(paste("Number of samples:", ncol(cts)))

print(paste("Sample names are :", paste(colnames(cts), collapse=", ")))

# Log-transform to stabilize variance
log_cts <- log2(cts + 1)

# Run k-means on genes
km_genes <- kmeans(log_cts, centers = 3, nstart = 25)

# Show how many genes per cluster
table(km_genes$cluster)

# Show which genes belongs to which cluster
km_genes$cluster

# Order genes by their assigned cluster
ordered_genes <- rownames(log_cts)[order(km_genes$cluster)]

# Build annotation for clusters
annotation_row <- data.frame(Cluster = factor(km_genes$cluster[ordered_genes]))
rownames(annotation_row) <- ordered_genes

# Heatmap with genes grouped by cluster
install.packages("pheatmap")

library(pheatmap)

pheatmap(log_cts[ordered_genes, ],
         scale = "row",
         cluster_rows = FALSE,          # keep k-means order
         clustering_method = "complete",
         annotation_row = annotation_row,
         show_rownames = FALSE,
         main = "Genes grouped by k-means clusters (pasilla)")