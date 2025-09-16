#Elbow plot → visual inspection of within-cluster sum of squares.

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

# Log-transform to stabilize variance
log_cts <- log2(cts + 1)

#To store total within-cluster sum of squares
wss <- numeric()  # to store total within-cluster sum of squares
k_values <- 2:10  # try k from 1 to 10

for (k in k_values) {
  set.seed(2)
  km <- kmeans(log_cts, centers = k, nstart = 25)
  wss[k-1] <- km$tot.withinss
}

# Plot Elbow curve
plot(k_values, wss, type = "b", pch = 19,
     xlab = "Number of clusters K",
     ylab = "Total within-cluster sum of squares",
     main = "Elbow method for selecting K")