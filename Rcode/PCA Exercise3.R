#Removing all the elements in the workspace
rm(list = ls())

#Installs the pasilla package (you only need to do this once).
BiocManager::install("pasilla")

#Loads the package into your session.
library (pasilla)

# For prcomp() function
library(stats)

# For visualization
library(ggplot2)    

#Finds the path to the pasilla_gene_counts.tsv file that comes with the package.
pasCts <- system.file ("extdata", "pasilla_gene_counts.tsv", package ="pasilla", mustWork = TRUE )

#Loads the gene count table into R, with genes as row names.
cts <- as.matrix ( read.csv(pasCts, sep="\t", row.names="gene_id"))

#Transpose of the matrix
cts<- t(cts)

#Removing columns with all zero elements 
cts <- cts[, colSums(cts != 0) > 0]

#Standardize the Data
standardized_data <- scale(cts)

#Perform PCA
#Apply PCA to the standardized data using the prcomp() function:

pca_result <- prcomp(standardized_data)

#Exploring PCA Results
#Let's explore the results of PCA:

#To see the proportion of variance explained by each principal component:

summary(pca_result)

#This output provides the standard deviations of each principal component and the proportion of total variance explained.

#To access the principal components themselves:

pca_components <- pca_result$rotation

#Visualize PCA Results

#Visualize the results of PCA using a scree plot to understand the amount of variance explained by each principal component:

scree_data <- data.frame(
  PC = 1:length(pca_result$sdev),
  VarianceExplained = pca_result$sdev^2 / sum(pca_result$sdev^2)
)

ggplot(scree_data, aes(x = PC, y = VarianceExplained)) +
  geom_bar(stat = "identity") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance Explained")
