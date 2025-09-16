#Removing all the elements in the workspace
rm(list = ls())

#Step 1: Load the Required Packages

# For prcomp() function
library(stats)

# For visualization
library(ggplot2)    

#Step 2: Load and Explore the Dataset
#The "USArrests" dataset is built into R. Load it and explore the first few rows:
  
data("USArrests")

head(USArrests)

#Step 3: Standardize the Data
#PCA requires data standardization. Standardize the variables in the dataset:
  
standardized_data <- scale(USArrests)

#Step 4: Perform PCA
#Apply PCA to the standardized data using the prcomp() function:
  
pca_result <- prcomp(standardized_data)

#Step 5: Exploring PCA Results
#Let's explore the results of PCA:

#To see the proportion of variance explained by each principal component:

summary(pca_result)

#This output provides the standard deviations of each principal component and the proportion of total variance explained.

#To access the principal components themselves:

pca_components <- pca_result$rotation

#Step 6: Visualize PCA Results

#Visualize the results of PCA using a scree plot to understand the amount of variance explained by each principal component:
  
  scree_data <- data.frame(
    PC = 1:length(pca_result$sdev),
    VarianceExplained = pca_result$sdev^2 / sum(pca_result$sdev^2)
  )

ggplot(scree_data, aes(x = PC, y = VarianceExplained)) +
  geom_bar(stat = "identity") +
  labs(title = "Scree Plot", x = "Principal Component", y = "Proportion of Variance Explained")

#Step 7: Interpretation of Results

#The summary(pca_result) output helps you understand the proportion of variance explained by each principal component. 
#This informs you about which components capture the most variability in your data.

#The pca_components matrix contains the loading scores for each variable on each principal component. 
#Larger absolute values indicate stronger contributions to that component.

#The scree plot provides an overview of how much variance is captured by each principal component. 
#The "elbow point" suggests how many components to retain to explain a significant portion of the data's variability.



