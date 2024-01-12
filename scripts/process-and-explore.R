# load required packages
library(ggfortify)
library(ggplot2)
library(gridExtra)

# set working directory
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

# Define pca and visualization function
perform_pca_visualization <- function(data, metadata, samplenames, groups, version, outdir) {
#' Perform PCA and visualize the results.
#'
#' This function performs PCA on the provided data and then visualizes it using a scatter plot to examine 
#' variation or any outliers and clusters.
#'
#' Inputs:
#'   \code{data}: \code{data.frame} object containing normalized data.
#'   \code{metadata}: \code{data.frame} object containing summary information about our data.
#'   The number of rows in metadata must match the number of columns in data.
#'   \code{sample_names}: The column name in metadata that contains the unique names of samples.
#'   \code{groups}: The column name from metadata that contains the group or class or type to which each sample belongs.
#'   \code{version}: Version number to append to the plot name.
#'   \code{outdir}: Output directory.
#'
#' Output:
#'   Saved PCA plots; with and without labels.
#'
#' @param data data.frame object containing normalized data.
#' @param metadata data.frame object containing summary information about our data.
#' @param sample_names The column name in metadata that contains the unique names of samples.
#' @param groups The column name from metadata that contains the group or class or type to which each sample belongs.
#' @param version Version number to append to the plot name.
#' @param outdir Output directory.
#'
#' @return None
#'
#' @examples
#' \dontrun{
#' data <- read.csv("data/normalized_rpkm_count_data.tsv", header = TRUE, row.names = 1, sep = "\t")
#' metadata <- read.csv("data/metadata.tsv", header = TRUE, row.names = 1, sep = "\t")
#' perform_pca_visualization(data, metadata, sample_names = "SampleNames", groups = "Group", version = 1, outdir = "plots")
  # Perform PCA
  pca_result <- prcomp(data, scale. = TRUE, center = TRUE)
  
  # Extract PC scores
  pca_scores <- data.frame(PC1 = pca_result$x[, 1],
                           PC2 = pca_result$x[, 2],
                           Group = metadata[[groups]],
                           SampleNames = metadata[[samplenames]]
  )
  
  # Extract percentage variation from the PCA summary
  pc1_percentage <- round(summary(pca_result)$importance[2, 1]*100, 2)
  pc2_percentage <- round(summary(pca_result)$importance[2, 2]*100, 2)
  
  # Visualize without sample names
  plot1 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2) +
    labs(title = "PCA Scatter Plot",
         x = sprintf("PC1 (%.2f%%)", pc1_percentage),
         y = sprintf("PC2 (%.2f%%)", pc2_percentage)) +
    theme_gray() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Save the plot
  ggsave(filename = sprintf("%s/PCA-Scatter-Plot-%d.png", outdir, version), plot = plot1, device = "png", width = 8, height = 6)
  
  # Visualize with sample names
  plot2 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = Group, label = SampleNames)) +
    geom_point(size = 1) +
    labs(title = "PCA Scatter Plot",
         x = sprintf("PC1 (%.2f%%)", pc1_percentage),
         y = sprintf("PC2 (%.2f%%)", pc2_percentage)) +
    geom_text(size = 3, vjust = 2) +
    theme_gray() +
    theme(plot.title = element_text(hjust = 0.5))
  # Save the plot
  ggsave(filename = sprintf("%s/PCA-Scatter-Plot-with-labels%d.png", outdir, version), plot = plot2 , device = "png", width = 8, height = 6)
 
  combined_plot <- grid.arrange(plot1, plot2, ncol = 2, widths = c(1.5,1.5), heights = c(1))

  # Save the combined plot
  ggsave(filename = sprintf("%s/combined-PCA-Scatter-Plot-%d.png", outdir, version), 
         plot = combined_plot, device = "png", width = 12, height = 6)
  
}

# define function to remove outliers
remove_outliers <- function(data, metadata, samplenames, outliernames) {
  #' Remove outliers from data and metadata
  #'
  #' This function removes specified outliers from both the data matrix and metadata.
  #'
  #' @param data A data.frame object containing the data matrix.
  #' @param metadata A data.frame object containing metadata information.
  #' @param samplenames The column name in metadata that contains the unique names of samples.
  #' @param outliernames A character vector containing sample names to be removed as outliers.
  #'
  #' @return A list containing the cleaned data and metadata.
  #'
  #' @examples
  #' \dontrun{
  #' data <- read.csv("data/normalized_rpkm_count_data.tsv", header = TRUE, row.names = 1, sep = "\t")
  #' metadata <- read.csv("data/metadata.tsv", header = TRUE, row.names = 1, sep = "\t")
  #' outliers <- c("Sample1", "Sample2", "Sample3")
  #' result <- remove_outliers(data, metadata, samplenames = "SampleNames", outliernames = outliers)
  #' cleaned_data <- result$data
  #' cleaned_metadata <- result$metadata
  #' 

    # Filter data based on sample names
    cleaned_data <- data[!(rownames(data) %in% outliernames), , drop = FALSE]
    
    # Filter metadata based on sample names
    cleaned_metadata <- subset(metadata, !(metadata[[samplenames]] %in% outliernames))
    
    # Return cleaned data and metadata as a list
    return(list(data = cleaned_data, metadata = cleaned_metadata))
}

# Load data and metadata
data <- read.csv("data/normalized_rpkm_count_data.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.csv("data/metadata.tsv", header = TRUE, row.names = 1, sep = "\t")

# Transpose the data
data <- t(data)

# perform pca and visualize
perform_pca_visualization(data=data, metadata=metadata, samplenames="SampleNamesRPKM", 
                          groups="Group", version=1, outdir="plots")

# remove outliers and observe again
results_outliers_removed <- remove_outliers(data, metadata, "SampleNamesRPKM", c("C5_rpkm", "A13_rpkm", "C2_rpkm"))
data2 <- results_outliers_removed$data
metadata2 <- results_outliers_removed$metadata
perform_pca_visualization(data=data2, metadata=metadata2, samplenames="SampleNamesRPKM", 
                          groups="Group", version=2, outdir="plots")

# remove outliers and observe again
results_outliers_removed <- remove_outliers(data2, metadata2, "SampleNamesRPKM", c("C6_rpkm"))
data3 <- results_outliers_removed$data
metadata3 <- results_outliers_removed$metadata
perform_pca_visualization(data=data3, metadata=metadata3, samplenames="SampleNamesRPKM", 
                          groups="Group", version=3, outdir="plots")

# save the raw versions of the final data with outliers removed
# Load data and metadata
data_raw <- read.csv("data/preprocessed_raw_count_data.tsv", header = TRUE, row.names = 1, sep = "\t")
data_raw <- t(data_raw)
results_outliers_removed_raw <- remove_outliers(data_raw, metadata, "SampleNamesRaw", c("C5_raw", "A13_raw", "C2_raw", "C6_raw"))
data_raw <- results_outliers_removed_raw$data
metadata_final <- results_outliers_removed_raw$metadata

# save trimmed data and metadata
data_raw <- t(data_raw)
gene_symbols <- data.frame(GeneSymbol = rownames(data_raw), stringsAsFactors = FALSE)
data_raw <- cbind(gene_symbols, data_raw)
write.table(data_raw, "data/final_preprocessed_raw_count_data.tsv", sep = "\t", row.names = FALSE)
write.table(metadata_final, "data/metadata_final.tsv", sep = "\t", row.names = TRUE)
