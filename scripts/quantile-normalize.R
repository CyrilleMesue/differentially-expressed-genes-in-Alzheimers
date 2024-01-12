# Install required packages if not already installed
BiocManager::install("preprocessCore")

# Load required libraries
library(preprocessCore)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(tibble)

# Set working directory
setwd(dirname(dirname(rstudioapi::getActiveDocumentContext()$path)))

# Load data and metadata
data <- read.csv("data/preprocessed_rpkm_count_data.tsv", header = TRUE, row.names = 1, sep = "\t")
metadata <- read.csv("data/metadata.tsv", header = TRUE, row.names = 1, sep = "\t")

# Example data generation (replace this with your actual data)
set.seed(123)
rpkm_data <- as.matrix(data)

# Log-scale transformation
log_rpkm_data <- log2(rpkm_data + 1)  # Adding 1 to avoid log(0)

# Quantile normalization
qn_rpkm_data <- normalize.quantiles(log_rpkm_data)

# Arrange data for plotting
df_before <- data.frame(Samples = melt(log_rpkm_data)[["Var2"]],
                        Value = melt(log_rpkm_data)[["value"]])

row.names(qn_rpkm_data) <- row.names(data)
colnames(qn_rpkm_data) <- colnames(data)
df_after <- data.frame(Samples = melt(qn_rpkm_data)[["Var2"]],
                       Value = melt(qn_rpkm_data)[["value"]])

# Boxplots
boxplot_before <- ggplot(df_before, aes(x = Samples, y = Value)) +
  geom_boxplot() +
  labs(title = "Boxplot Before Normalization",
       y = "Expression Value") +
  theme_gray() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

boxplot_after <- ggplot(df_after, aes(x = Samples, y = Value)) +
  geom_boxplot() +
  labs(title = "Boxplot After Normalization",
       y = "Expression Value") +
  theme_gray() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Q-Q plots
qqplot_before <- ggplot(df_before, aes(sample = Value)) +
  stat_qq() +
  stat_qq_line() + 
  labs(title = "Q-Q Plot Before Normalization",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_gray() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

qqplot_after <- ggplot(df_after, aes(sample = Value)) +
  stat_qq() +
  stat_qq_line() +
  labs(title = "Q-Q Plot After Normalization",
       x = "Theoretical Quantiles",
       y = "Sample Quantiles") +
  theme_gray() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Arrange plots in a 2x2 grid
combined_plot <- grid.arrange(
  boxplot_before + theme(plot.title = element_text(hjust = 0.5),
                         axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)),
  boxplot_after + theme(plot.title = element_text(hjust = 0.5),
                        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)), 
  qqplot_before + theme(plot.title = element_text(hjust = 0.5)), 
  qqplot_after + theme(plot.title = element_text(hjust = 0.5)), ncol = 2)

# Save the combined plot
ggsave("plots/combined_normalization_plot.png", plot = combined_plot, width = 12, height = 8)

# Save normalized data
gene_symbols <- data.frame(GeneSymbol = rownames(qn_rpkm_data), stringsAsFactors = FALSE)
qn_rpkm_data <- cbind(gene_symbols, qn_rpkm_data)

write.table(qn_rpkm_data, "data/normalized_rpkm_count_data.tsv", sep = "\t", row.names = FALSE)
