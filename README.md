# RNA-Seq Data Analysis Pipeline

## Overview

This repository contains the code and output files for a comprehensive RNA-Seq data analysis pipeline. The pipeline includes steps for preprocessing raw count data, conducting differential gene expression analysis using DESeq2, and performing gene set enrichment analysis (GSEA) on the obtained results. The analysis focuses on identifying differentially expressed genes and enriched biological pathways.

## Project Structure

The project is organized into the following directories:

- **data:** Contains the raw and preprocessed RNA-Seq data files, as well as metadata information.
  
- **output_files:** Subdivided into several subdirectories, each containing specific outputs from the analysis.

  - **DGE:** Differential gene expression analysis results and associated plots.
  
  - **GSEA:** Gene set enrichment analysis results and visualization plots.
  
  - **KEGG_GSEA:** KEGG pathway enrichment analysis results and associated plots.
  
  - **plots:** Various plots generated during the analysis, including PCA plots and normalization plots.

- **scripts:** Houses the R and Python scripts used for data preprocessing, analysis, and visualization.

## Quick Start

To reproduce the analysis, follow these steps:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/CyrilleMesue/differentially-expressed-genes-in-Alzheimers.git
   cd differentially-expressed-genes-in-Alzheimers
   ```

2. **Install dependencies:**
   Ensure you have R and the required packages installed. You can find a list of dependencies R in the r_requirements.txt file. Dependencies can be installed as follows.

   _**For Linux Installation**_
     ```
      sudo apt install libxml2-dev
      sudo apt install libcurl4-openssl-dev

     ```

   _**Python Dependencies**_
      ```bash
      pip install -r requirements.txt
      ```

   _**R Dependencies**_
      ```R
      if (!require("BiocManager", quietly = TRUE))
          install.packages("BiocManager")

      BiocManager::install("DESeq2")
      BiocManager::install("clusterProfiler")
      BiocManager::install('cowplot')
      install.packages("pheatmap")
      install.packages("RColorBrewer")
      install.packages("ggrepel")
      install.packages("preprocessCore")
      install.packages("enrichplot")
      install.packages("pathview")
      install.packages("ggfortify")
      install.packages("gridExtra")
      install.packages("tidyverse")
      install.packages("shinydashboard")
      ```

3. **Run the analysis:**
   Execute the R scripts in the `scripts` directory in the following order:
   - `preprocess-data.py`
   - `differential-gene-expression-analysis.R`
   - `gene-set-enrichment-analysis.R`

4. **Explore results:**
   The results will be generated in the `output_files` directory, organized into relevant subdirectories.

## Results Summary

Here is a brief summary of key results:

- **Differential Gene Expression Analysis:**
  - Volcano plots, MA plots, Heatmaps, PCA plots, and Cook's distance plots.
  - Differentially expressed genes identified in `Deseq2-results-all.tsv`.

- **Gene Set Enrichment Analysis (GSEA):**
  - Enrichment maps, dot plots, and cnet plots.
  - Results in the `GSEA` directory.

- **KEGG Pathway Enrichment Analysis:**
  - Enrichment maps, dot plots, and cnet plots specific to KEGG pathways.
  - Results in the `KEGG_GSEA` directory.

## Additional Information

For additional details, refer to the individual scripts in the `scripts` directory. Feel free to explore the code, modify parameters, and adapt the pipeline to your specific needs.


## Contributors 

<table>
  <tr>
    <td align="center"><a href="https://github.com/CyrilleMesue"><img src="https://avatars.githubusercontent.com/CyrilleMesue" width="100px;" alt=""/><br /><sub><b>Cyrille M. NJUME</b></sub></a><br /></td>
  </tr>
</table>

## References 

[1] Love, M. I. (2023, November 14). Analyzing RNA-seq data with DESeq2. Bioconductor.org. https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

[2] Gene Set Enrichment Analysis with ClusterProfiler. (2019, June 14). NGS Analysis. https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

[3] Yu, G. (2015). Chapter 15 Visualization of functional enrichment result | Biomedical Knowledge Mining using GOSemSim and clusterProfiler. Yulab-Smu.top. https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html


### Contact

For any feedback or queries, please reach out to [cyrillemesue@gmail.com](mailto:cyrillemesue@gmail.com).
