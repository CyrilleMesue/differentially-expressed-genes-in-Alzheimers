#!/usr/bin/env python

import pandas as pd
import gzip

# Replace 'path/to/downloaded/GSE53697_RNAseq_AD.txt.gz' with the actual path to your downloaded file
file_path = 'data/GSE53697_RNAseq_AD.txt.gz'

# Read the gzipped file into a pandas DataFrame
with gzip.open(file_path, 'rt') as file:
    df = pd.read_csv(file, sep='\t')

# Separate columns into raw and rpkm
raw_columns = [col for col in df.columns if 'raw' in col]
rpkm_columns = [col for col in df.columns if 'rpkm' in col]

# Create dataframes with GeneID and GeneSymbol, and corresponding raw and rpkm columns
raw_count_df = df[['GeneSymbol'] + raw_columns]
rpkm_count_df = df[['GeneSymbol'] + rpkm_columns]


# Save dataframes to TSV files
raw_count_df.to_csv('data/raw_count_data.tsv', sep='\t', index=False)
rpkm_count_df.to_csv('data/rpkm_count_data.tsv', sep='\t', index=False)

print(f"Number of genes in raw count data: {len(raw_count_df)}")
print(f"Number of genes after RPKM data: {len(rpkm_count_df)}")
print("Done!!!")
exit()