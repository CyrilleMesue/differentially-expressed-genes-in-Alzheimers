#!/usr/bin/env python

import pandas as pd
import numpy as np

# Load raw count and RPKM data
raw_count_df = pd.read_csv('data/raw_count_data.tsv', sep='\t')
rpkm_count_df = pd.read_csv('data/rpkm_count_data.tsv', sep='\t')

# Remove duplicates based on GeneSymbol column
raw_count_df = raw_count_df.drop_duplicates(keep='first')
rpkm_count_df = rpkm_count_df.drop_duplicates(keep='first')
raw_count_df = raw_count_df.drop_duplicates(subset='GeneSymbol', keep='first')
rpkm_count_df = rpkm_count_df.drop_duplicates(subset='GeneSymbol', keep='first')

# Remove missing values
raw_count_df = raw_count_df.dropna().reset_index(drop = True)
rpkm_count_df = rpkm_count_df.dropna().reset_index(drop = True)


# Remove genes with zero values across all samples
raw_count_df = raw_count_df.loc[:, (raw_count_df != 0).any(axis=0)]
rpkm_count_df = rpkm_count_df.loc[:, (rpkm_count_df != 0).any(axis=0)]

# Remove genes with negative values
raw_count_df = raw_count_df.loc[(raw_count_df.iloc[:,1:] >= 0).all(axis=1), :]
rpkm_count_df = rpkm_count_df.loc[(rpkm_count_df.iloc[:,1:] >= 0).all(axis=1), :]

# Remove genes with very low summation
threshold_raw = 5
threshold_rpkm = 1

raw_sum = raw_count_df.iloc[:, 1:].sum(axis=1)
rpkm_sum = rpkm_count_df.iloc[:, 1:].sum(axis=1)

raw_count_df = raw_count_df[(raw_sum > threshold_raw)]
rpkm_count_df = rpkm_count_df[(rpkm_sum > threshold_rpkm)]

# Save preprocessed data to new TSV files
raw_count_df.to_csv('data/preprocessed_raw_count_data.tsv', sep='\t', index=False)
rpkm_count_df.to_csv('data/preprocessed_rpkm_count_data.tsv', sep='\t', index=False)

# Display the number of remaining genes
print(f"Number of genes after preprocessing (raw count): {len(raw_count_df)}")
print(f"Number of genes after preprocessing (RPKM): {len(rpkm_count_df)}")
print("Done!!!")
exit()


