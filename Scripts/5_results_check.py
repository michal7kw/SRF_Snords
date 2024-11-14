# Converted from 5_results_check.ipynb

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import mygene
import os
pd.set_option('display.max_columns', None)

# Set the working directory
working_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords"
os.chdir(working_dir)
print(f"Current working directory: {os.getcwd()}")

from functions import check_dexseq_results


# # Examind the data

# Example usage
results_file = "output_v38/dexseq_results_PW1_vs_combined_controls.csv"
df = pd.read_csv(results_file)

df.head()

hist = plt.hist(df['log2fold_treated_control'], bins=50)


summary_stats = check_dexseq_results(results_file)

# Print summary
print("\n=== Analysis Summary ===")
for key, value in summary_stats.items():
    print(f"{key}: {value}")

# # Select significant results

df_sig = df[df['padj'] < 0.05]

print(df.shape)
print(df_sig.shape)

df = df_sig

# # Remove overlapping genes

df[df['groupID'].str.contains('\+')].head()

df = df[~df['groupID'].str.contains('\+')]

df.head()

# # Check problematic genes

# List of problematic genes from the error message
# problematic_genes = ['ENSG00000285404.1', 'ENSG00000100150.19', 
#                     'ENSG00000128245.15', 'ENSG00000252909.1']

# Create masks for each condition
extreme_fc_mask = (df['log2fold_treated_control'].abs() > 5)
high_disp_mask = (df['dispersion'] > 10)
missing_vals_mask = (
    df['log2fold_treated_control'].isna() |
    df['pvalue'].isna() |
    df['padj'].isna()
)
extreme_stat_mask = (df['stat'].abs() > 10000)
# problematic_genes_mask = df['groupID'].str.contains('|'.join(problematic_genes), regex=True)

# Print counts for each condition
print(f"Records with extreme fold changes (>5): {extreme_fc_mask.sum()}")
print(f"Records with high dispersion (>10): {high_disp_mask.sum()}")
print(f"Records with missing values: {missing_vals_mask.sum()}")
print(f"Records with extreme test statistics (>10000): {extreme_stat_mask.sum()}")
# print(f"Records from problematic genes: {problematic_genes_mask.sum()}")

df.columns

# # Examine extreme test statistics

# Print genes with extreme test statistics
# print(list(df[extreme_stat_mask]['groupID'].unique()))
df[extreme_stat_mask].head()

df.sort_values('stat', ascending=False)[:10]

n_sel = 100 
stat_sel = list(df.sort_values('stat', ascending=False)[:n_sel]['groupID'][~df.sort_values('stat', ascending=False)[:n_sel]['groupID'].str.contains('\+')].unique())
print(stat_sel)


# Initialize mygene client
mg = mygene.MyGeneInfo()

# Remove version numbers from ENSEMBL IDs
stat_sel_no_version = [id.split('.')[0] for id in stat_sel]

# Query the gene symbols
# results = mg.querymany(stat_sel_no_version, scopes='ensembl.gene', fields='symbol', species='human')
results = mg.querymany(stat_sel_no_version, 
                      scopes='ensembl.gene', 
                      fields='symbol', 
                      species='human',
                      assembly='GRCh38')

# Create a dictionary mapping ENSEMBL IDs to symbols
gene_map = {res['query']: res.get('symbol', 'Not found') for res in results}

gene_symbols = [gene_map[ensembl_id.split('.')[0]] for ensembl_id in stat_sel]

ensembl_to_symbol = dict(zip(stat_sel, gene_symbols))

# Print results
for ensembl, symbol in ensembl_to_symbol.items():
    print(f"{ensembl}: {symbol}")

# Examine records with extreme test statistics
extreme_stat_records = df[extreme_stat_mask][['groupID', 'featureID', 'stat', 'pvalue', 'padj', 'log2fold_treated_control', 'dispersion']]
print("\nSample of records with extreme test statistics:")
extreme_stat_records.head()


# # Examine extreme fold changes

# Print genes with extreme fold changes
# print(list(df[extreme_fc_mask]['groupID'].unique()))
df[extreme_fc_mask].head()

df.sort_values('log2fold_treated_control', key=abs, ascending=False)[:10]


n_sel = 100 
fc_sel = list(df.sort_values('log2fold_treated_control', key=abs, ascending=False)[:n_sel]['groupID'][~df.sort_values('log2fold_treated_control', key=abs, ascending=False)[:n_sel]['groupID'].str.contains('\+')].unique())
print(fc_sel)



# Initialize mygene client
mg = mygene.MyGeneInfo()

# Remove version numbers from ENSEMBL IDs
fc_sel_no_version = [id.split('.')[0] for id in fc_sel]

# Query the gene symbols
# results = mg.querymany(fc_sel_no_version, scopes='ensembl.gene', fields='symbol', species='human')
results = mg.querymany(fc_sel_no_version, 
                      scopes='ensembl.gene', 
                      fields='symbol', 
                      species='human',
                      assembly='GRCh38')

# Create a dictionary mapping ENSEMBL IDs to symbols
gene_map = {res['query']: res.get('symbol', 'Not found') for res in results}

gene_symbols = [gene_map[ensembl_id.split('.')[0]] for ensembl_id in fc_sel]

ensembl_to_symbol = dict(zip(fc_sel, gene_symbols))

# Print results
for ensembl, symbol in ensembl_to_symbol.items():
    print(f"{ensembl}: {symbol}")

# # Save filtered results

clean_df = df[
    (df['log2fold_treated_control'].notna()) &
    (df['pvalue'].notna()) &
    (df['padj'].notna())
]

# Print summary of filtering
print(f"Original number of records: {len(df)}")
print(f"Number of records after cleaning: {len(clean_df)}")
print(f"Number of records removed: {len(df) - len(clean_df)}")

# Save cleaned results to new CSV file
output_file = "output/dexseq_results_PW1_vs_combined_controls_cleaned_permisive.csv"
clean_df.to_csv(output_file, index=False)
print(f"\nCleaned results saved to: {output_file}")

# Print summary statistics of cleaned data
print(clean_df['log2fold_treated_control'].describe())

extreme_stat_records['groupID'].value_counts()

clean_df = df[
    # Remove extreme fold changes (keeping values between -5 and 5)
    (df['log2fold_treated_control'].abs() <= 5) &
    # Remove high dispersion
    (df['dispersion'] <= 10) &
    # Remove missing values
    (df['log2fold_treated_control'].notna()) &
    (df['pvalue'].notna()) &
    (df['padj'].notna()) &
    # Remove extreme test statistics
    (df['stat'].abs() <= 10000)
]

# Print summary of filtering
print(f"Original number of records: {len(df)}")
print(f"Number of records after cleaning: {len(clean_df)}")
print(f"Number of records removed: {len(df) - len(clean_df)}")

# Save cleaned results to new CSV file
output_file = "output/dexseq_results_PW1_vs_combined_controls_cleaned.csv"
clean_df.to_csv(output_file, index=False)
print(f"\nCleaned results saved to: {output_file}")

# Print summary statistics of cleaned data
print(clean_df['log2fold_treated_control'].describe())



