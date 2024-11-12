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

# # Local functions

def check_dexseq_results(results_file: str, output_dir: str = "qc_plots"):
    """
    Perform quality checks on DEXSeq results file.
    
    Args:
        results_file: Path to DEXSeq results CSV file
        output_dir: Directory to save QC plots
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Load results
    print("Loading results file...")
    df = pd.read_csv(results_file)
    
    # Basic statistics
    print("\n=== Basic Statistics ===")
    print(f"Total number of tests: {len(df)}")
    print(f"Number of NA p-values: {df['pvalue'].isna().sum()}")
    print(f"Number of NA adjusted p-values: {df['padj'].isna().sum()}")
    print(f"Number of significant results (padj < 0.05): {(df['padj'] < 0.05).sum()}")
    print(f"Number of significant results (padj < 0.1): {(df['padj'] < 0.1).sum()}")
    
    # Check for extreme fold changes
    fc_stats = df['log2fold_treated_control'].describe()
    print("\n=== Fold Change Statistics ===")
    print(fc_stats)
    
    # Identify potential problematic results
    print("\n=== Potential Issues ===")
    problematic = df[
        (df['log2fold_treated_control'].abs() > 5) |  # Extreme fold changes
        (df['dispersion'] > 10) |                     # High dispersion
        (df['pvalue'].isna()) |                      # Missing p-values
        (df['stat'].abs() > 10000)                   # Extreme test statistics
    ]
    print(f"Number of potentially problematic results: {len(problematic)}")
    if len(problematic) > 0:
        print("\nSample of problematic results:")
        print(problematic[['groupID', 'featureID', 'log2fold_treated_control', 'dispersion', 'pvalue', 'stat']].head())
    
    # Create plots
    plt.style.use('default')
    
    # 1. P-value distribution
    plt.figure(figsize=(10, 6))
    plt.hist(df['pvalue'].dropna(), bins=50, edgecolor='black')
    plt.title('P-value Distribution')
    plt.xlabel('P-value')
    plt.ylabel('Frequency')
    plt.savefig(f"{output_dir}/pvalue_distribution.png")
    plt.close()
    
    # 2. Volcano plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df['log2fold_treated_control'], 
               -np.log10(df['pvalue']),
               alpha=0.5)
    plt.title('Volcano Plot')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-log10(p-value)')
    plt.savefig(f"{output_dir}/volcano_plot.png")
    plt.close()
    
    # 3. MA plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df['exonBaseMean'],
               df['log2fold_treated_control'],
               alpha=0.5)
    plt.xscale('log')
    plt.title('MA Plot')
    plt.xlabel('Mean Expression')
    plt.ylabel('Log2 Fold Change')
    plt.savefig(f"{output_dir}/ma_plot.png")
    plt.close()
    
    # 4. Dispersion plot
    plt.figure(figsize=(10, 6))
    plt.scatter(df['exonBaseMean'],
               df['dispersion'],
               alpha=0.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.title('Dispersion Plot')
    plt.xlabel('Mean Expression')
    plt.ylabel('Dispersion')
    plt.savefig(f"{output_dir}/dispersion_plot.png")
    plt.close()
    
    # Check for the problematic gene/exon mentioned in the error
    problem_genes = ['ENSG00000285404.1', 'ENSG00000100150.19', 
                    'ENSG00000128245.15', 'ENSG00000252909.1']
    
    print("\n=== Checking Problematic Genes ===")
    for gene in problem_genes:
        gene_results = df[df['groupID'].str.contains(gene, na=False)]
        if len(gene_results) > 0:
            print(f"\nResults for {gene}:")
            print(gene_results[['featureID', 'log2fold_treated_control', 
                              'pvalue', 'padj', 'dispersion']].head())
    
    # Save problematic results to file
    if len(problematic) > 0:
        problematic.to_csv(f"{output_dir}/problematic_results.csv")
        print(f"\nProblematic results saved to {output_dir}/problematic_results.csv")
    
    # Return summary statistics
    return {
        'total_tests': len(df),
        'significant_005': (df['padj'] < 0.05).sum(),
        'significant_01': (df['padj'] < 0.1).sum(),
        'na_pvalues': df['pvalue'].isna().sum(),
        'problematic_count': len(problematic),
        'median_dispersion': df['dispersion'].median(),
        'median_fold_change': df['log2fold_treated_control'].median()
    }

# # Examind the data

# Example usage
results_file = "output/dexseq_results_PW1_vs_combined_controls.csv"
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

# # Select problematic genes

# List of problematic genes from the error message
problematic_genes = ['ENSG00000285404.1', 'ENSG00000100150.19', 
                    'ENSG00000128245.15', 'ENSG00000252909.1']

# Create masks for each condition
extreme_fc_mask = (df['log2fold_treated_control'].abs() > 5)
high_disp_mask = (df['dispersion'] > 10)
missing_vals_mask = (
    df['log2fold_treated_control'].isna() |
    df['pvalue'].isna() |
    df['padj'].isna()
)
extreme_stat_mask = (df['stat'].abs() > 10000)
problematic_genes_mask = df['groupID'].str.contains('|'.join(problematic_genes), regex=True)

# Print counts for each condition
print(f"Records with extreme fold changes (>5): {extreme_fc_mask.sum()}")
print(f"Records with high dispersion (>10): {high_disp_mask.sum()}")
print(f"Records with missing values: {missing_vals_mask.sum()}")
print(f"Records with extreme test statistics (>10000): {extreme_stat_mask.sum()}")
print(f"Records from problematic genes: {problematic_genes_mask.sum()}")

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
results = mg.querymany(stat_sel_no_version, scopes='ensembl.gene', fields='symbol', species='human')

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


# - log2fold change is tiny: 0.042537 (about 3% increase)
# - Very low dispersion: 0.001905
# - High stat value: 14819.648036
# - This means: The change is small, but it's happening very consistently across all replicates

# ## Examine extreme fold changes

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
results = mg.querymany(fc_sel_no_version, scopes='ensembl.gene', fields='symbol', species='human')

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



