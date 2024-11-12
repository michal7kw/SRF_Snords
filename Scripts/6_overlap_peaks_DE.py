# Converted from 6_overlap_peaks_DE.ipynb

import pandas as pd
import numpy as np
from intervaltree import IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Tuple
import os

# Set the working directory
working_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords"
os.chdir(working_dir)
print(f"Current working directory: {os.getcwd()}")

def load_and_process_data(peaks_file: str, dexseq_file: str, padj_threshold: float = 0.05) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load and process the peaks and DEXSeq data files.
    
    Args:
        peaks_file: Path to peaks CSV file
        dexseq_file: Path to DEXSeq results CSV file
        padj_threshold: Adjusted p-value threshold for significance
    
    Returns:
        Tuple of processed peaks and DEXSeq DataFrames
    """
    # Load data
    peaks_df = pd.read_csv(peaks_file)
    dexseq_df = pd.read_csv(dexseq_file)
    
    # Filter significant DEXSeq results
    dexseq_df = dexseq_df[dexseq_df['padj'] < padj_threshold].copy()
    
    return peaks_df, dexseq_df


def create_interval_tree(peaks_df: pd.DataFrame) -> IntervalTree:
    """
    Create an interval tree from peaks data for efficient overlap detection.
    """
    tree = IntervalTree()
    for _, row in peaks_df.iterrows():
        tree.addi(row['start'], row['end'], {
            'seqnames': row['seqnames'],
            'peak_id': row.name,
            'peak_name': row['peak_name']
        })
    return tree

def find_overlaps(dexseq_df: pd.DataFrame, peaks_tree: IntervalTree) -> List[dict]:
    """
    Find overlaps between DEXSeq exons and peaks.
    """
    overlaps = []
    
    for _, exon in dexseq_df.iterrows():
        chr_name = exon['genomicData.seqnames']
        start = exon['genomicData.start']
        end = exon['genomicData.end']
        
        # Find overlapping peaks
        overlapping = peaks_tree.overlap(start, end)
        
        for overlap in overlapping:
            if overlap.data['seqnames'] == chr_name:
                overlaps.append({
                    'exon_id': exon['featureID'],
                    'peak_id': overlap.data['peak_id'],
                    'peak_name': overlap.data['peak_name'],
                    'dexseq_name': exon['dexseq_name'],
                    'chromosome': chr_name,
                    'overlap_start': max(start, overlap.begin),
                    'overlap_end': min(end, overlap.end),
                    'overlap_length': min(end, overlap.end) - max(start, overlap.begin),
                    'exon_log2fc': exon['log2fold_treated_control'],
                    'exon_padj': exon['padj']
                })
    
    return overlaps

def plot_overlap_lengths(overlaps_df: pd.DataFrame, ax):
    """Plot distribution of overlap lengths."""
    sns.histplot(data=overlaps_df, x='overlap_length', bins=30, ax=ax)
    ax.set_title('Distribution of Overlap Lengths')
    ax.set_xlabel('Overlap Length (bp)')

def plot_fc_correlation(overlaps_df: pd.DataFrame, peaks_df: pd.DataFrame, ax):
    """Plot correlation between exon and peak fold changes."""
    peak_fc = peaks_df.loc[overlaps_df['peak_id']]['L2FC'].values
    sns.scatterplot(data=overlaps_df, x='exon_log2fc', y=peak_fc, ax=ax)
    ax.set_title('Exon Log2FC vs Peak Log2FC')
    ax.set_xlabel('Exon Log2FC')
    ax.set_ylabel('Peak Log2FC')
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))

def plot_chromosome_distribution(overlaps_df: pd.DataFrame, ax):
    """Plot distribution of overlaps across chromosomes."""
    overlaps_df['chromosome'].value_counts().plot(kind='bar', ax=ax)
    ax.set_title('Overlaps by Chromosome')
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Number of Overlaps')

def visualize_overlaps(overlaps_df: pd.DataFrame, peaks_df: pd.DataFrame, dexseq_df: pd.DataFrame):
    """
    Create visualizations for the overlap analysis.
    """
    # Create figure with multiple subplots
    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(2, 2)
    
    # 1. Distribution of overlap lengths
    ax1 = fig.add_subplot(gs[0, 0])
    plot_overlap_lengths(overlaps_df, ax1)
    
    # 2. Scatter plot of exon log2FC vs peak fold change
    ax2 = fig.add_subplot(gs[0, 1])
    plot_fc_correlation(overlaps_df, peaks_df, ax2)
    
    # 3. Chromosome distribution of overlaps
    ax3 = fig.add_subplot(gs[1, 0])
    plot_chromosome_distribution(overlaps_df, ax3)
    
    plt.tight_layout()
    return fig

# Example usage
peaks_file = "./DATA/Peak.csv"
dexseq_file = "./output/dexseq_results_PW1_vs_combined_controls_cleaned.csv"
output_prefix = "./output/overlap_analysis"

# overlaps_df = main(peaks_file, dexseq_file, output_prefix)

# Load and process data
peaks_df, dexseq_df = load_and_process_data(peaks_file, dexseq_file)

pd.set_option('display.max_columns', None)
peaks_df.head()

peaks_df = peaks_df.rename(columns={"Row.names": "peak_name"})

print(list(peaks_df['seqnames'][:5]))
print(list(peaks_df.iloc[:5].index))

dexseq_df.head()

dexseq_df = dexseq_df.rename(columns={"Unnamed: 0": "dexseq_name"})

dexseq_df.columns

print(list(dexseq_df['genomicData.seqnames'][:5]))
print(list(dexseq_df['featureID'][:5]))


# Create interval tree for efficient overlap detection
peaks_tree = create_interval_tree(peaks_df)

list(peaks_tree.items())[:10]



# Find overlaps
overlaps = find_overlaps(dexseq_df, peaks_tree)
overlaps_df = pd.DataFrame(overlaps)

overlaps_df.head()

overlaps_df.shape

# Save overlaps to file
overlaps_df.to_csv(f"{output_prefix}_overlaps.csv", index=False)

fig, ax = plt.subplots(figsize=(10, 6))
plot_overlap_lengths(overlaps_df, ax)


fig, ax = plt.subplots(figsize=(10, 6))
plot_fc_correlation(overlaps_df, peaks_df, ax)

fig, ax = plt.subplots(figsize=(10, 6))
plot_chromosome_distribution(overlaps_df, ax)

total_peaks = len(peaks_df)
total_exons = len(dexseq_df)
peaks_with_overlaps = len(overlaps_df['peak_id'].unique())
exons_with_overlaps = len(overlaps_df['exon_id'].unique())

print(f"Total Peaks: {total_peaks}")
print(f"Total Diff. Expressed Exons: {total_exons}")
print(f"Peaks with Overlaps: {peaks_with_overlaps}")
print(f"Exons with Overlaps: {exons_with_overlaps}")
print(f"Total Overlap Events: {len(overlaps_df)}")

def find_overlaps2(dexseq_df: pd.DataFrame, peaks_tree: IntervalTree) -> List[dict]:
    """
    Find overlaps and nearby peaks for DEXSeq exons.
    Also calculates distance to nearest peak for non-overlapping cases.
    """
    results = []
    
    for index, exon in dexseq_df.iterrows():
        chr_name = exon['genomicData.seqnames']
        exon_start = exon['genomicData.start']
        exon_end = exon['genomicData.end']
        exon_center = (exon_start + exon_end) / 2
        
        # Extend search window by 10kb in each direction
        search_start = exon_start - 10000
        search_end = exon_end + 10000
        
        # Find peaks in extended window
        nearby = peaks_tree.overlap(search_start, search_end)
        if index == 0 and len(nearby) > 0:
            print(list(nearby))
        
        for peak in nearby:
            if peak.data['seqnames'] == chr_name:
                peak_center = (peak.begin + peak.end) / 2
                
                # Calculate distance (negative if peak is upstream, positive if downstream)
                distance = peak_center - exon_center
                
                # Calculate overlap (if any)
                overlap_start = max(exon_start, peak.begin)
                overlap_end = min(exon_end, peak.end)
                overlap_length = max(0, overlap_end - overlap_start)
                
                results.append({
                    'exon_id': exon['featureID'],
                    'peak_id': peak.data['peak_id'],
                    'peak_name': peak.data['peak_name'],
                    'dexseq_name': exon['dexseq_name'],
                    'chromosome': chr_name,
                    'distance_to_peak': distance,
                    'overlap_length': overlap_length,
                    'exon_log2fc': exon['log2fold_treated_control'],
                    'exon_padj': exon['padj']
                })
    
    return results

def plot_distance_vs_fc(overlaps_df: pd.DataFrame, ax):
    """Plot relationship between peak distance and exon fold change."""
    # Convert distance to kb for better visualization
    overlaps_df['distance_kb'] = overlaps_df['distance_to_peak'] / 1000
    
    # Create scatter plot
    sns.scatterplot(
        data=overlaps_df,
        x='distance_kb',
        y='exon_log2fc',
        alpha=0.5,
        ax=ax
    )
    
    # Add trend line
    sns.regplot(
        data=overlaps_df,
        x='distance_kb',
        y='exon_log2fc',
        scatter=False,
        color='red',
        ax=ax
    )
    
    ax.set_title('Exon Log2FC vs Distance to Peak')
    ax.set_xlabel('Distance to Peak (kb)')
    ax.set_ylabel('Exon Log2FC')
    
    # Add vertical line at x=0 to mark the exon position
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    # Add horizontal line at y=0 to mark no change in expression
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

def visualize_distance_analysis(overlaps_df: pd.DataFrame):
    """
    Create visualizations for the distance analysis.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 1. Distance vs Fold Change scatter plot
    plot_distance_vs_fc(overlaps_df, ax1)
    
    # 2. Distribution of distances
    sns.histplot(
        data=overlaps_df,
        x='distance_to_peak',
        bins=50,
        ax=ax2
    )
    ax2.set_title('Distribution of Peak Distances')
    ax2.set_xlabel('Distance to Peak (bp)')
    
    plt.tight_layout()
    return fig

# # Find overlaps and distances
overlaps = find_overlaps2(dexseq_df, peaks_tree)
overlaps_df = pd.DataFrame(overlaps)

overlaps_df.head()

# Create visualization
fig = visualize_distance_analysis(overlaps_df)
plt.show()



