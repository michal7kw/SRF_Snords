# Standard library imports
import os
import signal
import time
from typing import Any, Dict, List, Optional, Tuple
import matplotlib.pyplot as plt
import numpy as np
import requests
import json
import requests
import json
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import numpy as np
from numba import jit, prange
import pandas as pd
import multiprocessing
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
import asyncio
import aiohttp
from typing import Dict, List, Tuple, Optional, Set
import time
from collections import defaultdict
import sqlite3
from functools import partial
import hashlib
import os
import requests
import json
import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import requests
import json
from multiprocessing import Pool
import pandas as pd
from typing import List, Dict
from intervaltree import IntervalTree
import numpy as np
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import requests
import json

# Third party imports
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
import matplotlib.pyplot as plt
import multiprocessing
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq
from intervaltree import IntervalTree

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

def visualize_top_alignments(enriched_overlaps_df: pd.DataFrame, n: int = 5):
    """
    Visualize the top n alignments from the enriched overlaps dataframe.
    
    Args:
        enriched_overlaps_df (pd.DataFrame): DataFrame containing alignment scores and gene information
        n (int): Number of top alignments to display (default=5)
    """
    # Get top n rows by sum_als1_score
    top_alignments = enriched_overlaps_df.nlargest(n, 'sum_als1_score')
    
    print(f"Top {n} alignments by ALS1 score:\n")
    
    for _, row in top_alignments.iterrows():
        # Print gene name and basic information
        print(f"Gene: {row['gene_name']} ({row['geneID']})")
        print(f"Exon: {row['exon_id']}")
        print(f"Peak: {row['peak_name']}")
        
        # Print scores
        print(f"ALS1 Scores:")
        if row['previous_exon_als1_score'] > 0:
            print(f"  Previous exon: {row['previous_exon_als1_score']:.1f}")
        print(f"  Current exon:  {row['current_exon_als1_score']:.1f}")
        if row['next_exon_als1_score'] > 0:
            print(f"  Next exon:     {row['next_exon_als1_score']:.1f}")
        print(f"  Total:         {row['sum_als1_score']:.1f}")
        
        print(f"ALS2 Scores:")
        if row['previous_exon_als2_score'] > 0:
            print(f"  Previous exon: {row['previous_exon_als2_score']:.1f}")
        print(f"  Current exon:  {row['current_exon_als2_score']:.1f}")
        if row['next_exon_als2_score'] > 0:
            print(f"  Next exon:     {row['next_exon_als2_score']:.1f}")
        print(f"  Total:         {row['sum_als2_score']:.1f}")
        
        # Print additional information
        print(f"Distance to peak: {row['distance_kb']:.2f} kb")
        if row['overlap_length'] > 0:
            print(f"Overlap length: {row['overlap_length']} bp")
        print(f"Log2 fold change: {row['exon_log2fc']:.3f}")
        print("-" * 50 + "\n")

############################### Find Overlaps #########################################################################################

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
            # IMPORTANT: Only include if chromosomes match AND there is actual overlap
            if (overlap.data['seqnames'] == chr_name and 
                max(start, overlap.begin) < min(end, overlap.end)):  # This ensures real overlap
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

def find_overlaps3(dexseq_df: pd.DataFrame, peaks_tree: IntervalTree, max_distance: int = 2000) -> List[dict]:
    results = []
    
    # Group exons by gene
    for gene_name, gene_exons in dexseq_df.groupby('dexseq_name'):
        # Get gene boundaries
        gene_start = gene_exons['genomicData.start'].min()
        gene_end = gene_exons['genomicData.end'].max()
        
        # Sort exons by position
        gene_exons = gene_exons.sort_values('genomicData.start')
        
        # Process each exon
        for idx, exon in gene_exons.iterrows():
            chr_name = exon['genomicData.seqnames']
            exon_start = exon['genomicData.start']
            exon_end = exon['genomicData.end']
            exon_center = (exon_start + exon_end) / 2
            
            # Find all exons within max_distance
            nearby_exons = gene_exons[
                (gene_exons.index != idx) &  # Exclude current exon
                (
                    # Exon start or end is within max_distance
                    ((gene_exons['genomicData.start'] >= exon_start - max_distance) &
                     (gene_exons['genomicData.start'] <= exon_end + max_distance)) |
                    ((gene_exons['genomicData.end'] >= exon_start - max_distance) &
                     (gene_exons['genomicData.end'] <= exon_end + max_distance))
                )
            ]
            
            # Define search window to include all nearby exons, but constrained by gene boundaries
            if len(nearby_exons) > 0:
                search_start = max(gene_start, min(nearby_exons['genomicData.start'].min(), exon_start))
                search_end = min(gene_end, max(nearby_exons['genomicData.end'].max(), exon_end))
            else:
                search_start = max(gene_start, exon_start)
                search_end = min(gene_end, exon_end)
            
            # Find peaks in the search window
            nearby_peaks = peaks_tree.overlap(search_start, search_end)
            
            for peak in nearby_peaks:
                if peak.data['seqnames'] == chr_name:
                    # Check if peak is within max_distance from exon center
                    peak_center = (peak.begin + peak.end) / 2
                    distance = peak_center - exon_center
                    
                    if abs(distance) <= max_distance:  # Only include peaks within max_distance
                        # Calculate overlap with current exon
                        overlap_start = max(exon_start, peak.begin)
                        overlap_end = min(exon_end, peak.end)
                        overlap_length = max(0, overlap_end - overlap_start)
                        
                        # Get information about nearby exons
                        nearby_exons_info = [{
                            'exon_id': ne['featureID'],
                            'distance': min(
                                abs(ne['genomicData.start'] - exon_end),
                                abs(ne['genomicData.end'] - exon_start)
                            )
                        } for _, ne in nearby_exons.iterrows()]
                        
                        results.append({
                            'exon_id': exon['featureID'],
                            'peak_id': peak.data['peak_id'],
                            'peak_name': peak.data['peak_name'],
                            'dexseq_name': exon['dexseq_name'],
                            'chromosome': chr_name,
                            'distance_to_peak': distance,
                            'overlap_length': overlap_length,
                            'exon_log2fc': exon['log2fold_treated_control'],
                            'exon_padj': exon['padj'],
                            'nearby_exons': nearby_exons_info,
                            'n_nearby_exons': len(nearby_exons_info)
                        })
    
    return results

#########################################################################################################################################

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

def get_exon_sequence(chromosome: str, start: int, end: int) -> Optional[str]:
    """
    Fetch genomic sequence for given coordinates using Ensembl REST API.
    """
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?version=38"  # Added version parameter
    max_retries = 3

    for attempt in range(max_retries):
        try:
            r = requests.get(server + ext, 
                           headers={"Content-Type": "text/plain"},
                           timeout=30)
            if r.ok:
                sequence = r.text.strip()
                if sequence:
                    return sequence
            time.sleep(1)
        except requests.exceptions.RequestException as e:
            print(f"Request failed for {chromosome}:{start}-{end}, attempt {attempt + 1}: {str(e)}")
    return None

def align_sequence(query: Seq, target: Seq) -> Dict[str, Any]:
    """
    Perform local alignment and return the best score and alignment.
    Added input validation and better error handling.
    """
    if not query or not target:
        return {"score": 0, "alignment": None}
        
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -100
    aligner.extend_gap_score = -100
    
    try:
        alignments = list(aligner.align(query, target))  # Convert iterator to list
        if alignments:
            best_alignment = max(alignments, key=lambda a: a.score)
            return {"score": best_alignment.score, "alignment": best_alignment}
        return {"score": 0, "alignment": None}
    except Exception as e:
        print(f"Error during alignment: {str(e)}")
        return {"score": 0, "alignment": None}

def analyze_als_complementarity(row: pd.Series, dexseq_df: pd.DataFrame, 
                              als1: str = 'CTTTTCCAAGGAATGTT', 
                              als2: str = 'TTCCGATGAGAATGACG') -> dict:
    """
    Analyze ALS complementarity for an exon and its neighbors using local alignment.
    Now tracks original exon and neighbors separately.
    """
    if not isinstance(row, pd.Series) or not isinstance(dexseq_df, pd.DataFrame):
        print("Invalid input types")
        return {}

    try:
        current_gene = row['dexseq_name'].split(':')[0]
        current_exon_num = int(row['exon_id'].replace('E', ''))
        
        # Find neighboring exons from the same gene
        gene_exons = dexseq_df[dexseq_df['dexseq_name'].str.startswith(current_gene)].copy()
        if gene_exons.empty:
            return {}
            
        gene_exons['exon_num'] = gene_exons['featureID'].str.replace('E', '').astype(int)
        gene_exons = gene_exons.sort_values('exon_num')
        
        # Categorize exons as current, previous, or next
        sequences = {
            'current': None,
            'previous': None,
            'next': None
        }
        
        for idx, exon in gene_exons.iterrows():
            exon_num_diff = exon['exon_num'] - current_exon_num
            if exon_num_diff == 0:
                category = 'current'
            elif exon_num_diff == -1:
                category = 'previous'
            elif exon_num_diff == 1:
                category = 'next'
            else:
                continue
                
            seq = get_exon_sequence(
                str(exon['genomicData.seqnames']),
                int(exon['genomicData.start']),
                int(exon['genomicData.end'])
            )
            if seq:
                sequences[category] = {
                    'sequence': seq.strip().upper(),
                    'exon_id': exon['featureID']
                }
        
        results = {
            'current_exon': None,
            'previous_exon': None,
            'next_exon': None
        }
        
        # Process each exon category
        for category, seq_data in sequences.items():
            if seq_data:
                try:
                    forward_seq = Seq(seq_data['sequence'])
                    als1_seq = Seq(als1.upper())
                    als2_seq = Seq(als2.upper())
                    
                    forward_als1 = align_sequence(forward_seq, als1_seq)
                    forward_als2 = align_sequence(forward_seq, als2_seq)
                    
                    results[f'{category}_exon'] = {
                        'exon_id': seq_data['exon_id'],
                        'als1_score': float(forward_als1["score"]),
                        'als2_score': float(forward_als2["score"]),
                        'als1_alignment': str(forward_als1["alignment"]) if forward_als1["alignment"] else "",
                        'als2_alignment': str(forward_als2["alignment"]) if forward_als2["alignment"] else "",
                        'sequence': seq_data['sequence'],
                        'sequence_length': len(seq_data['sequence'])
                    }
                except Exception as e:
                    print(f"Error processing {category} exon: {str(e)}")
                    
        return results
        
    except Exception as e:
        print(f"Error in analyze_als_complementarity: {str(e)}")
        return {}

############################### analyze overlaps with als ##################################################################################################

def process_chunk(chunk_data: tuple) -> list:
    """
    Process a chunk of overlaps data with separate tracking of original and neighboring exons.
    """
    try:
        chunk, dexseq_df = chunk_data
        results = []
        
        print(f"Processing chunk of size {len(chunk)}")
        for idx, row in chunk.iterrows():
            try:
                als_analysis = analyze_als_complementarity(row, dexseq_df)
                
                row_dict = row.to_dict()
                
                # Add scores for current exon and neighbors separately
                if als_analysis.get('current_exon'):
                    row_dict.update({
                        'current_exon_id': als_analysis['current_exon']['exon_id'],
                        'current_exon_als1_score': als_analysis['current_exon']['als1_score'],
                        'current_exon_als2_score': als_analysis['current_exon']['als2_score'],
                        'current_exon_sequence_length': als_analysis['current_exon']['sequence_length']
                    })
                
                if als_analysis.get('previous_exon'):
                    row_dict.update({
                        'previous_exon_id': als_analysis['previous_exon']['exon_id'],
                        'previous_exon_als1_score': als_analysis['previous_exon']['als1_score'],
                        'previous_exon_als2_score': als_analysis['previous_exon']['als2_score'],
                        'previous_exon_sequence_length': als_analysis['previous_exon']['sequence_length']
                    })
                
                if als_analysis.get('next_exon'):
                    row_dict.update({
                        'next_exon_id': als_analysis['next_exon']['exon_id'],
                        'next_exon_als1_score': als_analysis['next_exon']['als1_score'],
                        'next_exon_als2_score': als_analysis['next_exon']['als2_score'],
                        'next_exon_sequence_length': als_analysis['next_exon']['sequence_length']
                    })
                
                results.append(row_dict)
                
            except Exception as e:
                print(f"Error processing row with exon {row['exon_id']}: {str(e)}")
                results.append(row.to_dict())
        
        return results
        
    except Exception as e:
        print(f"Error in process_chunk: {str(e)}")
        return []
    
def analyze_overlaps_with_als(overlaps_df: pd.DataFrame, dexseq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze ALS complementarity for all overlapping exons and their neighbors using parallel processing.
    Added input validation and better error handling.
    """
    if overlaps_df.empty or dexseq_df.empty:
        print("Empty input DataFrame(s)")
        return pd.DataFrame()
        
    try:
        print("Starting ALS complementarity analysis...")
        
        # Determine number of cores to use (leave one free for system)
        n_cores = max(1, multiprocessing.cpu_count() - 1)
        print(f"Using {n_cores} CPU cores for parallel processing")
        
        # Split the dataframe into chunks
        chunk_size = max(1, len(overlaps_df) // n_cores)
        chunks = [overlaps_df.iloc[i:i + chunk_size] for i in range(0, len(overlaps_df), chunk_size)]
        print(f"Split data into {len(chunks)} chunks of approximately {chunk_size} rows each")
        
        # Prepare data for parallel processing
        chunk_data = [(chunk, dexseq_df) for chunk in chunks]
        
        # Process chunks in parallel with timeout
        all_results = []
        print("Starting parallel processing...")
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            chunk_results = list(executor.map(process_chunk, chunk_data, timeout=3600))  # 1 hour timeout
            for results in chunk_results:
                if results:  # Check if results is not empty
                    all_results.extend(results)
        
        if not all_results:
            print("No results generated from analysis")
            return pd.DataFrame()
            
        print(f"Analysis complete. Processed {len(all_results)} total entries")
        return pd.DataFrame(all_results)
        
    except Exception as e:
        print(f"Error in analyze_overlaps_with_als: {str(e)}")
        return pd.DataFrame()

############################### analyze overlaps with als #########################################################################################
def analyze_overlaps_with_als_optimized(overlaps_df: pd.DataFrame, dexseq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze ALS complementarity for all overlapping exons and their neighbors using parallel processing.
    Optimized for high-core-count systems.
    """
    if overlaps_df.empty or dexseq_df.empty:
        print("Empty input DataFrame(s)")
        return pd.DataFrame()
        
    try:
        print("Starting ALS complementarity analysis...")
        
        # Use more cores but still leave some headroom for system processes
        n_cores = max(1, multiprocessing.cpu_count() - 2)  # Changed from -1 to -2
        print(f"Using {n_cores} CPU cores for parallel processing")
        
        # Optimize chunk size for better load balancing
        # Smaller chunks allow better distribution across cores
        optimal_chunk_size = max(1, min(
            len(overlaps_df) // (n_cores * 4),  # Increased parallelization
            100  # Cap maximum chunk size for better load balancing
        ))
        chunks = [overlaps_df.iloc[i:i + optimal_chunk_size] 
                 for i in range(0, len(overlaps_df), optimal_chunk_size)]
        
        print(f"Split data into {len(chunks)} chunks of approximately {optimal_chunk_size} rows each")
        
        # Prepare data for parallel processing
        chunk_data = [(chunk, dexseq_df) for chunk in chunks]
        
        # Process chunks in parallel with optimized settings
        all_results = []
        print("Starting parallel processing...")
        
        # Use a context manager for better resource management
        ctx = multiprocessing.get_context('spawn')  # More stable than fork
        with ProcessPoolExecutor(
            max_workers=n_cores,
            mp_context=ctx,
            initializer=_init_worker  # Add worker initialization
        ) as executor:
            # Use as_completed for better progress tracking and fault tolerance
            futures = {executor.submit(process_chunk_optimized, data): i 
                      for i, data in enumerate(chunk_data)}
            
            # Process results as they complete
            for future in concurrent.futures.as_completed(futures):
                chunk_idx = futures[future]
                try:
                    results = future.result(timeout=3600)  # 1 hour timeout per chunk
                    if results:
                        all_results.extend(results)
                    print(f"Completed chunk {chunk_idx + 1}/{len(chunks)}")
                except Exception as e:
                    print(f"Error processing chunk {chunk_idx}: {str(e)}")
        
        if not all_results:
            print("No results generated from analysis")
            return pd.DataFrame()
            
        print(f"Analysis complete. Processed {len(all_results)} total entries")
        return pd.DataFrame(all_results)
        
    except Exception as e:
        print(f"Error in analyze_overlaps_with_als: {str(e)}")
        return pd.DataFrame()

def _init_worker():
    """
    Initialize worker process to ignore SIGINT.
    This prevents keyboard interrupts from affecting worker processes.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def process_chunk_optimized(chunk_data: tuple) -> list:
    """
    Process a chunk of overlaps data with optimized memory handling.
    """
    try:
        chunk, dexseq_df = chunk_data
        results = []
        
        # Pre-filter dexseq_df to only include relevant genes for this chunk
        genes_in_chunk = {row['dexseq_name'].split(':')[0] 
                         for _, row in chunk.iterrows()}
        filtered_dexseq = dexseq_df[
            dexseq_df['dexseq_name'].str.split(':').str[0].isin(genes_in_chunk)
        ].copy()
        
        for idx, row in chunk.iterrows():
            try:
                als_analysis = analyze_als_complementarity(row, filtered_dexseq)
                
                # Convert row to dict only once
                row_dict = row.to_dict()
                
                # Update dictionary with a single operation
                updates = {}
                
                for exon_type in ['current', 'previous', 'next']:
                    if als_analysis.get(f'{exon_type}_exon'):
                        exon_data = als_analysis[f'{exon_type}_exon']
                        updates.update({
                            f'{exon_type}_exon_id': exon_data['exon_id'],
                            f'{exon_type}_exon_als1_score': exon_data['als1_score'],
                            f'{exon_type}_exon_als2_score': exon_data['als2_score'],
                            f'{exon_type}_exon_sequence_length': exon_data['sequence_length']
                        })
                
                row_dict.update(updates)
                results.append(row_dict)
                
            except Exception as e:
                print(f"Error processing row with exon {row['exon_id']}: {str(e)}")
                results.append(row.to_dict())
        
        return results
        
    except Exception as e:
        print(f"Error in process_chunk: {str(e)}")
        return []

############################### Find Overlaps Alternative ##########################################################################################
def process_gene_group(args) -> List[dict]:
    """
    Helper function to process a single gene group.
    Takes a tuple of (gene_name, gene_exons_df, peaks_tree, max_distance)
    """
    gene_name, gene_exons, peaks_tree, max_distance = args
    results = []
    
    # Get gene boundaries
    gene_start = gene_exons['genomicData.start'].min()
    gene_end = gene_exons['genomicData.end'].max()
    
    # Sort exons by position
    gene_exons = gene_exons.sort_values('genomicData.start')
    
    # Process each exon
    for idx, exon in gene_exons.iterrows():
        chr_name = exon['genomicData.seqnames']
        exon_start = exon['genomicData.start']
        exon_end = exon['genomicData.end']
        exon_center = (exon_start + exon_end) / 2
        
        # Find all exons within max_distance
        nearby_exons = gene_exons[
            (gene_exons.index != idx) &
            (
                ((gene_exons['genomicData.start'] >= exon_start - max_distance) &
                 (gene_exons['genomicData.start'] <= exon_end + max_distance)) |
                ((gene_exons['genomicData.end'] >= exon_start - max_distance) &
                 (gene_exons['genomicData.end'] <= exon_end + max_distance))
            )
        ]
        
        # Define search window
        if len(nearby_exons) > 0:
            search_start = max(gene_start, min(nearby_exons['genomicData.start'].min(), exon_start))
            search_end = min(gene_end, max(nearby_exons['genomicData.end'].max(), exon_end))
        else:
            search_start = max(gene_start, exon_start)
            search_end = min(gene_end, exon_end)
        
        # Find peaks in the search window
        nearby_peaks = peaks_tree.overlap(search_start, search_end)
        
        for peak in nearby_peaks:
            if peak.data['seqnames'] == chr_name:
                peak_center = (peak.begin + peak.end) / 2
                distance = peak_center - exon_center
                
                if abs(distance) <= max_distance:
                    overlap_start = max(exon_start, peak.begin)
                    overlap_end = min(exon_end, peak.end)
                    overlap_length = max(0, overlap_end - overlap_start)
                    
                    nearby_exons_info = [{
                        'exon_id': ne['featureID'],
                        'distance': min(
                            abs(ne['genomicData.start'] - exon_end),
                            abs(ne['genomicData.end'] - exon_start)
                        )
                    } for _, ne in nearby_exons.iterrows()]
                    
                    results.append({
                        'exon_id': exon['featureID'],
                        'peak_id': peak.data['peak_id'],
                        'peak_name': peak.data['peak_name'],
                        'dexseq_name': exon['dexseq_name'],
                        'chromosome': chr_name,
                        'distance_to_peak': distance,
                        'overlap_length': overlap_length,
                        'exon_log2fc': exon['log2fold_treated_control'],
                        'exon_padj': exon['padj'],
                        'nearby_exons': nearby_exons_info,
                        'n_nearby_exons': len(nearby_exons_info)
                    })
    
    return results

def find_overlaps3_parallel(dexseq_df: pd.DataFrame, peaks_tree: IntervalTree, 
                          max_distance: int = 2000, n_processes: int = None) -> List[dict]:
    """
    Parallelized version of find_overlaps3 function.
    
    Args:
        dexseq_df: DataFrame containing exon information
        peaks_tree: IntervalTree containing peak information
        max_distance: Maximum distance to consider for overlaps
        n_processes: Number of processes to use (defaults to CPU count)
    
    Returns:
        List of dictionaries containing overlap information
    """
    # Group exons by gene
    grouped = list(dexseq_df.groupby('dexseq_name'))
    
    # Create argument tuples for each gene group
    process_args = [(name, group, peaks_tree, max_distance) for name, group in grouped]
    
    # Create process pool and map the work
    with Pool(processes=n_processes) as pool:
        results_nested = pool.map(process_gene_group, process_args)
    
    # Flatten results
    results = [item for sublist in results_nested for item in sublist]
    
    return results

# Example usage:
# results = find_overlaps3_parallel(dexseq_df, peaks_tree, max_distance=2000, n_processes=4)

############################### analyze overlaps with als Alternative #########################################################################################
class SequenceCache:
    def __init__(self, cache_file="sequence_cache.db"):
        self.cache_file = cache_file
        self.setup_db()
        
    def setup_db(self):
        with sqlite3.connect(self.cache_file) as conn:
            conn.execute("""
                CREATE TABLE IF NOT EXISTS sequences (
                    key TEXT PRIMARY KEY,
                    sequence TEXT,
                    timestamp REAL
                )
            """)
            conn.execute("CREATE INDEX IF NOT EXISTS idx_key ON sequences(key)")
    
    def get_key(self, chrom: str, start: int, end: int) -> str:
        return f"{chrom}:{start}-{end}"
    
    def get(self, chrom: str, start: int, end: int) -> Optional[str]:
        key = self.get_key(chrom, start, end)
        try:
            with sqlite3.connect(self.cache_file) as conn:
                result = conn.execute(
                    "SELECT sequence FROM sequences WHERE key = ?", 
                    (key,)
                ).fetchone()
                return result[0] if result else None
        except:
            return None
    
    def set(self, chrom: str, start: int, end: int, sequence: str):
        key = self.get_key(chrom, start, end)
        try:
            with sqlite3.connect(self.cache_file) as conn:
                conn.execute(
                    "INSERT OR REPLACE INTO sequences (key, sequence, timestamp) VALUES (?, ?, ?)",
                    (key, sequence, time.time())
                )
        except:
            pass

class AsyncSequenceFetcher:
    def __init__(self, max_concurrent=50, cache_file="sequence_cache.db"):
        self.max_concurrent = max_concurrent
        self.cache = SequenceCache(cache_file)
        self.semaphore = asyncio.Semaphore(max_concurrent)
        
    async def fetch_sequence(self, session: aiohttp.ClientSession, chrom: str, start: int, end: int) -> Optional[str]:
        # Check cache first
        cached = self.cache.get(chrom, start, end)
        if cached:
            return cached
            
        url = f"https://rest.ensembl.org/sequence/region/human/{chrom}:{start}..{end}:1"
        headers = {"Content-Type": "text/plain"}
        
        async with self.semaphore:
            try:
                async with session.get(url, headers=headers) as response:
                    if response.status == 200:
                        sequence = await response.text()
                        sequence = sequence.strip()
                        # Cache the result
                        self.cache.set(chrom, start, end, sequence)
                        return sequence
            except:
                pass
        return None
        
    async def fetch_batch(self, coords: List[Tuple[str, int, int]]) -> Dict[str, str]:
        async with aiohttp.ClientSession() as session:
            tasks = [
                self.fetch_sequence(session, chrom, start, end)
                for chrom, start, end in coords
            ]
            sequences = await asyncio.gather(*tasks)
            return {
                f"{chrom}:{start}-{end}": seq
                for (chrom, start, end), seq in zip(coords, sequences)
                if seq
            }

@jit(nopython=True, parallel=True, fastmath=True)
def calculate_alignment_scores_batch(sequences: np.ndarray, als_seq: np.ndarray) -> np.ndarray:
    """
    Vectorized alignment score calculation for a batch of sequences
    """
    n_sequences = sequences.shape[0]
    scores = np.zeros(n_sequences, dtype=np.float32)
    
    for i in prange(n_sequences):
        seq = sequences[i]
        m, n = len(seq), len(als_seq)
        score_matrix = np.zeros((m + 1, n + 1), dtype=np.float32)
        
        for j in range(1, m + 1):
            for k in range(1, n + 1):
                match = score_matrix[j-1, k-1] + (1.0 if seq[j-1] == als_seq[k-1] else -1.0)
                delete = score_matrix[j-1, k] - 1.0
                insert = score_matrix[j, k-1] - 1.0
                score_matrix[j, k] = max(0.0, match, delete, insert)
                scores[i] = max(scores[i], score_matrix[j, k])
    
    return scores

def process_sequences(sequences: List[str], als1: str, als2: str) -> List[Dict[str, float]]:
    if not sequences:
        return []
    
    valid_sequences = [seq for seq in sequences if seq]
    if not valid_sequences:
        return []
        
    try:
        seq_arrays = np.array([list(seq.upper()) for seq in valid_sequences], dtype=np.uint8)
        als1_array = np.array(list(als1.upper()), dtype=np.uint8)
        als2_array = np.array(list(als2.upper()), dtype=np.uint8)
        
        als1_scores = calculate_alignment_scores_batch(seq_arrays, als1_array)
        als2_scores = calculate_alignment_scores_batch(seq_arrays, als2_array)
        
        return [
            {
                'als1_score': float(als1_scores[i]),
                'als2_score': float(als2_scores[i])
            }
            for i in range(len(valid_sequences))
        ]
    except Exception as e:
        print(f"Error processing sequences: {str(e)}")
        return []

def process_chunk_optimized_2(args: Tuple[pd.DataFrame, pd.DataFrame, str, str, str]) -> List[Dict]:
    chunk, dexseq_df, als1, als2, cache_file = args
    results = []
    
    try:
        # Join with dexseq_df to get coordinates
        merged_data = pd.merge(
            chunk,
            dexseq_df[['featureID', 'genomicData.seqnames', 'genomicData.start', 'genomicData.end']],
            left_on='exon_id',
            right_on='featureID',
            how='inner'
        )
        
        print(f"Processing chunk with {len(merged_data)} entries after merging")
        
        # Group by chromosome
        for chrom, group in merged_data.groupby('genomicData.seqnames'):
            coords = [
                (str(chrom), int(row['genomicData.start']), int(row['genomicData.end']))
                for _, row in group.iterrows()
            ]
            
            # Create async fetcher
            fetcher = AsyncSequenceFetcher(cache_file=cache_file)
            
            # Run async fetch in event loop
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
            sequences = loop.run_until_complete(fetcher.fetch_batch(coords))
            loop.close()
            
            # Process sequences while we have them
            for idx, row in group.iterrows():
                key = f"{chrom}:{int(row['genomicData.start'])}-{int(row['genomicData.end'])}"
                if key in sequences:
                    seq_scores = process_sequences([sequences[key]], als1, als2)
                    if seq_scores:
                        row_dict = row.to_dict()
                        row_dict.update(seq_scores[0])
                        results.append(row_dict)
    
    except Exception as e:
        print(f"Error in chunk processing: {str(e)}")
        import traceback
        traceback.print_exc()
    
    return results

def analyze_overlaps_with_als_optimized_2(overlaps_df: pd.DataFrame, 
                                      dexseq_df: pd.DataFrame,
                                      als1: str = 'CTTTTCCAAGGAATGTT',
                                      als2: str = 'TTCCGATGAGAATGACG',
                                      batch_size: int = 1000,
                                      cache_file: str = "sequence_cache.db") -> pd.DataFrame:
    if overlaps_df.empty or dexseq_df.empty:
        return pd.DataFrame()
    
    try:
        print(f"Input shapes - overlaps: {overlaps_df.shape}, dexseq: {dexseq_df.shape}")
        
        # Use more cores
        n_cores = max(1, multiprocessing.cpu_count())
        print(f"Using {n_cores} cores for processing")
        
        # Create larger chunks
        n_chunks = max(1, min(len(overlaps_df) // batch_size, n_cores * 4))
        chunks = np.array_split(overlaps_df, n_chunks)
        print(f"Split data into {len(chunks)} chunks")
        
        # Prepare args
        chunk_args = [(chunk, dexseq_df, als1, als2, cache_file) for chunk in chunks]
        
        # Process chunks
        all_results = []
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            for i, results in enumerate(executor.map(process_chunk_optimized_2, chunk_args)):
                all_results.extend(results)
                print(f"Processed chunk {i + 1}/{len(chunks)} with {len(results)} results")
        
        print(f"Analysis complete. Processed {len(all_results)} entries")
        return pd.DataFrame(all_results)
        
    except Exception as e:
        print(f"Error in analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        return pd.DataFrame()

# ##########################################################################################################################################################################################################

def examine_exon(db, gene_id="ENSG00000114744.9", exon_id="E001"):
    """
    Examine the gene's exon details
    """
    print(f"Searching for gene_id: {gene_id}, exon_id: {exon_id}")
    
    # Get all exonic parts
    exons = []
    for feature in db.features_of_type('exonic_part', order_by='start'):
        gene_ids = feature.attributes.get('gene_id', [])
        if gene_id in gene_ids:
            exons.append(feature)
    
    print(f"Found {len(exons)} exons for this gene")
    
    if not exons:
        print("No exons found for this gene ID")
        return
        
    # Sort exons by position
    exons.sort(key=lambda x: x.start)
    
    # Debug print first few exons
    # print("\nFirst few exons found:")
    # for i, exon in enumerate(exons[:3]):
    #     print(f"Exon {i}: ID={exon.id}, Number={exon.attributes.get('exonic_part_number', ['None'])}")
    
    found = False
    for i, exon in enumerate(exons):
        if exon_id in exon.id or exon_id in exon.attributes.get('exonic_part_number', []):
            found = True
            print(f"\nExon {exon_id} details:")
            print(f"Position: {i+1} out of {len(exons)} exons")
            print(f"Coordinates: {exon.seqid}:{exon.start}-{exon.end}")
            print(f"Length: {exon.end - exon.start} bp")
            print(f"Transcripts: {exon.attributes.get('transcripts', ['None'])}")
            print("\nSurrounding exons:")
            if i > 0:
                prev_exon = exons[i-1]
                print(f"Previous exon: {prev_exon.id}")
                print(f"Previous coordinates: {prev_exon.seqid}:{prev_exon.start}-{prev_exon.end}")
            if i < len(exons)-1:
                next_exon = exons[i+1]
                print(f"Next exon: {next_exon.id}")
                print(f"Next coordinates: {next_exon.seqid}:{next_exon.start}-{next_exon.end}")
            break
    
    if not found:
        print(f"\nNo exon found with ID {exon_id}")

def print_all_exons(db, gene_id="ENSG00000114744.9"):
    """
    Print details of all exons in the gene
    """
    # Get all exonic parts
    exons = []
    for feature in db.features_of_type('exonic_part', order_by='start'):
        if gene_id in feature.attributes.get('gene_id', []):
            exons.append(feature)
    
    # Sort exons by position
    exons.sort(key=lambda x: x.start)
    
    print(f"\nCOMM2 Gene Exon Details:")
    print(f"Total number of exons: {len(exons)}")
    print("\nExon coordinates:")
    print(f"{'Exon ID':<10} {'Coordinates':<30} {'Length':<10} {'Transcripts'}")
    print("-" * 80)
    
    for exon in exons:
        exon_number = exon.attributes.get('exonic_part_number', ['???'])[0]
        coordinates = f"{exon.seqid}:{exon.start}-{exon.end}"
        length = exon.end - exon.start
        transcripts = ','.join(exon.attributes.get('transcripts', ['None']))
        
        print(f"E{exon_number:<9} {coordinates:<30} {length:<10} {transcripts}")

def extract_exon_coordinates(db, gene_id):
    """
    Extract coordinates of all exons for a given gene in a tuple format.
    
    Args:
        db: gffutils database object
        gene_id (str): Gene ID (e.g., "ENSG00000114744.9")
        
    Returns:
        list: List of tuples containing (start, end, exon_id) coordinates for each exon
    """
    # Get all exonic parts
    exons = []
    for feature in db.features_of_type('exonic_part', order_by='start'):
        if gene_id in feature.attributes.get('gene_id', []):
            exon_id = f"E{feature.attributes.get('exonic_part_number', ['???'])[0]}"
            exons.append((int(feature.start), int(feature.end), exon_id))
    
    # Sort exons by position
    exons.sort(key=lambda x: x[0])
    
    return exons

def get_ensembl_exon_info(gene_id="ENSG00000114744", assembly="GRCh38", dexseq_parts=None):
    """
    Get canonical exon information from Ensembl and show overlaps with DexSeq parts
    
    Args:
        gene_id (str): Ensembl gene ID
        assembly (str): Genome assembly version
        dexseq_parts (list): List of tuples containing (start, end, *_) coordinates
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{gene_id}?expand=1"
    
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    
    if r.ok:
        decoded = r.json()
        print("\nEnsembl Transcript Information:")
        print(f"Gene: {decoded.get('display_name')} ({gene_id})")
        
        # Get transcript details
        transcripts = decoded.get('Transcript', [])
        print(f"\nFound {len(transcripts)} transcripts")
        
        # Track all exons across all transcripts
        all_exons = []
        for transcript in transcripts:
            transcript_id = transcript['id']
            sorted_exons = sorted(transcript['Exon'], key=lambda x: x['start'])
            for i, exon in enumerate(sorted_exons, 1):
                all_exons.append({
                    'transcript_id': transcript_id,
                    'exon_number': i,
                    'start': exon['start'],
                    'end': exon['end']
                })
        
        # Print canonical transcript structure first
        canonical = next((t for t in transcripts if t['is_canonical']), None)
        if canonical:
            print("\nCanonical transcript structure:")
            print(f"Transcript ID: {canonical['id']}")
            print(f"{'Exon Number':<12} {'Coordinates':<30} {'Length':<10}")
            print("-" * 80)
            
            sorted_exons = sorted(canonical['Exon'], key=lambda x: x['start'])
            for i, exon in enumerate(sorted_exons, 1):
                coords = f"chr{decoded['seq_region_name']}:{exon['start']}-{exon['end']}"
                length = exon['end'] - exon['start'] + 1
                print(f"{i:<12} {coords:<30} {length:<10}")
        
        # Compare DexSeq parts with all exons
        print("\nDexSeq parts and their overlaps:")
        print(f"{'DexSeq Part':<15} {'Coordinates':<30} {'Overlapping Biological Exons'}")
        print("-" * 80)
        
        for j, part in enumerate(dexseq_parts, 1):
            start, end = part[0], part[1]
            overlapping_exons = []
            
            for exon in all_exons:
                if (start <= exon['end'] and end >= exon['start']):
                    overlap = {
                        'transcript_id': exon['transcript_id'],
                        'exon_number': exon['exon_number'],
                        'coords': f"{exon['start']}-{exon['end']}"
                    }
                    if overlap not in overlapping_exons:
                        overlapping_exons.append(overlap)
            
            # Format overlapping exons information
            if overlapping_exons:
                exon_info = "; ".join([
                    f"Exon {e['exon_number']} ({e['transcript_id']}, {e['coords']})"
                    for e in overlapping_exons
                ])
            else:
                exon_info = "No overlaps"
            
            dexseq_coords = f"chr{decoded['seq_region_name']}:{start}-{end}"
            print(f"E{j:03d}{' ':<11} {dexseq_coords:<30} {exon_info}")
    
    else:
        print(f"Failed to fetch data: {r.status_code}")

def get_canonical_exons(gene_id="ENSG00000114744", assembly="GRCh38"):
    """
    Get canonical transcript exons from Ensembl and return as list of tuples.
    
    Args:
        gene_id (str): Ensembl gene ID
        assembly (str): Genome assembly version (default: GRCh38)
        
    Returns:
        list: List of tuples (start, end, label) for each canonical exon
    """
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{gene_id}?expand=1"
    
    try:
        r = requests.get(server + ext, headers={"Content-Type": "application/json"})
        r.raise_for_status()
        
        decoded = r.json()
        transcripts = decoded.get('Transcript', [])
        
        # Find canonical transcript
        canonical = next((t for t in transcripts if t['is_canonical']), None)
        if not canonical:
            print("No canonical transcript found")
            return []
            
        # Extract and sort exons
        exons = canonical.get('Exon', [])
        sorted_exons = sorted(exons, key=lambda x: x['start'])
        
        # Format as list of tuples
        result = []
        for i, exon in enumerate(sorted_exons, 1):
            start = exon['start']
            end = exon['end']
            label = f"Biological Exon {i}"
            result.append((start, end, label))
            
        return result
        
    except requests.exceptions.RequestException as e:
        print(f"Error accessing Ensembl API: {e}")
        return []
    except Exception as e:
        print(f"Unexpected error: {e}")
        return []

def visualize_gene_structure(biological_exons, dexseq_parts, peak, gene_name="Gene", figsize=(15, 6)):
    """
    Visualize gene structure showing biological exons, DexSeq parts, and peak.
    
    Parameters:
    -----------
    biological_exons : list of tuples
        List of (start, end, label) tuples for biological exons
    dexseq_parts : list of tuples
        List of (start, end, label) tuples for DexSeq parts
    peak : tuple
        Tuple of (start, end, label) for the peak
    gene_name : str, optional
        Name of the gene for the plot title (default: "Gene")
    figsize : tuple, optional
        Figure size in inches (default: (15, 6))
    
    Returns:
    --------
    tuple
        (figure, axis) tuple containing the plot objects
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Plot biological exons
    for start, end, label in biological_exons:
        ax.broken_barh([(start, end-start)], (20, 4), 
                      color='blue', alpha=0.3, label=label)

    # Plot DexSeq parts
    for start, end, label in dexseq_parts:
        ax.broken_barh([(start, end-start)], (10, 4), 
                      color='green', alpha=0.5, label=label)

    # Plot peak
    ax.broken_barh([(peak[0], peak[1]-peak[0])], (30, 4), 
                  color='red', alpha=0.5, label=peak[2])

    # Customize plot
    ax.set_ylim(5, 40)
    ax.set_xlim(min(min(x[0] for x in biological_exons), min(x[0] for x in dexseq_parts))-100,
                max(max(x[1] for x in biological_exons), max(x[1] for x in dexseq_parts))+100)
    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Track')
    ax.grid(True)

    # Create custom legend
    biological_legend = plt.Rectangle((0,0),1,1, color='blue', alpha=0.3)
    dexseq_legend = plt.Rectangle((0,0),1,1, color='green', alpha=0.5)
    peak_legend = plt.Rectangle((0,0),1,1, color='red', alpha=0.5)

    ax.legend([biological_legend, dexseq_legend, peak_legend],
             ['Biological Exons', 'DexSeq Parts', 'Peak'],
             loc='center left', bbox_to_anchor=(1, 0.5))

    plt.title(f'{gene_name}: Biological Exons vs DexSeq Parts vs Peak')
    plt.tight_layout()

    # Calculate and print intersections
    intersections = []
    for start, end, label in biological_exons:
        if peak[0] <= end and peak[1] >= start:
            intersection_start = max(start, peak[0])
            intersection_end = min(end, peak[1])
            intersection_length = intersection_end - intersection_start
            intersections.append({
                'exon': label,
                'start': intersection_start,
                'end': intersection_end,
                'length': intersection_length
            })
            print(f"Peak overlaps with {label}: {intersection_start:,} - {intersection_end:,} ({intersection_length:,} bp)")

    return fig, ax, intersections

# Example usage:
# fig, ax, intersections = visualize_gene_structure(biological_exons, dexseq_parts, peak, gene_name="COMMD2")
# plt.show()

def visualize_gene_structure_with_specific_exon(biological_exons, dexseq_parts, peak, specific_exon=None, gene_name="Gene", figsize=(15, 6)):
    """
    Visualize gene structure showing biological exons, DexSeq parts, peak, and a specific exon.
    Also identifies and prints the closest biological exons to the specific exon.
    
    Parameters:
    -----------
    biological_exons : list of tuples
        List of (start, end, label) tuples for biological exons
    dexseq_parts : list of tuples
        List of (start, end, label) tuples for DexSeq parts
    peak : tuple
        Tuple of (start, end, label) for the peak
    specific_exon : tuple or None, optional
        Tuple of (start, end, label) for a specific exon to highlight
    gene_name : str, optional
        Name of the gene for the plot title (default: "Gene")
    figsize : tuple, optional
        Figure size in inches (default: (15, 6))
    
    Returns:
    --------
    tuple
        (figure, axis, intersections, closest_exons) tuple containing the plot objects, 
        intersection data, and information about closest biological exons
    """
    fig, ax = plt.subplots(figsize=figsize)

    # Plot biological exons
    for start, end, label in biological_exons:
        ax.broken_barh([(start, end-start)], (20, 4), 
                      color='blue', alpha=0.3, label=label)

    # Plot DexSeq parts
    for start, end, label in dexseq_parts:
        ax.broken_barh([(start, end-start)], (10, 4), 
                      color='green', alpha=0.5, label=label)

    # Plot peak
    ax.broken_barh([(peak[0], peak[1]-peak[0])], (30, 4), 
                  color='red', alpha=0.5, label=peak[2])

    # Find closest biological exons if specific_exon is provided
    closest_exons = {}
    if specific_exon:
        start, end, label = specific_exon
        specific_center = (start + end) / 2
        
        # Calculate distances to all biological exons
        distances = []
        for b_start, b_end, b_label in biological_exons:
            b_center = (b_start + b_end) / 2
            distance = abs(specific_center - b_center)
            distances.append({
                'label': b_label,
                'distance': distance,
                'start': b_start,
                'end': b_end,
                'center': b_center,
                'relative_position': 'upstream' if b_center < specific_center else 'downstream'
            })
        
        # Sort by distance
        distances.sort(key=lambda x: x['distance'])
        
        # Get closest upstream and downstream exons
        upstream_exons = [x for x in distances if x['relative_position'] == 'upstream']
        downstream_exons = [x for x in distances if x['relative_position'] == 'downstream']
        
        closest_exons = {
            'upstream': upstream_exons[0] if upstream_exons else None,
            'downstream': downstream_exons[0] if downstream_exons else None
        }
        
        # Print information about closest exons
        print(f"\nClosest biological exons to {label}:")
        if closest_exons['upstream']:
            up = closest_exons['upstream']
            distance_bp = abs(up['center'] - specific_center)
            print(f"Upstream: {up['label']}")
            print(f"  Distance: {distance_bp:,.0f} bp")
            print(f"  Coordinates: {up['start']:,} - {up['end']:,}")
            
        if closest_exons['downstream']:
            down = closest_exons['downstream']
            distance_bp = abs(down['center'] - specific_center)
            print(f"Downstream: {down['label']}")
            print(f"  Distance: {distance_bp:,.0f} bp")
            print(f"  Coordinates: {down['start']:,} - {down['end']:,}")

        # Plot specific exon
        ax.broken_barh([(start, end-start)], (0, 4), 
                      color='purple', alpha=0.7, label=label)
        # Add text label below the specific exon
        # ax.text(start, -2, label, fontsize=8, rotation=45, ha='right')

    # Customize plot
    ax.set_ylim(-3, 40)
    ax.set_xlim(min(min(x[0] for x in biological_exons), min(x[0] for x in dexseq_parts))-100,
                max(max(x[1] for x in biological_exons), max(x[1] for x in dexseq_parts))+100)
    ax.set_xlabel('Genomic Position')
    ax.set_ylabel('Track')
    ax.grid(True)

    # Create custom legend
    biological_legend = plt.Rectangle((0,0),1,1, color='blue', alpha=0.3)
    dexseq_legend = plt.Rectangle((0,0),1,1, color='green', alpha=0.5)
    peak_legend = plt.Rectangle((0,0),1,1, color='red', alpha=0.5)
    legend_elements = [biological_legend, dexseq_legend, peak_legend]
    legend_labels = ['Biological Exons', 'DexSeq Parts', 'Peak']
    
    if specific_exon:
        specific_legend = plt.Rectangle((0,0),1,1, color='purple', alpha=0.7)
        legend_elements.append(specific_legend)
        legend_labels.append('Specific Exon')

    ax.legend(legend_elements, legend_labels,
             loc='center left', bbox_to_anchor=(1, 0.5))

    plt.title(f'{gene_name}: Biological Exons vs DexSeq Parts vs Peak')
    plt.tight_layout()

    # Calculate and print intersections
    intersections = []
    for start, end, label in biological_exons:
        if peak[0] <= end and peak[1] >= start:
            intersection_start = max(start, peak[0])
            intersection_end = min(end, peak[1])
            intersection_length = intersection_end - intersection_start
            intersections.append({
                'exon': label,
                'start': intersection_start,
                'end': intersection_end,
                'length': intersection_length
            })
            print(f"\nPeak overlaps with {label}: {intersection_start:,} - {intersection_end:,} ({intersection_length:,} bp)")

    return fig, ax, intersections, closest_exons

def fetch_exon_coordinates(gene_id):
    """Fetch exon coordinates from Ensembl REST API"""
    base_url = "https://rest.ensembl.org"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    ext = f"/lookup/id/{gene_id}?expand=1"
    response = requests.get(base_url + ext, headers=headers)
    
    if response.status_code != 200:
        return []
    
    data = response.json()
    exon_coords = []
    
    # Collect unique exon coordinates from all transcripts
    for transcript in data.get('Transcript', []):
        for exon in transcript.get('Exon', []):
            coords = (exon['start'], exon['end'])
            if coords not in exon_coords:
                exon_coords.append(coords)
    
    return sorted(exon_coords)

def examine_gff_content(gff_file, gene_id, db):
    """
    Examine the content of the GFF file for a specific gene
    """    
    # First, let's look at all features in the database
    print("\nSample of features in the database:")
    for feature in db.all_features():
        print(f"\nFeature type: {feature.featuretype}")
        print(f"ID: {feature.id}")
        print(f"Coordinates: {feature.seqid}:{feature.start}-{feature.end}")
        print("Attributes:", dict(feature.attributes))
        break  # Just show the first one as an example
    
    # Try to find our specific exon
    print("\nTrying to find specific features...")
    try:
        # Try exact ID
        feature = db[f"{gene_id}:E063"]
        print("Found by exact ID!")
    except:
        print("Not found by exact ID")
    
    # Look at all features in the relevant chromosome region
    print("\nFeatures in the relevant region:")
    for feature in db.region(seqid='chr1', start=70147000, end=70149000):
        print(f"\nFeature type: {feature.featuretype}")
        print(f"ID: {feature.id}")
        print(f"Coordinates: {feature.seqid}:{feature.start}-{feature.end}")
        print("Attributes:", dict(feature.attributes))

def peek_gff_file(gff_file, num_lines=10):
    """
    Show the first few lines of the GFF file
    """
    with open(gff_file, 'r') as f:
        lines = []
        for line in f:
            if not line.startswith('#'):  # Skip comment lines
                lines.append(line.strip())
                if len(lines) >= num_lines:
                    break
    
    print("Sample lines from GFF file:")
    for line in lines:
        print(line)

def get_sequence_for_coordinates(chromosome, start, end):
    """Get sequence from Ensembl REST API for given coordinates"""
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}:{end}:1"
    
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    
    if r.status_code != 200:
        return None
        
    return r.json()['seq']

def get_gene_info(chromosome, start, end):
    """Get gene information from Ensembl REST API for given coordinates"""
    server = "https://rest.ensembl.org"
    ext = f"/overlap/region/human/{chromosome}:{start}:{end}?feature=gene"
    
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})
    
    if r.status_code != 200:
        return None
        
    return r.json()

def search_ensembl_region(sequence, chromosome, start, end, margin=5000):
    """
    Search for sequence matches in Ensembl transcripts within specified region
    
    Args:
        sequence (str): Query sequence
        chromosome (str): Chromosome number
        start (int): Start position
        end (int): End position
        margin (int): Additional base pairs to search on each side
    """
    # Ensembl REST API endpoint
    server = "https://rest.ensembl.org"
    
    # Get transcripts in the region
    ext = f"/overlap/region/human/{chromosome}:{start-margin}-{end+margin}?"
    params = {
        "feature": "transcript",
        "content-type": "application/json"
    }
    
    try:
        r = requests.get(server + ext, params=params)
        r.raise_for_status()  # Raise exception for bad status codes
        transcripts = r.json()
        
        if not transcripts:
            print(f"No transcripts found in region chr{chromosome}:{start-margin}-{end+margin}")
            return
            
        print(f"Found {len(transcripts)} transcripts in region. Checking for sequence matches...\n")
        
        # For each transcript, get sequence and check for matches
        for transcript in transcripts:
            transcript_id = transcript.get('transcript_id', transcript.get('id', 'Unknown ID'))
            
            # Get transcript sequence
            ext = f"/sequence/id/{transcript_id}?"
            params = {
                "content-type": "application/json",
                "type": "cdna"
            }
            
            r = requests.get(server + ext, params=params)
            r.raise_for_status()
            transcript_seq = r.json().get('seq', '')
            
            # Check for sequence match
            if sequence in transcript_seq:
                print(f"Match found in transcript: {transcript_id}")
                # Using .get() method with default values to avoid KeyError
                print(f"Gene name: {transcript.get('external_name', 'Unknown')}")
                print(f"Gene ID: {transcript.get('parent', 'Unknown')}")
                print(f"Transcript biotype: {transcript.get('biotype', 'Unknown')}")
                print(f"Transcript position: chr{chromosome}:{transcript.get('start', 'Unknown')}-{transcript.get('end', 'Unknown')}")
                print("-----")
            
    except requests.exceptions.RequestException as e:
        print(f"Error accessing Ensembl API: {e}")
        return
    except json.JSONDecodeError as e:
        print(f"Error parsing API response: {e}")
        return
    except Exception as e:
        print(f"Unexpected error: {e}")
        return

def search_ensembl_region2(sequence, chromosome, start, end, margin=5000):
    """
    Search for sequence matches in Ensembl GRCh38 transcripts within specified region
    
    Args:
        sequence (str): Query sequence
        chromosome (str): Chromosome number
        start (int): Start position
        end (int): End position
        margin (int): Additional base pairs to search on each side
    """
    # Ensembl GRCh38 REST API endpoint
    server = "https://grch37.rest.ensembl.org" if use_grch37 else "https://rest.ensembl.org"
    
    # Get transcripts in the region
    ext = f"/overlap/region/human/{chromosome}:{start-margin}-{end+margin}?"
    params = {
        "feature": "transcript",
        "content-type": "application/json",
        "assembly": "GRCh38.14"  # Explicitly specify GRCh38
    }
    
    try:
        r = requests.get(server + ext, params=params)
        r.raise_for_status()
        transcripts = r.json()
        
        if not transcripts:
            print(f"No transcripts found in region chr{chromosome}:{start-margin}-{end+margin} (GRCh38)")
            return
            
        print(f"Found {len(transcripts)} transcripts in GRCh38 region. Checking for sequence matches...\n")
        
        # For each transcript, get sequence and check for matches
        for transcript in transcripts:
            transcript_id = transcript.get('transcript_id', transcript.get('id', 'Unknown ID'))
            
            # Get transcript sequence
            ext = f"/sequence/id/{transcript_id}?"
            params = {
                "content-type": "application/json",
                "type": "cdna"
            }
            
            r = requests.get(server + ext, params=params)
            r.raise_for_status()
            transcript_seq = r.json().get('seq', '')
            
            # Check for sequence match (both forward and reverse complement)
            rev_comp = str(Seq(sequence).reverse_complement())
            if sequence in transcript_seq or rev_comp in transcript_seq:
                print(f"Match found in transcript: {transcript_id}")
                print(f"Gene name: {transcript.get('external_name', 'Unknown')}")
                print(f"Gene ID: {transcript.get('parent', 'Unknown')}")
                print(f"Transcript biotype: {transcript.get('biotype', 'Unknown')}")
                print(f"Transcript position (GRCh38): chr{chromosome}:{transcript.get('start', 'Unknown')}-{transcript.get('end', 'Unknown')}")
                print(f"Strand: {'+' if transcript.get('strand', 0) > 0 else '-'}")
                print("-----")
            
    except requests.exceptions.RequestException as e:
        print(f"Error accessing Ensembl API: {e}")
        return
    except json.JSONDecodeError as e:
        print(f"Error parsing API response: {e}")
        return
    except Exception as e:
        print(f"Unexpected error: {e}")
        return

# ##########################################################################################################################################################################################################

# def analyze_als_complementarity(row: pd.Series, dexseq_df: pd.DataFrame, 
#                               als1: str = 'CTTTTCCAAGGAATGTT', 
#                               als2: str = 'TTCCGATGAGAATGACG') -> dict:
#     """
#     Analyze ALS complementarity for an exon and its neighbors using local alignment.
#     """
#     print(f"Analyzing exon {row['exon_id']} from gene {row['dexseq_name'].split(':')[0]}")
    
#     # Get current exon info
#     current_gene = row['dexseq_name'].split(':')[0]
#     current_exon_num = int(row['exon_id'].replace('E', ''))
    
#     # Find neighboring exons from the same gene
#     gene_exons = dexseq_df[dexseq_df['dexseq_name'].str.startswith(current_gene)].copy()
#     gene_exons['exon_num'] = gene_exons['featureID'].str.replace('E', '').astype(int)
#     gene_exons = gene_exons.sort_values('exon_num')
    
#     # Get sequences for current and neighboring exons
#     sequences = {}
#     for idx, exon in gene_exons.iterrows():
#         if abs(exon['exon_num'] - current_exon_num) <= 1:  # Current and immediate neighbors
#             seq = get_exon_sequence(
#                 exon['genomicData.seqnames'],
#                 exon['genomicData.start'],
#                 exon['genomicData.end']
#             )
#             if seq:
#                 sequences[exon['featureID']] = seq.strip()

#     def align_sequence(query, target):
#         """Perform local alignment and return the best score."""
#         aligner = PairwiseAligner()
#         aligner.mode = 'local'
#         aligner.match_score = 1
#         aligner.mismatch_score = -1
#         aligner.open_gap_score = -100
#         aligner.extend_gap_score = -100
#         try:
#             alignments = aligner.align(query, target)
#             best_alignment = max(alignments, key=lambda a: a.score)
#             return {"score": best_alignment.score, "alignment": best_alignment}
#         except Exception as e:
#             print(f"Error during alignment: {str(e)}")
#             return {"score": 0, "alignment": None}
    
#     results = {}
#     for exon_id, seq in sequences.items():
#         try:
#             # Convert sequences to Seq objects - only forward sequence
#             forward_seq = Seq(seq.upper())
#             als1_seq = Seq(als1.upper())
#             als2_seq = Seq(als2.upper())
            
#             # Align with ALS1 - only forward alignment
#             forward_als1 = align_sequence(forward_seq, als1_seq)
#             best_als1 = forward_als1
            
#             # Align with ALS2 - only forward alignment
#             forward_als2 = align_sequence(forward_seq, als2_seq)
#             best_als2 = forward_als2
            
#             results[exon_id] = {
#                 'als1_score': best_als1["score"],
#                 'als2_score': best_als2["score"],
#                 'als1_alignment': str(best_als1["alignment"]) if best_als1["alignment"] else "",
#                 'als2_alignment': str(best_als2["alignment"]) if best_als2["alignment"] else "",
#                 'sequence': seq
#             }
            
#         except Exception as e:
#             print(f"Error processing sequence for exon {exon_id}: {str(e)}")
#             results[exon_id] = {
#                 'als1_score': 0,
#                 'als2_score': 0,
#                 'als1_alignment': "",
#                 'als2_alignment': "",
#                 'sequence': seq
#             }
    
#     return results

# def process_chunk(chunk_data):
#     """
#     Process a chunk of overlaps data in parallel
#     """
#     chunk, dexseq_df = chunk_data
#     print(f"Processing chunk of size {len(chunk)}")
#     results = []
    
#     for idx, row in chunk.iterrows():
#         try:
#             als_analysis = analyze_als_complementarity(row, dexseq_df)
            
#             # Prepare the result dictionary
#             row_dict = row.to_dict()
            
#             # Add ALS scores if available
#             if als_analysis and row['exon_id'] in als_analysis:
#                 result = als_analysis[row['exon_id']]
#                 row_dict.update({
#                     'current_exon_als1_score': result['als1_score'],
#                     'current_exon_als2_score': result['als2_score'],
#                     'neighboring_exons_analysis': als_analysis
#                 })
#             else:
#                 row_dict.update({
#                     'current_exon_als1_score': 0,
#                     'current_exon_als2_score': 0,
#                     'neighboring_exons_analysis': {}
#                 })
            
#             results.append(row_dict)
            
#         except Exception as e:
#             print(f"Error processing row with exon {row['exon_id']}: {str(e)}")
#             row_dict = row.to_dict()
#             row_dict.update({
#                 'current_exon_als1_score': 0,
#                 'current_exon_als2_score': 0,
#                 'neighboring_exons_analysis': {}
#             })
#             results.append(row_dict)
    
#     print(f"Completed chunk processing with {len(results)} results")
#     return results

# def analyze_overlaps_with_als(overlaps_df: pd.DataFrame, dexseq_df: pd.DataFrame) -> pd.DataFrame:
#     """
#     Analyze ALS complementarity for all overlapping exons and their neighbors using parallel processing.
#     """
#     print("Starting ALS complementarity analysis...")
    
#     # Determine number of cores to use (leave one free for system)
#     n_cores = multiprocessing.cpu_count() - 1
#     print(f"Using {n_cores} CPU cores for parallel processing")
    
#     # Split the dataframe into chunks
#     chunk_size = len(overlaps_df) // n_cores
#     if chunk_size == 0:
#         chunk_size = 1
#     chunks = [overlaps_df.iloc[i:i + chunk_size] for i in range(0, len(overlaps_df), chunk_size)]
#     print(f"Split data into {len(chunks)} chunks of approximately {chunk_size} rows each")
    
#     # Prepare data for parallel processing
#     chunk_data = [(chunk, dexseq_df) for chunk in chunks]
    
#     # Process chunks in parallel
#     all_results = []
#     print("Starting parallel processing...")
#     with ProcessPoolExecutor(max_workers=n_cores) as executor:
#         chunk_results = list(executor.map(process_chunk, chunk_data))
#         for results in chunk_results:
#             all_results.extend(results)
    
#     print(f"Analysis complete. Processed {len(all_results)} total entries")
#     return pd.DataFrame(all_results)

# ############################### Implementation 3 ##################################################################################################
# def analyze_als_complementarity(row: pd.Series, dexseq_df: pd.DataFrame, 
#                               als1: str = 'CTTTTCCAAGGAATGTT', 
#                               als2: str = 'TTCCGATGAGAATGACG') -> dict:
#     """
#     Analyze ALS complementarity for an exon and its neighbors using local alignment.
#     Added input validation and better sequence handling.
#     """
#     if not isinstance(row, pd.Series) or not isinstance(dexseq_df, pd.DataFrame):
#         print("Invalid input types")
#         return {}
        
#     try:
#         # print(f"Analyzing exon {row['exon_id']} from gene {row['dexseq_name'].split(':')[0]}")
        
#         # Get current exon info
#         current_gene = row['dexseq_name'].split(':')[0]
#         current_exon_num = int(row['exon_id'].replace('E', ''))
        
#         # Find neighboring exons from the same gene
#         gene_exons = dexseq_df[dexseq_df['dexseq_name'].str.startswith(current_gene)].copy()
#         if gene_exons.empty:
#             print(f"No exons found for gene {current_gene}")
#             return {}
            
#         print(f"Gene {current_gene} has {len(gene_exons)} exons")
#         print(gene_exons.head())

#         gene_exons['exon_num'] = gene_exons['featureID'].str.replace('E', '').astype(int)
#         gene_exons = gene_exons.sort_values('exon_num')
        
#         # Get sequences for current and neighboring exons
#         sequences = {}
#         for idx, exon in gene_exons.iterrows():
#             if abs(exon['exon_num'] - current_exon_num) <= 1:
#                 seq = get_exon_sequence(
#                     str(exon['genomicData.seqnames']),  # Ensure string type
#                     int(exon['genomicData.start']),     # Ensure int type
#                     int(exon['genomicData.end'])        # Ensure int type
#                 )
#                 if seq:
#                     sequences[exon['featureID']] = seq.strip().upper()  # Normalize sequence
        
#         if not sequences:
#             print(f"No valid sequences found for exon {row['exon_id']}")
#             return {}
        
#         results = {}
#         for exon_id, seq in sequences.items():
#             try:
#                 # Convert sequences to Seq objects with validation
#                 try:
#                     forward_seq = Seq(seq)
#                     als1_seq = Seq(als1.upper())
#                     als2_seq = Seq(als2.upper())
#                 except ValueError as e:
#                     print(f"Invalid sequence for exon {exon_id}: {str(e)}")
#                     continue
                
#                 # Align with ALS1
#                 forward_als1 = align_sequence(forward_seq, als1_seq)
#                 best_als1 = forward_als1
                
#                 # Align with ALS2
#                 forward_als2 = align_sequence(forward_seq, als2_seq)
#                 best_als2 = forward_als2
                
#                 results[exon_id] = {
#                     'als1_score': float(best_als1["score"]),  # Ensure float type
#                     'als2_score': float(best_als2["score"]),  # Ensure float type
#                     'als1_alignment': str(best_als1["alignment"]) if best_als1["alignment"] else "",
#                     'als2_alignment': str(best_als2["alignment"]) if best_als2["alignment"] else "",
#                     'sequence': seq,
#                     'sequence_length': len(seq)  # Add sequence length for validation
#                 }
                
#             except Exception as e:
#                 print(f"Error processing sequence for exon {exon_id}: {str(e)}")
#                 continue
                
#         return results
        
#     except Exception as e:
#         print(f"Error in analyze_als_complementarity: {str(e)}")
#         return {}

# def process_chunk(chunk_data: tuple) -> list:
#     """
#     Process a chunk of overlaps data in parallel with enhanced error handling.
#     """
#     try:
#         chunk, dexseq_df = chunk_data
#         print(f"Processing chunk of size {len(chunk)}")
#         results = []
        
#         for idx, row in chunk.iterrows():
#             try:
#                 als_analysis = analyze_als_complementarity(row, dexseq_df)
                
#                 # Prepare the result dictionary
#                 row_dict = row.to_dict()
                
#                 # Add ALS scores if available
#                 if als_analysis and row['exon_id'] in als_analysis:
#                     result = als_analysis[row['exon_id']]
#                     row_dict.update({
#                         'current_exon_als1_score': result['als1_score'],
#                         'current_exon_als2_score': result['als2_score'],
#                         'neighboring_exons_analysis': als_analysis
#                     })
#                 else:
#                     row_dict.update({
#                         'current_exon_als1_score': 0.0,  # Consistent float type
#                         'current_exon_als2_score': 0.0,  # Consistent float type
#                         'neighboring_exons_analysis': {}
#                     })
                
#                 results.append(row_dict)
                
#             except Exception as e:
#                 print(f"Error processing row with exon {row['exon_id']}: {str(e)}")
#                 row_dict = row.to_dict()
#                 row_dict.update({
#                     'current_exon_als1_score': 0.0,
#                     'current_exon_als2_score': 0.0,
#                     'neighboring_exons_analysis': {}
#                 })
#                 results.append(row_dict)
        
#         print(f"Completed chunk processing with {len(results)} results")
#         return results
        
#     except Exception as e:
#         print(f"Error in process_chunk: {str(e)}")
#         return []

# def analyze_overlaps_with_als(overlaps_df: pd.DataFrame, dexseq_df: pd.DataFrame) -> pd.DataFrame:
#     """
#     Analyze ALS complementarity for all overlapping exons and their neighbors using parallel processing.
#     Added input validation and better error handling.
#     """
#     if overlaps_df.empty or dexseq_df.empty:
#         print("Empty input DataFrame(s)")
#         return pd.DataFrame()
        
#     try:
#         print("Starting ALS complementarity analysis...")
        
#         # Determine number of cores to use (leave one free for system)
#         n_cores = max(1, multiprocessing.cpu_count() - 1)
#         print(f"Using {n_cores} CPU cores for parallel processing")
        
#         # Split the dataframe into chunks
#         chunk_size = max(1, len(overlaps_df) // n_cores)
#         chunks = [overlaps_df.iloc[i:i + chunk_size] for i in range(0, len(overlaps_df), chunk_size)]
#         print(f"Split data into {len(chunks)} chunks of approximately {chunk_size} rows each")
        
#         # Prepare data for parallel processing
#         chunk_data = [(chunk, dexseq_df) for chunk in chunks]
        
#         # Process chunks in parallel with timeout
#         all_results = []
#         print("Starting parallel processing...")
#         with ProcessPoolExecutor(max_workers=n_cores) as executor:
#             chunk_results = list(executor.map(process_chunk, chunk_data, timeout=3600))  # 1 hour timeout
#             for results in chunk_results:
#                 if results:  # Check if results is not empty
#                     all_results.extend(results)
        
#         if not all_results:
#             print("No results generated from analysis")
#             return pd.DataFrame()
            
#         print(f"Analysis complete. Processed {len(all_results)} total entries")
#         return pd.DataFrame(all_results)
        
#     except Exception as e:
#         print(f"Error in analyze_overlaps_with_als: {str(e)}")
#         return pd.DataFrame()
    

############################### Legacy ##################################################################################################
# def get_exon_sequence(chromosome: str, start: int, end: int) -> str:
#     """
#     Fetch genomic sequence for given coordinates using Ensembl REST API.
#     """
#     server = "https://rest.ensembl.org"
#     ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?"
    
#     r = requests.get(server + ext, headers={"Content-Type": "text/plain"})
#     if r.ok:
#         return r.text
#     return None


# def get_exon_sequence(chromosome: str, start: int, end: int) -> Optional[str]:
#     """
#     Fetch genomic sequence for given coordinates using Ensembl REST API.
#     Includes retry logic and better error handling.
#     """
#     server = "https://rest.ensembl.org"
#     ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?"
#     max_retries = 3
    
#     for attempt in range(max_retries):
#         try:
#             r = requests.get(server + ext, 
#                            headers={"Content-Type": "text/plain"},
#                            timeout=30)  # Add timeout
#             if r.ok:
#                 sequence = r.text.strip()
#                 if sequence:  # Check if sequence is not empty
#                     return sequence
#             time.sleep(1)  # Add delay between retries
#         except requests.exceptions.RequestException as e:
#             print(f"Request failed for {chromosome}:{start}-{end}, attempt {attempt + 1}: {str(e)}")
#             if attempt == max_retries - 1:
#                 print(f"Failed to fetch sequence after {max_retries} attempts")
#     return None