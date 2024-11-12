# Standard library imports
import os
import signal
import time
from typing import Any, Dict, List, Optional, Tuple

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
    
############################### Implementation 1 ##################################################################################################
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

############################### Implementation 1 Optimized #########################################################################################
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


# ############################### Implementation 2 ##################################################################################################
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