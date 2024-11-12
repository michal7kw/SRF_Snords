# Converted from 6_overlap_peaks_DE_permisive.ipynb

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
dexseq_file = "./output/dexseq_results_PW1_vs_combined_controls_cleaned_permisive.csv"
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

dexseq_df.head()

# Create visualization
fig = visualize_distance_analysis(overlaps_df)
plt.show()

# # Complementarity analysis

ASE1 = 'AACATTCCTTGGAAAAG'
ASE2 = 'CGTCATTCTCATCGGAA'

cASE1 = 'CTTTTCCAAGGAATGTT'
cASE2 = 'TTCCGATGAGAATGACG'

%%script false --no-raise-error
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import requests
from concurrent.futures import ProcessPoolExecutor
import multiprocessing

def get_exon_sequence(chromosome: str, start: int, end: int) -> str:
    """
    Fetch genomic sequence for given coordinates using Ensembl REST API.
    """
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?"
    
    r = requests.get(server + ext, headers={"Content-Type": "text/plain"})
    if r.ok:
        return r.text
    return None

def analyze_als_complementarity(row: pd.Series, dexseq_df: pd.DataFrame, 
                              als1: str = 'CTTTTCCAAGGAATGTT', 
                              als2: str = 'TTCCGATGAGAATGACG') -> dict:
    """
    Analyze ALS complementarity for an exon and its neighbors using local alignment.
    """
    print(f"Analyzing exon {row['exon_id']} from gene {row['dexseq_name'].split(':')[0]}")
    
    # Get current exon info
    current_gene = row['dexseq_name'].split(':')[0]
    current_exon_num = int(row['exon_id'].replace('E', ''))
    
    # Find neighboring exons from the same gene
    gene_exons = dexseq_df[dexseq_df['dexseq_name'].str.startswith(current_gene)].copy()
    gene_exons['exon_num'] = gene_exons['featureID'].str.replace('E', '').astype(int)
    gene_exons = gene_exons.sort_values('exon_num')
    
    # Get sequences for current and neighboring exons
    sequences = {}
    for idx, exon in gene_exons.iterrows():
        if abs(exon['exon_num'] - current_exon_num) <= 1:  # Current and immediate neighbors
            seq = get_exon_sequence(
                exon['genomicData.seqnames'],
                exon['genomicData.start'],
                exon['genomicData.end']
            )
            if seq:
                sequences[exon['featureID']] = seq.strip()

    def align_sequence(query, target):
        """Perform local alignment and return the best score."""
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.open_gap_score = -100
        aligner.extend_gap_score = -100
        try:
            alignments = aligner.align(query, target)
            best_alignment = max(alignments, key=lambda a: a.score)
            return {"score": best_alignment.score, "alignment": best_alignment}
        except Exception as e:
            print(f"Error during alignment: {str(e)}")
            return {"score": 0, "alignment": None}
    
    results = {}
    for exon_id, seq in sequences.items():
        try:
            # Convert sequences to Seq objects - only forward sequence
            forward_seq = Seq(seq.upper())
            als1_seq = Seq(als1.upper())
            als2_seq = Seq(als2.upper())
            
            # Align with ALS1 - only forward alignment
            forward_als1 = align_sequence(forward_seq, als1_seq)
            best_als1 = forward_als1
            
            # Align with ALS2 - only forward alignment
            forward_als2 = align_sequence(forward_seq, als2_seq)
            best_als2 = forward_als2
            
            results[exon_id] = {
                'als1_score': best_als1["score"],
                'als2_score': best_als2["score"],
                'als1_alignment': str(best_als1["alignment"]) if best_als1["alignment"] else "",
                'als2_alignment': str(best_als2["alignment"]) if best_als2["alignment"] else "",
                'sequence': seq
            }
            
        except Exception as e:
            print(f"Error processing sequence for exon {exon_id}: {str(e)}")
            results[exon_id] = {
                'als1_score': 0,
                'als2_score': 0,
                'als1_alignment': "",
                'als2_alignment': "",
                'sequence': seq
            }
    
    return results

def process_chunk(chunk_data):
    """
    Process a chunk of overlaps data in parallel
    """
    chunk, dexseq_df = chunk_data
    print(f"Processing chunk of size {len(chunk)}")
    results = []
    
    for idx, row in chunk.iterrows():
        try:
            als_analysis = analyze_als_complementarity(row, dexseq_df)
            
            # Prepare the result dictionary
            row_dict = row.to_dict()
            
            # Add ALS scores if available
            if als_analysis and row['exon_id'] in als_analysis:
                result = als_analysis[row['exon_id']]
                row_dict.update({
                    'current_exon_als1_score': result['als1_score'],
                    'current_exon_als2_score': result['als2_score'],
                    'neighboring_exons_analysis': als_analysis
                })
            else:
                row_dict.update({
                    'current_exon_als1_score': 0,
                    'current_exon_als2_score': 0,
                    'neighboring_exons_analysis': {}
                })
            
            results.append(row_dict)
            
        except Exception as e:
            print(f"Error processing row with exon {row['exon_id']}: {str(e)}")
            row_dict = row.to_dict()
            row_dict.update({
                'current_exon_als1_score': 0,
                'current_exon_als2_score': 0,
                'neighboring_exons_analysis': {}
            })
            results.append(row_dict)
    
    print(f"Completed chunk processing with {len(results)} results")
    return results


def analyze_overlaps_with_als(overlaps_df: pd.DataFrame, dexseq_df: pd.DataFrame) -> pd.DataFrame:
    """
    Analyze ALS complementarity for all overlapping exons and their neighbors using parallel processing.
    """
    print("Starting ALS complementarity analysis...")
    
    # Determine number of cores to use (leave one free for system)
    n_cores = multiprocessing.cpu_count() - 1
    print(f"Using {n_cores} CPU cores for parallel processing")
    
    # Split the dataframe into chunks
    chunk_size = len(overlaps_df) // n_cores
    if chunk_size == 0:
        chunk_size = 1
    chunks = [overlaps_df.iloc[i:i + chunk_size] for i in range(0, len(overlaps_df), chunk_size)]
    print(f"Split data into {len(chunks)} chunks of approximately {chunk_size} rows each")
    
    # Prepare data for parallel processing
    chunk_data = [(chunk, dexseq_df) for chunk in chunks]
    
    # Process chunks in parallel
    all_results = []
    print("Starting parallel processing...")
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        chunk_results = list(executor.map(process_chunk, chunk_data))
        for results in chunk_results:
            all_results.extend(results)
    
    print(f"Analysis complete. Processed {len(all_results)} total entries")
    return pd.DataFrame(all_results)

# After running the analysis:
enriched_overlaps_df = analyze_overlaps_with_als(overlaps_df, dexseq_df)

%%script false --no-raise-error
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
import requests
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import pandas as pd
from typing import Dict, Optional, Any
import time

def get_exon_sequence(chromosome: str, start: int, end: int) -> Optional[str]:
    """
    Fetch genomic sequence for given coordinates using Ensembl REST API.
    Includes retry logic and better error handling.
    """
    server = "https://rest.ensembl.org"
    ext = f"/sequence/region/human/{chromosome}:{start}..{end}:1?"
    max_retries = 3
    
    for attempt in range(max_retries):
        try:
            r = requests.get(server + ext, 
                           headers={"Content-Type": "text/plain"},
                           timeout=30)  # Add timeout
            if r.ok:
                sequence = r.text.strip()
                if sequence:  # Check if sequence is not empty
                    return sequence
            time.sleep(1)  # Add delay between retries
        except requests.exceptions.RequestException as e:
            print(f"Request failed for {chromosome}:{start}-{end}, attempt {attempt + 1}: {str(e)}")
            if attempt == max_retries - 1:
                print(f"Failed to fetch sequence after {max_retries} attempts")
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
    Added input validation and better sequence handling.
    """
    if not isinstance(row, pd.Series) or not isinstance(dexseq_df, pd.DataFrame):
        print("Invalid input types")
        return {}
        
    try:
        # print(f"Analyzing exon {row['exon_id']} from gene {row['dexseq_name'].split(':')[0]}")
        
        # Get current exon info
        current_gene = row['dexseq_name'].split(':')[0]
        current_exon_num = int(row['exon_id'].replace('E', ''))
        
        # Find neighboring exons from the same gene
        gene_exons = dexseq_df[dexseq_df['dexseq_name'].str.startswith(current_gene)].copy()
        if gene_exons.empty:
            print(f"No exons found for gene {current_gene}")
            return {}
            
        print(f"Gene {current_gene} has {len(gene_exons)} exons")
        print(gene_exons.head())

        gene_exons['exon_num'] = gene_exons['featureID'].str.replace('E', '').astype(int)
        gene_exons = gene_exons.sort_values('exon_num')
        
        # Get sequences for current and neighboring exons
        sequences = {}
        for idx, exon in gene_exons.iterrows():
            if abs(exon['exon_num'] - current_exon_num) <= 1:
                seq = get_exon_sequence(
                    str(exon['genomicData.seqnames']),  # Ensure string type
                    int(exon['genomicData.start']),     # Ensure int type
                    int(exon['genomicData.end'])        # Ensure int type
                )
                if seq:
                    sequences[exon['featureID']] = seq.strip().upper()  # Normalize sequence
        
        if not sequences:
            print(f"No valid sequences found for exon {row['exon_id']}")
            return {}
        
        results = {}
        for exon_id, seq in sequences.items():
            try:
                # Convert sequences to Seq objects with validation
                try:
                    forward_seq = Seq(seq)
                    als1_seq = Seq(als1.upper())
                    als2_seq = Seq(als2.upper())
                except ValueError as e:
                    print(f"Invalid sequence for exon {exon_id}: {str(e)}")
                    continue
                
                # Align with ALS1
                forward_als1 = align_sequence(forward_seq, als1_seq)
                best_als1 = forward_als1
                
                # Align with ALS2
                forward_als2 = align_sequence(forward_seq, als2_seq)
                best_als2 = forward_als2
                
                results[exon_id] = {
                    'als1_score': float(best_als1["score"]),  # Ensure float type
                    'als2_score': float(best_als2["score"]),  # Ensure float type
                    'als1_alignment': str(best_als1["alignment"]) if best_als1["alignment"] else "",
                    'als2_alignment': str(best_als2["alignment"]) if best_als2["alignment"] else "",
                    'sequence': seq,
                    'sequence_length': len(seq)  # Add sequence length for validation
                }
                
            except Exception as e:
                print(f"Error processing sequence for exon {exon_id}: {str(e)}")
                continue
                
        return results
        
    except Exception as e:
        print(f"Error in analyze_als_complementarity: {str(e)}")
        return {}

def process_chunk(chunk_data: tuple) -> list:
    """
    Process a chunk of overlaps data in parallel with enhanced error handling.
    """
    try:
        chunk, dexseq_df = chunk_data
        print(f"Processing chunk of size {len(chunk)}")
        results = []
        
        for idx, row in chunk.iterrows():
            try:
                als_analysis = analyze_als_complementarity(row, dexseq_df)
                
                # Prepare the result dictionary
                row_dict = row.to_dict()
                
                # Add ALS scores if available
                if als_analysis and row['exon_id'] in als_analysis:
                    result = als_analysis[row['exon_id']]
                    row_dict.update({
                        'current_exon_als1_score': result['als1_score'],
                        'current_exon_als2_score': result['als2_score'],
                        'neighboring_exons_analysis': als_analysis
                    })
                else:
                    row_dict.update({
                        'current_exon_als1_score': 0.0,  # Consistent float type
                        'current_exon_als2_score': 0.0,  # Consistent float type
                        'neighboring_exons_analysis': {}
                    })
                
                results.append(row_dict)
                
            except Exception as e:
                print(f"Error processing row with exon {row['exon_id']}: {str(e)}")
                row_dict = row.to_dict()
                row_dict.update({
                    'current_exon_als1_score': 0.0,
                    'current_exon_als2_score': 0.0,
                    'neighboring_exons_analysis': {}
                })
                results.append(row_dict)
        
        print(f"Completed chunk processing with {len(results)} results")
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
    
enriched_overlaps_df2 = analyze_overlaps_with_als(overlaps_df, dexseq_df)

%%script false --no-raise-error
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

enriched_overlaps_df3 = analyze_overlaps_with_als(overlaps_df, dexseq_df)

# # Load the data back

# Save the DataFrame
output_path = "enriched_overlaps_with_als.csv"
pickle_path = "enriched_overlaps_with_als.pkl"

# # Save as CSV
# enriched_overlaps_df3.to_csv(output_path, index=False)

# # Save as pickle (better for preserving data types and complex objects)
# enriched_overlaps_df3.to_pickle(pickle_path)

# Load the data back
# From CSV
# loaded_df_csv = pd.read_csv(output_path)

# From pickle
enriched_overlaps_df3 = pd.read_pickle(pickle_path)

# Verify the data loaded correctly
# print(f"Original DataFrame shape: {enriched_overlaps_df3.shape}")
# print(f"Loaded CSV DataFrame shape: {loaded_df_csv.shape}")
# print(f"Loaded pickle DataFrame shape: {loaded_df_pickle.shape}")

# # Check if all columns were preserved
# print("\nColumns in original:", enriched_overlaps_df3.columns.tolist())
# print("Columns in loaded pickle:", loaded_df_pickle.columns.tolist())

# # Quick comparison of data types
# print("\nData types comparison:")
# print("Original:\n", enriched_overlaps_df3.dtypes)
# print("\nLoaded from pickle:\n", loaded_df_pickle.dtypes)

# Load the enriched overlaps dataframe from CSV
enriched_overlaps_df = pd.read_csv('enriched_overlaps_results.csv')

# Convert the string representation of dictionary back to dictionary
enriched_overlaps_df['neighboring_exons_analysis'] = enriched_overlaps_df['neighboring_exons_analysis'].apply(eval)

print("Loaded enriched overlaps results from enriched_overlaps_results.csv")

# Load the enriched overlaps dataframe from CSV
enriched_overlaps_df2 = pd.read_csv('enriched_overlaps_results2.csv')

# Convert the string representation of dictionary back to dictionary
enriched_overlaps_df2['neighboring_exons_analysis'] = enriched_overlaps_df2['neighboring_exons_analysis'].apply(eval)

print("Loaded enriched overlaps results from enriched_overlaps_results2.csv")

# # Save the enriched overlaps dataframe to a CSV file
# enriched_overlaps_df.to_csv('enriched_overlaps_results.csv', index=False)
# print("Saved enriched overlaps results to enriched_overlaps_results.csv")

# # Save the enriched overlaps dataframe to a CSV file
# enriched_overlaps_df2.to_csv('enriched_overlaps_results2.csv', index=False)
# print("Saved enriched overlaps results to enriched_overlaps_results2.csv")

enriched_overlaps_df.head()

enriched_overlaps_df2.head()

enriched_overlaps_df3.head()

print(enriched_overlaps_df3.shape)
print(enriched_overlaps_df2.shape)
print(enriched_overlaps_df.shape)

# Compare basic statistics
print("First version stats:")
print(enriched_overlaps_df['current_exon_als1_score'].describe())
print("\nSecond version stats:")
print(enriched_overlaps_df2['current_exon_als1_score'].describe())

# Check for differences in scores
diff_mask = (enriched_overlaps_df['current_exon_als1_score'] != 
             enriched_overlaps_df2['current_exon_als1_score'])
print("\nRows with different scores:", sum(diff_mask))

def visualize_top_alignments(enriched_overlaps_df: pd.DataFrame, n: int = 5):
    """
    Visualize the top n alignments for both ALS sequences.
    """
    # Sort by ALS1 score
    top_als1 = enriched_overlaps_df.nlargest(n, 'current_exon_als1_score')
    print("Top ALS1 alignments:")
    for _, row in top_als1.iterrows():
        print(f"\nGene: {row['dexseq_name'].split(':')[0]}, Exon: {row['exon_id']}")
        print(row['neighboring_exons_analysis'][row['exon_id']]['als1_alignment'])
    
    print("\n" + "="*50 + "\n")
    
    # Sort by ALS2 score
    top_als2 = enriched_overlaps_df.nlargest(n, 'current_exon_als2_score')
    print("Top ALS2 alignments:")
    for _, row in top_als2.iterrows():
        print(f"\nGene: {row['dexseq_name'].split(':')[0]}, Exon: {row['exon_id']}")
        print(row['neighboring_exons_analysis'][row['exon_id']]['als2_alignment'])

visualize_top_alignments(enriched_overlaps_df2)

# Visualize ALS scores distribution
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

sns.histplot(data=enriched_overlaps_df, x='current_exon_als1_score', ax=ax1)
ax1.set_title('Distribution of ALS1 Complementarity Scores')
ax1.set_xlabel('ALS1 Score')

sns.histplot(data=enriched_overlaps_df, x='current_exon_als2_score', ax=ax2)
ax2.set_title('Distribution of ALS2 Complementarity Scores')
ax2.set_xlabel('ALS2 Score')

plt.tight_layout()
plt.show()

# Visualize ALS scores distribution
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

sns.histplot(data=enriched_overlaps_df2, x='current_exon_als1_score', ax=ax1)
ax1.set_title('Distribution of ALS1 Complementarity Scores')
ax1.set_xlabel('ALS1 Score')

sns.histplot(data=enriched_overlaps_df2, x='current_exon_als2_score', ax=ax2)
ax2.set_title('Distribution of ALS2 Complementarity Scores')
ax2.set_xlabel('ALS2 Score')

plt.tight_layout()
plt.show()

# Visualize ALS scores distribution
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

sns.histplot(data=enriched_overlaps_df3, x='current_exon_als1_score', ax=ax1)
ax1.set_title('Distribution of ALS1 Complementarity Scores')
ax1.set_xlabel('ALS1 Score')

sns.histplot(data=enriched_overlaps_df3, x='next_exon_als1_score', ax=ax2)
ax2.set_title('Distribution of ALS2 Complementarity Scores')
ax2.set_xlabel('ALS2 Score')

plt.tight_layout()
plt.show()

enriched_overlaps_df3["groupID"] = enriched_overlaps_df3["dexseq_name"].str.split(':').str[0]

enriched_overlaps_df3.head()

enriched_overlaps_df3["geneID"] = enriched_overlaps_df3["groupID"].str.split('.').str[0]

enriched_overlaps_df3.head()

import mygene

# Initialize mygene client
mg = mygene.MyGeneInfo()

# Remove version numbers from ENSEMBL IDs
sel_genes_no_version = enriched_overlaps_df3["geneID"]

# Query the gene symbols
results = mg.querymany(sel_genes_no_version, scopes='ensembl.gene', fields='symbol', species='human')


# Add gene names to enriched_overlaps_df3
enriched_overlaps_df3["gene_name"] = [item.get('symbol', '') for item in results]

enriched_overlaps_df3.head()

enriched_overlaps_df3.loc[enriched_overlaps_df3["current_exon_als1_score"].isna(), "current_exon_als1_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["current_exon_als2_score"].isna(), "current_exon_als2_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["next_exon_als1_score"].isna(), "next_exon_als1_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["next_exon_als2_score"].isna(), "next_exon_als2_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["previous_exon_als1_score"].isna(), "previous_exon_als1_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["previous_exon_als2_score"].isna(), "previous_exon_als2_score"] = 0

enriched_overlaps_df3["sum_als1_score"] = enriched_overlaps_df3["current_exon_als1_score"] + enriched_overlaps_df3["next_exon_als1_score"] + enriched_overlaps_df3["previous_exon_als1_score"]
enriched_overlaps_df3["sum_als2_score"] = enriched_overlaps_df3["current_exon_als2_score"] + enriched_overlaps_df3["next_exon_als2_score"] + enriched_overlaps_df3["previous_exon_als2_score"]


enriched_overlaps_df3.to_excel('enriched_overlaps_df3_with_sum_scores.xlsx', index=False)

enriched_overlaps_df2["groupID"] = enriched_overlaps_df2["dexseq_name"].str.split(':').str[0]

enriched_overlaps_df2 = enriched_overlaps_df2.sort_values('current_exon_als1_score', ascending=False)

enriched_overlaps_df2.head()

enriched_overlaps_df2.to_excel('enriched_overlaps_df2.xlsx', index=False)

n_sel = 100 
sel_genes = enriched_overlaps_df2["groupID"][:n_sel]

import mygene

# Initialize mygene client
mg = mygene.MyGeneInfo()

# Remove version numbers from ENSEMBL IDs
sel_genes_no_version = [id.split('.')[0] for id in sel_genes]

# Query the gene symbols
results = mg.querymany(sel_genes_no_version, scopes='ensembl.gene', fields='symbol', species='human')

# Create a dictionary mapping ENSEMBL IDs to symbols
gene_map = {res['query']: res.get('symbol', 'Not found') for res in results}

gene_symbols = [gene_map[ensembl_id.split('.')[0]] for ensembl_id in sel_genes]

ensembl_to_symbol = dict(zip(sel_genes, gene_symbols))

# # Print results
# for ensembl, symbol in ensembl_to_symbol.items():
#     print(f"{ensembl}: {symbol}")

print(gene_symbols)

