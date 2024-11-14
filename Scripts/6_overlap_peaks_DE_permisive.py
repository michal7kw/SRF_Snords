# Converted from 6_overlap_peaks_DE_permisive.ipynb

# # Environment

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

import importlib
import functions
importlib.reload(functions)
from functions import *

# Example usage
peaks_file = "./DATA/Peak.csv"
dexseq_file = "./output_v38/dexseq_results_PW1_vs_combined_controls.csv"
output_prefix = "./output_v38/overlap_analysis"

# # Load and process data

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

enriched_overlaps_df3 = analyze_overlaps_with_als_optimized(overlaps_df, dexseq_df)

# enriched_overlaps_df3 = analyze_overlaps_with_als(overlaps_df, dexseq_df)

# # Load the data back

output_path = "output_v38/enriched_overlaps_with_als.csv"
pickle_path = "output_v38/enriched_overlaps_with_als.pkl"

%%script false --no-raise-error
# Save as pickle (better for preserving data types and complex objects)
enriched_overlaps_df3.to_pickle(pickle_path)

# Save as CSV
enriched_overlaps_df3.to_csv(output_path, index=False)

# Load the data back

# From CSV
# loaded_df_csv = pd.read_csv(output_path)

# From pickle
enriched_overlaps_df3 = pd.read_pickle(pickle_path)

enriched_overlaps_df3.shape

enriched_overlaps_df3.head()

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

enriched_overlaps_df3_backup = enriched_overlaps_df3.copy()

%%script false --no-raise-error
enriched_overlaps_df3 = enriched_overlaps_df3_backup.copy()

enriched_overlaps_df3.loc[enriched_overlaps_df3["current_exon_als1_score"].isna(), "current_exon_als1_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["current_exon_als2_score"].isna(), "current_exon_als2_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["next_exon_als1_score"].isna(), "next_exon_als1_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["next_exon_als2_score"].isna(), "next_exon_als2_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["previous_exon_als1_score"].isna(), "previous_exon_als1_score"] = 0
enriched_overlaps_df3.loc[enriched_overlaps_df3["previous_exon_als2_score"].isna(), "previous_exon_als2_score"] = 0

enriched_overlaps_df3["sum_als1_score"] = enriched_overlaps_df3["current_exon_als1_score"] + enriched_overlaps_df3["next_exon_als1_score"] + enriched_overlaps_df3["previous_exon_als1_score"]
enriched_overlaps_df3["sum_als2_score"] = enriched_overlaps_df3["current_exon_als2_score"] + enriched_overlaps_df3["next_exon_als2_score"] + enriched_overlaps_df3["previous_exon_als2_score"]


enriched_overlaps_df3 = enriched_overlaps_df3[(enriched_overlaps_df3["sum_als1_score"] > 0) | (enriched_overlaps_df3["sum_als2_score"] > 0)]

enriched_overlaps_df3.shape

enriched_overlaps_df3.sort_values(by="sum_als1_score", ascending=False, inplace=True)

enriched_overlaps_df3.shape

enriched_overlaps_df3.head()

enriched_overlaps_df3.shape

enriched_overlaps_df3.to_excel('enriched_overlaps_df3_with_sum_scores.xlsx', index=False)

enriched_overlaps_df3 = pd.read_excel('enriched_overlaps_df3_with_sum_scores.xlsx')

pd.set_option('display.max_columns', None)
enriched_overlaps_df3[(enriched_overlaps_df3.gene_name == 'LRRC7') & (enriched_overlaps_df3.exon_id == 'E063')]


# Original peak data
#
# chr1:70148486-70148679,chr1,70148486,70148679,194,"Exon (ENST00000370952.4/55631, exon 14 of 15)",1,70107750,70148005,40256,57554,ENST00000588515.5,40736,ENSG00000033122,LRRC7,leucine rich repeat containing 7,,,controlledV5_peak_43,1,15,14,7.5,2.906890596
#

import matplotlib.pyplot as plt
import numpy as np

# Define ranges
range1 = (69665487, 70249510) # LRRC7 gene coordinates
range2 = (70026415, 70026802) # peak coordinates

fig, ax = plt.subplots(figsize=(12, 4))

# Plot ranges as horizontal bars
ax.broken_barh([(range1[0], range1[1]-range1[0])], (10, 9), color='blue', alpha=0.5, label='LRRC7 gene')
ax.broken_barh([(range2[0], range2[1]-range2[0])], (20, 9), color='green', alpha=0.5, label='peak')

# Calculate and plot intersection
intersection_start = max(range1[0], range2[0])
intersection_end = min(range1[1], range2[1])
ax.broken_barh([(intersection_start, intersection_end-intersection_start)], (15, 9), 
               color='red', alpha=0.5, label='Intersection')

# Customize plot
ax.set_ylim(5, 35)
ax.set_xlim(range1[0]-10000, range1[1]+10000)
ax.set_xlabel('Genomic Position')
ax.grid(True)
ax.legend()

# Add range labels
for i, (start, end) in enumerate([(range1[0], range1[1]), (range2[0], range2[1])]):
    y = 10 if i == 0 else 20
    plt.text(start, y-3, f'{start:,}', verticalalignment='top')
    plt.text(end, y-3, f'{end:,}', verticalalignment='top')

plt.title('LRRC7 Gene and Peak Intersection')
plt.tight_layout()
plt.show()

# Calculate intersection size
intersection_size = intersection_end - intersection_start
print(f"Intersection range: {intersection_start:,} - {intersection_end:,} ({intersection_size:,} positions)")

import requests
import json

def get_transcript_info(gene_id):
    """
    Fetch and display information about all transcripts for a gene
    """
    base_url = "https://rest.ensembl.org"
    headers = {
        "Content-Type": "application/json",
        "Accept": "application/json"
    }
    
    ext = f"/lookup/id/{gene_id}?expand=1"
    response = requests.get(base_url + ext, headers=headers)
    
    if response.status_code != 200:
        raise Exception(f"Failed to fetch gene data: {response.status_code}")
    
    data = response.json()
    transcripts = data.get('Transcript', [])
    
    print(f"Gene: {gene_id}")
    print(f"Number of transcripts: {len(transcripts)}")
    print("\nTranscript details:")
    
    for transcript in transcripts:
        exon_count = len(transcript.get('Exon', []))
        print(f"\nTranscript ID: {transcript['id']}")
        print(f"Exon count: {exon_count}")
        if exon_count > 0:
            first_exon = transcript['Exon'][0]
            last_exon = transcript['Exon'][-1]
            print(f"First exon coordinates: {first_exon['start']}-{first_exon['end']}")
            print(f"Last exon coordinates: {last_exon['start']}-{last_exon['end']}")

# Get information about all transcripts
get_transcript_info("ENSG00000033122")

import matplotlib.pyplot as plt
import numpy as np
import requests
import json

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

# Define ranges
range1 = (69665487, 70249510)  # LRRC7 gene coordinates
range2 = (70026415, 70026802)  # peak coordinates

# Fetch exon coordinates
exon_coords = fetch_exon_coordinates("ENSG00000033122")

fig, ax = plt.subplots(figsize=(12, 6))

# Plot gene range
ax.broken_barh([(range1[0], range1[1]-range1[0])], (10, 9), 
               color='blue', alpha=0.3, label='LRRC7 gene')

# Plot exons
first_exon = True
for start, end in exon_coords:
    ax.broken_barh([(start, end-start)], (10, 9), 
                  color='green', alpha=0.7, 
                  label='exon' if first_exon else None)
    first_exon = False

# Plot peak
ax.broken_barh([(range2[0], range2[1]-range2[0])], (20, 9), 
               color='red', alpha=0.5, label='peak')

# Customize plot
ax.set_ylim(5, 35)
ax.set_xlim(range1[0]-10000, range1[1]+10000)
ax.set_xlabel('Genomic Position')
ax.grid(True)
ax.legend()

# Add range labels
for i, (start, end) in enumerate([(range1[0], range1[1]), (range2[0], range2[1])]):
    y = 10 if i == 0 else 20
    plt.text(start, y-3, f'{start:,}', verticalalignment='top')
    plt.text(end, y-3, f'{end:,}', verticalalignment='top')

plt.title('LRRC7 Gene, Exons, and Peak Intersection')
plt.tight_layout()
plt.show()

# Calculate intersection size
intersection_size = intersection_end - intersection_start
print(f"Intersection range: {intersection_start:,} - {intersection_end:,} ({intersection_size:,} positions)")
print(f"Number of exons: {len(exon_coords)}")

visualize_top_alignments(enriched_overlaps_df3)

