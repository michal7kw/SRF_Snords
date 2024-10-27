#%% Imports and Setup
# Standard library imports
import os
from typing import Set, List, Optional, Dict, Tuple

# Third-party imports
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector

# Enable automatic conversion between pandas and R dataframes
pandas2ri.activate()

# Import R packages
dexseq = importr('DEXSeq')

#%% Function to extract GFF features
def get_gff_features(dexseq_gff: str) -> Set[str]:
    """Extract features from DEXSeq-formatted GFF file."""
    if not os.path.exists(dexseq_gff):
        raise FileNotFoundError(f"DEXSeq annotation file not found: {dexseq_gff}")
        
    gff_features = set()
    
    with open(dexseq_gff, 'r') as f:
        for line in f:
            if 'exonic_part' in line:
                fields = line.strip().split('\t')
                attrs = dict(item.strip().split(' ', 1) for item in fields[8].strip().split(';'))
                gene_id = attrs['gene_id'].strip('"')
                exon_num = attrs['exonic_part_number'].strip('"')
                feature_id = f"{gene_id}:E{exon_num}"
                gff_features.add(feature_id)
    
    print(f"Loaded {len(gff_features)} features from GFF file")
    print("Sample features:", list(gff_features)[:5])
    return gff_features

#%% Function to process count file
def process_count_file(file_path: str, output_dir: str, gff_features: Set[str]) -> Optional[str]:
    """Process DEXSeq count file with debugging."""
    os.makedirs(output_dir, exist_ok=True)
    basename = os.path.basename(file_path)
    output_path = os.path.join(output_dir, f"processed_{basename}")
    
    print(f"\nProcessing file: {basename}")
    
    try:
        # Debug: Check file content
        print(f"Reading file: {file_path}")
        with open(file_path, 'r') as f:
            first_lines = [next(f) for _ in range(5)]
        print("First 5 lines of input file:")
        for line in first_lines:
            print(line.strip())

        count_dict = {}
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if not line.startswith('_'):
                    try:
                        feature_id, count = line.strip().split('\t')
                        feature_id = feature_id.strip('"')
                        
                        if '":"' in feature_id:
                            gene_id, exon_num = feature_id.split('":"')
                            gene_id = gene_id.strip('"')
                            exon_num = exon_num.strip('"')
                            feature_id = f"{gene_id}:E{exon_num}"
                            
                            if feature_id in gff_features:
                                count_dict[feature_id] = int(count)
                            else:
                                print(f"Feature not found in GFF: {feature_id}")
                    except Exception as e:
                        print(f"Error processing line {line_num}: {line.strip()}")
                        print(f"Error details: {str(e)}")
        
        print(f"Processed {len(count_dict)} features")
        
        with open(output_path, 'w') as f:
            for feature in sorted(gff_features):
                count = count_dict.get(feature, 0)
                f.write(f"{feature}\t{count}\n")
        
        if os.path.exists(output_path):
            print(f"Output file created: {output_path}")
            with open(output_path, 'r') as f:
                first_lines = [next(f) for _ in range(5)]
            print("First 5 lines of output file:")
            for line in first_lines:
                print(line.strip())
        
        verified = verify_processed_file(output_path, gff_features)
        print(f"File verification: {'Passed' if verified else 'Failed'}")
        
        return output_path if verified else None
        
    except Exception as e:
        print(f"Error processing {basename}: {str(e)}")
        return None

#%% Function to verify processed file
def verify_processed_file(file_path: str, gff_features: Set[str]) -> bool:
    """Verify that processed file matches GFF features exactly."""
    try:
        file_features = set()
        with open(file_path, 'r') as f:
            for line in f:
                feature_id = line.strip().split('\t')[0]
                file_features.add(feature_id)
        
        return file_features == gff_features
        
    except Exception as e:
        print(f"Error verifying {file_path}: {str(e)}")
        return False

#%% Function to create DEXSeq dataset
def create_dexseq_dataset(sample_info: pd.DataFrame, processed_files: List[str], dexseq_gff: str) -> ro.vectors.Vector:
    """Create DEXSeqDataSet with proper formatting."""
    try:
        sample_data = pd.DataFrame({
            'sample': sample_info['sample'],
            'condition': sample_info['condition']
        }).sort_values('sample')
        
        print("Sample data for DEXSeq:")
        print(sample_data)
        
        with ro.default_converter + pandas2ri.converter:
            sample_data_r = ro.conversion.py2rpy(sample_data)
        
        dxd = dexseq.DEXSeqDataSetFromHTSeq(
            countfiles=StrVector(processed_files),
            sampleData=sample_data_r,
            design=Formula('~ sample + exon + condition:exon'),
            flattenedfile=dexseq_gff
        )
        
        dxd = dexseq.estimateSizeFactors(dxd)
        dxd = dexseq.estimateDispersions(dxd)
        
        return dxd
        
    except Exception as e:
        print(f"Error creating DEXSeq dataset: {str(e)}")
        raise

#%% Configuration
CONFIG = {
    'dexseq_gff': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.dexseq.gff",
    'data_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output",
    'output_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/output"
}

#%% Test individual steps
# 1. Get GFF features
try:
    gff_features = get_gff_features(CONFIG['dexseq_gff'])
except Exception as e:
    print(f"Error getting GFF features: {e}")
    gff_features = None

#%% 2. Process a single test file
if gff_features is not None:
    test_file = os.path.join(CONFIG['data_dir'], 'EDO_1.dexeq_counts')
    processed_file = process_count_file(test_file, CONFIG['output_dir'], gff_features)
    print(f"Test file processing result: {processed_file}")

#%% 3. Process all files
if gff_features is not None:
    count_files = [f for f in os.listdir(CONFIG['data_dir']) if f.endswith('.dexeq_counts')]
    processed_files = []
    
    print(f"\nFound {len(count_files)} count files:")
    for f in count_files:
        print(f"- {f}")
        
    for file in count_files:
        file_path = os.path.join(CONFIG['data_dir'], file)
        if processed_file := process_count_file(file_path, CONFIG['output_dir'], gff_features):
            processed_files.append(processed_file)
            print(f"Successfully processed: {file}")
        else:
            print(f"Failed to process: {file}")

#%% 4. Create and process sample information
if gff_features is not None and processed_files:
    sample_info = pd.DataFrame([
        {
            'sample': os.path.basename(f).replace('.dexeq_counts', ''),
            'condition': 'EDO' if 'EDO' in f else 'ND1'
        }
        for f in processed_files
        if any(cond in f for cond in ['EDO', 'ND1'])
    ]).sort_values('sample').reset_index(drop=True)
    
    print("\nSample information:")
    print(sample_info)

#%% 5. Create DEXSeq dataset
if gff_features is not None and processed_files and not sample_info.empty:
    try:
        dxd = create_dexseq_dataset(sample_info, processed_files, CONFIG['dexseq_gff'])
        print("Successfully created DEXSeq dataset!")
    except Exception as e:
        print(f"Failed to create DEXSeq dataset: {e}")
