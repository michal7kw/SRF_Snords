#%% Imports
import os
from typing import Set, List, Optional, Dict, Tuple
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.conversion import localconverter

# Enable automatic conversion between pandas and R dataframes
pandas2ri.activate()

# Import R packages
dexseq = importr('DEXSeq')

#%% Feature extraction function
def get_gff_features(gff_file: str) -> Tuple[Set[str], Dict[str, str]]:
    """Extract features from GFF file and create mapping between versioned and unversioned IDs."""
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file not found: {gff_file}")
        
    gff_features = set()
    id_mapping = {}
    
    with open(gff_file, 'r') as f:
        for line in f:
            if 'exonic_part' in line:
                fields = line.strip().split('\t')
                attrs = dict(item.strip().split(' ', 1) for item in fields[8].strip().split(';'))
                gene_id = attrs['gene_id'].strip('"')
                exon_num = attrs['exonic_part_number'].strip('"')
                
                # Store the exact GFF format
                feature_id = f"{gene_id}:E{exon_num}"
                gff_features.add(feature_id)
                
                # Create mapping from base ID to GFF format
                base_id = gene_id.split('.')[0]
                base_feature_id = f"{base_id}:E{exon_num}"
                id_mapping[base_feature_id] = feature_id
                
                # Store first few examples for debugging
                if len(gff_features) <= 5:
                    print(f"GFF Feature example: {feature_id}")
                    print(f"Base Feature example: {base_feature_id}")
    
    print(f"Loaded {len(gff_features)} features from GFF file")
    return gff_features, id_mapping

#%% Process count file function
def process_count_file(file_path: str, output_dir: str, gff_features: Set[str], 
                      id_mapping: Dict[str, str]) -> Optional[str]:
    """Process count file with exact DEXSeq feature ID format."""
    os.makedirs(output_dir, exist_ok=True)
    basename = os.path.basename(file_path)
    output_path = os.path.join(output_dir, f"processed_{basename}")
    
    try:
        # First, create a mapping of base IDs to counts
        count_dict = {}
        with open(file_path, 'r') as f:
            for line in f:
                if not line.startswith('_'):
                    feature_id, count = line.strip().split('\t')
                    if '":"' in feature_id:
                        gene_id, exon_num = feature_id.strip('"').split('":"')
                        # Create base feature ID without version
                        base_id = gene_id.split('.')[0]
                        base_feature_id = f"{base_id}:E{exon_num.zfill(3)}"
                        
                        # Store count for this base feature
                        if base_feature_id in id_mapping:
                            gff_feature_id = id_mapping[base_feature_id]
                            count_dict[gff_feature_id] = int(count)
        
        if not count_dict:
            print(f"No matching features found in {basename}")
            return None
            
        print(f"Found {len(count_dict)} matching features in {basename}")
        
        # Write output file using exact GFF feature IDs and their counts
        with open(output_path, 'w') as f:
            for gff_feature in sorted(gff_features):
                count = count_dict.get(gff_feature, 0)
                f.write(f"{gff_feature}\t{count}\n")
        
        # Debug: Show first few lines of output file
        print(f"\nFirst few lines of processed file {output_path}:")
        with open(output_path, 'r') as f:
            for i, line in enumerate(f):
                if i < 5:
                    print(f"Output line: {line.strip()}")
                else:
                    break
        
        return output_path
        
    except Exception as e:
        print(f"Error processing {basename}: {str(e)}")
        return None

#%% Create DEXSeq dataset function
def create_dexseq_dataset(sample_info: pd.DataFrame, processed_files: List[str], 
                         gff_file: str) -> ro.vectors.Vector:
    """Create DEXSeq dataset with proper R conversion."""
    try:
        # Ensure sample order matches count files
        sample_data = pd.DataFrame({
            'sample': sample_info['sample'],
            'condition': sample_info['condition']
        }).sort_values('sample')
        
        print("\nSample data:")
        print(sample_data)
        
        # Convert pandas DataFrame to R DataFrame
        with localconverter(ro.default_converter + pandas2ri.converter):
            sample_data_r = ro.conversion.py2rpy(sample_data)
        
        print("\nCreating DEXSeqDataSet...")
        dxd = dexseq.DEXSeqDataSetFromHTSeq(
            countfiles=StrVector(processed_files),
            sampleData=sample_data_r,
            design=Formula('~ sample + exon + condition:exon'),
            flattenedfile=gff_file
        )
        
        print("Estimating size factors...")
        dxd = dexseq.estimateSizeFactors(dxd)
        
        print("Estimating dispersions...")
        dxd = dexseq.estimateDispersions(dxd)
        
        return dxd
        
    except Exception as e:
        print(f"Error creating DEXSeq dataset: {str(e)}")
        raise

#%% Run analysis function
def run_analysis(gff_file: str, data_dir: str, output_dir: str) -> Optional[ro.vectors.Vector]:
    """Run DEXSeq analysis with version-aware feature handling."""
    try:
        print("Loading GFF features...")
        gff_features, id_mapping = get_gff_features(gff_file)
        
        print("\nProcessing count files...")
        count_files = [f for f in os.listdir(data_dir) if f.endswith('.dexeq_counts')]
        processed_files = []
        
        for file in count_files:
            file_path = os.path.join(data_dir, file)
            if processed_file := process_count_file(file_path, output_dir, gff_features, id_mapping):
                processed_files.append(processed_file)
                print(f"Successfully processed: {file}")
        
        if not processed_files:
            raise RuntimeError("No files were processed successfully")
        
        # Create sample information for EDO vs ND1 comparison
        sample_info = pd.DataFrame([
            {
                'sample': os.path.basename(f).replace('processed_', '').replace('.dexeq_counts', ''),
                'condition': 'EDO' if 'EDO' in f else 'ND1'
            }
            for f in processed_files
            if any(cond in f for cond in ['EDO', 'ND1'])
        ]).sort_values('sample').reset_index(drop=True)
        
        print("\nSample information:")
        print(sample_info)
        
        # Create DEXSeq dataset
        return create_dexseq_dataset(sample_info, processed_files, gff_file)
        
    except Exception as e:
        print(f"Analysis failed: {str(e)}")
        return None

#%% Run the analysis
if __name__ == "__main__":
    # Configuration
    CONFIG = {
        'gff_file': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.dexseq.gff",
        'data_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output",
        'output_dir': "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/output"
    }
    
    # Run analysis
    dxd = run_analysis(**CONFIG)



