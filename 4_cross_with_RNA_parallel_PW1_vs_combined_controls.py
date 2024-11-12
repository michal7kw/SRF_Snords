# %%
import os
from typing import Set, List, Optional, Dict, Tuple
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.conversion import localconverter

# Set environment variables for threading
os.environ["OPENBLAS_NUM_THREADS"] = "16"
os.environ["OMP_NUM_THREADS"] = "16"

# Enable automatic conversion between pandas and R dataframes
pandas2ri.activate()

# Import R packages
dexseq = importr('DEXSeq')
deseq2 = importr('DESeq2')
stats = importr('stats')
base = importr('base')

def clean_gff_file(gff_file: str) -> str:
    """Clean GFF file and return path to cleaned file."""
    clean_file = gff_file + ".clean"
    if not os.path.exists(clean_file):
        with open(gff_file, 'r') as f, open(clean_file, 'w') as out:
            for line in f:
                if not line.startswith('#'):
                    out.write(line)
    return clean_file

def clean_count_file(count_file: str) -> str:
    """Clean count file and return path to cleaned file."""
    clean_file = count_file + ".clean"
    if not os.path.exists(clean_file):
        with open(count_file, 'r') as f, open(clean_file, 'w') as out:
            for line in f:
                if not line.startswith('_'):
                    out.write(line)
    return clean_file

def load_and_combine_count_data(data_dir: str, control_groups: List[str], treatment_group: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load count data and combine control replicates into a single control group.
    
    Args:
        data_dir: Directory containing count files
        control_groups: List of control group names to combine
        treatment_group: Name of treatment group
        
    Returns:
        Tuple of (combined counts DataFrame, sample table DataFrame)
    """
    count_files = {}
    sample_data = []
    
    # Load control group files and create sample data
    for ctrl_idx, group in enumerate(control_groups, start=1):
        for rep in range(1, 4):  # Assuming 3 replicates per group
            sample = f"{group}_{rep}"
            file_path = os.path.join(data_dir, f"{sample}.dexeq_counts")
            
            if os.path.exists(file_path):
                count_files[sample] = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
                sample_data.append({
                    'sample': sample,
                    'condition': 'control',
                    'individual': str(rep),
                    'original_group': group
                })
            else:
                raise FileNotFoundError(f"Count file not found: {file_path}")
    
    # Load treatment group files
    for rep in range(1, 4):  # Assuming 3 replicates
        sample = f"{treatment_group}_{rep}"
        file_path = os.path.join(data_dir, f"{sample}.dexeq_counts")
        
        if os.path.exists(file_path):
            count_files[sample] = pd.read_csv(file_path, sep='\t', header=0, index_col=0)
            sample_data.append({
                'sample': sample,
                'condition': 'treated',
                'individual': str(rep),
                'original_group': treatment_group
            })
        else:
            raise FileNotFoundError(f"Count file not found: {file_path}")
    
    # Combine all count files
    counts_df = pd.concat(count_files.values(), axis=1)
    counts_df.columns = count_files.keys()
    
    # Create sample table
    sample_table = pd.DataFrame(sample_data)
    sample_table.set_index('sample', inplace=True)
    
    # Ensure sample table matches count data columns
    sample_table = sample_table.loc[counts_df.columns]
    
    return counts_df, sample_table

def convert_dexseq_results_to_df(dexseq_results):
    """Convert DEXSeqResults object to pandas DataFrame."""
    try:
        # Convert DEXSeqResults to a data.frame in R
        r_df = ro.r('''
            function(res) {
                df <- as.data.frame(res)
                # Remove any complex columns that might cause conversion issues
                df$genomicData <- NULL
                df$countData <- NULL
                df$transcripts <- NULL
                return(df)
            }
        ''')(dexseq_results)
        
        # Convert R data.frame to pandas DataFrame
        with localconverter(ro.default_converter + pandas2ri.converter):
            results_df = ro.conversion.rpy2py(r_df)
            
        return results_df
    
    except Exception as e:
        print(f"Error converting DEXSeq results: {str(e)}")
        raise

def run_dexseq_analysis(counts_df: pd.DataFrame, 
                       sample_table: pd.DataFrame, 
                       gff_file: str,
                       data_dir: str,
                       min_count: int = 10,
                       n_cores: int = 32) -> ro.vectors.ListVector:
    """Run DEXSeq differential exon usage analysis with parallel processing."""
    
    # Clean GFF file and prepare count files
    clean_gff = clean_gff_file(gff_file)
    count_files = []
    for sample in counts_df.columns:
        original_file = os.path.join(data_dir, f"{sample}.dexeq_counts")
        if not os.path.exists(original_file):
            raise FileNotFoundError(f"Count file not found: {original_file}")
        clean_file = clean_count_file(original_file)
        count_files.append(clean_file)
    
    # Create a clean sample table with only necessary columns
    clean_sample_table = sample_table[['condition', 'individual']].copy()
    
    # Add nested individual identifier
    clean_sample_table['ind.n'] = pd.Categorical(
        [int(x) for x in clean_sample_table['individual']],
        categories=[1, 2, 3]
    )
    
    # Ensure all relevant columns are categorical
    clean_sample_table['condition'] = pd.Categorical(
        clean_sample_table['condition'],
        categories=['control', 'treated']
    )
    clean_sample_table['individual'] = pd.Categorical(
        clean_sample_table['individual'],
        categories=['1', '2', '3']
    )
    
    print("\nAnalysis Setup:")
    print(f"Number of samples: {len(count_files)}")
    print(f"Sample table shape: {clean_sample_table.shape}")
    print(f"Using {n_cores} CPU cores")
    print("\nSample table contents:")
    print(clean_sample_table)
    
    # Convert pandas objects to R
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_sample_table = ro.conversion.py2rpy(clean_sample_table)
        
        # Explicitly convert columns to R factors
        ro.r('''
            convert_to_factor <- function(df) {
                for(col in colnames(df)) {
                    df[[col]] <- as.factor(df[[col]])
                }
                return(df)
            }
        ''')
        r_sample_table = ro.r['convert_to_factor'](r_sample_table)    

    try:
        r_count_files = StrVector(count_files)
        
        print("\nCreating DEXSeqDataSet...")
        initial_design = Formula('~ individual + exon')
        
        # Set up parallel processing
        ro.r('''
            library(BiocParallel)
            setup_parallel <- function(cores) {
                return(MulticoreParam(workers = cores))
            }
        ''')
        
        bpparam = ro.r['setup_parallel'](n_cores)
        
        # Create DEXSeqDataSet
        dxd = dexseq.DEXSeqDataSetFromHTSeq(
            countfiles=r_count_files,
            sampleData=r_sample_table,
            design=initial_design,
            flattenedfile=clean_gff
        )
        
        # Filter low count exons
        print("\nFiltering low count exons...")
        dxd = ro.r('''
            function(dds, min_count) {
                exon_means <- rowMeans(counts(dds))
                keep <- exon_means >= min_count
                dds_filtered <- dds[keep,]
                message(sprintf("Kept %d out of %d exons after filtering", 
                              sum(keep), length(keep)))
                return(dds_filtered)
            }
        ''')(dxd, min_count)
        
        # Size factors estimation
        print("\nEstimating size factors...")
        dxd = ro.r('''
            function(dds) {
                if("sizeFactor" %in% colnames(colData(dds))) {
                    colData(dds)$sizeFactor <- NULL
                }
                
                tryCatch({
                    dds <- estimateSizeFactors(dds)
                }, error = function(e) {
                    message("Using alternative size factor estimation...")
                    geometric_means <- apply(counts(dds), 1, function(row) {
                        exp(mean(log(row[row > 0])))
                    })
                    size_factors <- apply(counts(dds), 2, function(col) {
                        median((col / geometric_means)[geometric_means > 0])
                    })
                    sizeFactors(dds) <- size_factors
                })
                return(dds)
            }
        ''')(dxd)
        
        # Create nested design formula
        nested_design = Formula('~ condition + condition:ind.n + condition:exon')
        
        # Update design
        dxd = ro.r('''
            function(dds, design) {
                cd <- colData(dds)
                for(col in colnames(cd)) {
                    if(col != "sizeFactor") {
                        cd[[col]] <- as.factor(cd[[col]])
                    }
                }
                colData(dds) <- cd
                design(dds) <- design
                return(dds)
            }
        ''')(dxd, nested_design)
        
        # Dispersion estimation
        print("\nEstimating dispersions...")
        dxd = ro.r('''
            function(dds, BPPARAM) {
                tryCatch({
                    dds <- estimateDispersions(dds, BPPARAM=BPPARAM)
                }, error = function(e) {
                    message("Using robust dispersion estimation...")
                    dds <- estimateDispersions(dds, fitType="local", BPPARAM=BPPARAM)
                })
                return(dds)
            }
        ''')(dxd, bpparam)
        
        # Test for differential exon usage
        print("\nTesting for differential exon usage...")
        reduced_design = Formula('~ condition + condition:ind.n')
        
        dxd = ro.r('''
            function(dds, reduced, full, BPPARAM) {
                options(warn=-1)
                dds <- tryCatch({
                    testForDEU(dds, reducedModel=reduced, fullModel=full, BPPARAM=BPPARAM)
                }, error = function(e) {
                    message("Error in DEU testing, trying alternative approach...")
                    testForDEU(dds, reducedModel=reduced, fullModel=full,
                              BPPARAM=BPPARAM)
                })
                options(warn=0)
                return(dds)
            }
        ''')(dxd, reduced_design, nested_design, bpparam)
        
        # Estimate fold changes
        print("\nEstimating exon fold changes...")
        dxd = ro.r('''
            function(dds, fitExpToVar, BPPARAM) {
                options(warn=-1)
                dds <- tryCatch({
                    estimateExonFoldChanges(dds, fitExpToVar=fitExpToVar, BPPARAM=BPPARAM)
                }, error = function(e) {
                    message("Using alternative fold change estimation...")
                    estimateExonFoldChanges(dds, fitExpToVar=fitExpToVar,
                                          BPPARAM=BPPARAM)
                })
                options(warn=0)
                return(dds)
            }
        ''')(dxd, "condition", bpparam)
        
        # Extract results
        print("\nExtracting results...")
        results = ro.r('''
            function(dds) {
                res <- DEXSeqResults(dds)
                res <- res[!is.na(res$padj), ]
                return(res)
            }
        ''')(dxd)
        
        return results
        
    except Exception as e:
        print(f"\nError during DEXSeq analysis: {str(e)}")
        raise

# %%
# Define paths
data_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output_v38"
working_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords"
# gff_file = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.dexseq.gff"
gff_file = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/v38/gencode.v38.annotation.dexseq.gff"
output_dir = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/output_v38"

# Define control groups and treatment group
control_groups = ['EDO', 'ND1']
treatment_group = 'PW1'

os.chdir(working_dir)

# Load and combine count data with combined controls
counts_df, sample_table = load_and_combine_count_data(data_dir, control_groups, treatment_group)

print("\nSample Table Summary:")
print(sample_table)
print("\nCount Data Shape:", counts_df.shape)

# Run DEXSeq analysis
dexseq_results = run_dexseq_analysis(counts_df, sample_table, gff_file, data_dir, n_cores=16)

# Convert results to pandas DataFrame
results_df = convert_dexseq_results_to_df(dexseq_results)

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Save results
output_file = os.path.join(output_dir, 'dexseq_results_PW1_vs_combined_controls.csv')
results_df.to_csv(output_file)

print("\nAnalysis Summary:")
print(f"Total features tested: {len(results_df)}")
if 'padj' in results_df.columns:
    print(f"Significant features (padj < 0.1): {(results_df['padj'] < 0.1).sum()}")
    print(f"Significant features (padj < 0.05): {(results_df['padj'] < 0.05).sum()}")
else:
    print("Note: No adjusted p-values found in results")

print("\nColumns in results:")
print(results_df.columns.tolist())
