#%% Imports
import os
from collections import Counter
import re

#%% Function to safely split gene ID
def split_gene_id(gene_id: str) -> tuple:
    """Safely split gene ID into base and version."""
    try:
        # Split on last dot to handle cases with multiple dots
        parts = gene_id.rsplit('.', 1)
        if len(parts) == 2:
            return parts[0], parts[1]
    except Exception as e:
        print(f"Error splitting gene ID {gene_id}: {e}")
    return None, None

#%% Function to extract gene versions from count files
def get_count_file_versions(count_file: str) -> tuple:
    """Extract and analyze gene versions from count file."""
    versions = []
    gene_examples = {}
    all_genes = set()
    
    print(f"\nAnalyzing count file: {os.path.basename(count_file)}")
    try:
        with open(count_file, 'r') as f:
            for line in f:
                if not line.startswith('_'):
                    feature_id = line.strip().split('\t')[0]
                    if '":"' in feature_id:
                        gene_id = feature_id.strip('"').split('":"')[0]
                        base_id, version = split_gene_id(gene_id)
                        if base_id:
                            all_genes.add(base_id)
                            if version:
                                versions.append(version)
                                if len(gene_examples) < 5 and base_id not in gene_examples:
                                    gene_examples[base_id] = gene_id
    
        version_counts = Counter(versions)
        print("\nVersion distribution:")
        for version, count in version_counts.most_common():
            print(f"Version {version}: {count} genes")
        
        print("\nSample gene IDs:")
        for base_id, full_id in list(gene_examples.items())[:5]:
            print(f"{base_id} -> {full_id}")
            
        return version_counts, gene_examples, all_genes
    
    except Exception as e:
        print(f"Error analyzing file: {e}")
        return Counter(), {}, set()

#%% Function to parse DEXSeq GFF attributes
def parse_gff_attributes(attr_string: str) -> dict:
    """Parse GFF attribute string into a dictionary."""
    attrs = {}
    for attr in attr_string.strip().split(';'):
        if attr:
            key_value = attr.strip().split(' ', 1)
            if len(key_value) == 2:
                key, value = key_value
                attrs[key] = value.strip('"')
    return attrs

#%% Function to analyze DEXSeq GFF file
def analyze_gff_file(gff_file: str) -> tuple:
    """Analyze gene versions in DEXSeq GFF file."""
    versions = []
    gene_examples = {}
    problematic_ids = set()
    all_genes = set()
    
    print(f"\nAnalyzing GFF file: {os.path.basename(gff_file)}")
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if len(fields) >= 9:
                        attrs = parse_gff_attributes(fields[8])
                        if 'gene_id' in attrs:
                            gene_id = attrs['gene_id']
                            base_id, version = split_gene_id(gene_id)
                            if base_id:
                                all_genes.add(base_id)
                                if version and len(gene_examples) < 5:
                                    if base_id not in gene_examples:
                                        versions.append(version)
                                        gene_examples[base_id] = gene_id
        
        version_counts = Counter(versions)
        print("\nVersion distribution in sample:")
        for version, count in version_counts.most_common():
            print(f"Version {version}: {count} genes")
        
        print(f"\nTotal unique genes found: {len(all_genes)}")
        print("\nSample gene IDs:")
        for base_id, full_id in list(gene_examples.items())[:5]:
            print(f"{base_id} -> {full_id}")
            
        return version_counts, gene_examples, all_genes
    
    except Exception as e:
        print(f"Error analyzing GFF file: {e}")
        raise

#%% Function to check gene ID overlap
def check_id_overlap(gff_genes: set, count_genes: set):
    """Check if gene IDs match between GFF and count files."""
    print("\nChecking gene ID overlap:")
    
    common_genes = gff_genes & count_genes
    only_gff = gff_genes - count_genes
    only_counts = count_genes - gff_genes
    
    print(f"\nTotal genes in GFF: {len(gff_genes)}")
    print(f"Total genes in counts: {len(count_genes)}")
    print(f"Genes in both files: {len(common_genes)}")
    
    print("\nSample overlapping genes:")
    for gene in list(common_genes)[:5]:
        print(f"- {gene}")
    
    if len(only_counts) > 0:
        print("\nWarning: Some genes in count file not found in GFF!")
        print("Sample missing genes:")
        for gene in list(only_counts)[:5]:
            print(f"- {gene}")

#%% Run the analysis
if __name__ == "__main__":
    # Paths
    count_file = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output/EDO_1.dexeq_counts"
    current_gff = "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.dexseq.gff"
    
    # Analyze files
    count_versions, count_examples, count_genes = get_count_file_versions(count_file)
    gff_versions, gff_examples, gff_genes = analyze_gff_file(current_gff)
    
    # Check overlap
    check_id_overlap(gff_genes, count_genes)
    
    print("\nConclusion:")
    print("Using Gencode v31 throughout the pipeline.")
    print("Version numbers in file formats may differ but this doesn't affect the analysis.")


##############################################################################3

# %%
with open("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output/EDO_1.dexeq_counts") as f:
    # Print first 5 non-comment lines
    lines = [next(f) for _ in range(5)]
    print("First few lines of EDO_1.dexeq_counts:")
    for line in lines:
        print(line.strip())
# %%
with open("/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.dexseq.gff"   ) as f:
    # Print first 5 non-comment lines
    lines = [next(f) for _ in range(5)]
    print("First few lines of gencode.v31.basic.annotation.dexseq.gff:")
    for line in lines:
        print(line.strip())
# %%
