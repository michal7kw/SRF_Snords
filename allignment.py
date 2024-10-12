# sequence_utils.py
import logging
from typing import Optional
from Bio import Seq, SeqIO, Entrez

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def generate_antisense(sequence: str) -> str:
    """Generate the antisense sequence."""
    return str(Seq.Seq(sequence).reverse_complement())

def fetch_gene_sequence(gene_name: str, email: str) -> Optional[str]:
    """Fetch the gene sequence from NCBI."""
    Entrez.email = email
    try:
        handle = Entrez.esearch(db="nucleotide", term=f"{gene_name}[Gene Name] AND human[Organism]")
        record = Entrez.read(handle)
        if record["IdList"]:
            gene_id = record["IdList"][0]
            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
            gene_record = SeqIO.read(handle, "genbank")
            return str(gene_record.seq)
    except Exception as e:
        logger.error(f"Error fetching sequence for {gene_name}: {str(e)}")
    return None

# alignment_utils.py
from typing import List
import parasail

def align_sequence(query: str, target: str, allowed_mismatches: int = 1) -> int:
    """Perform local alignment and return the best score."""
    result = parasail.sw_trace_striped_16(query, target, 10, 1, parasail.dnafull)
    if result.score >= len(query) - allowed_mismatches:
        return result.score
    return 0

def check_alignment(antisense_seq: str, gene_list: List[str], email: str, 
                    allowed_mismatches: int, gene_sequences: Optional[dict] = None) -> List[str]:
    """Check alignment of antisense sequence with the given list of genes."""
    results = []
    for gene in gene_list:
        if gene_sequences is not None:
            gene_seq = gene_sequences.get(gene)
        else:
            gene_seq = fetch_gene_sequence(gene, email)
        if gene_seq:
            score = align_sequence(antisense_seq, gene_seq, allowed_mismatches)
            if score == len(antisense_seq):  # Perfect match
                results.append(gene)
            logger.info(f"Score for {gene}: {score}")
        else:
            logger.warning(f"Gene sequence not found for {gene}")
    return results

# main.py
import os
import json
import configparser
from sequence_utils import generate_antisense, fetch_gene_sequence
from alignment_utils import check_alignment

def main():
    # Load configuration
    config = configparser.ConfigParser()
    config.read('config.ini')

    email = config['NCBI']['email']
    working_dir = config['Paths']['working_directory']

    # Set working directory
    os.chdir(working_dir)
    logger.info(f"Current working directory: {os.getcwd()}")

    # Define sequences and genes
    ASE1 = 'AACATTCCTTGGAAAAG'
    ASE2 = 'CGTCATTCTCATCGGAA'
    selected = ["GEMIN4", "YME1L1", "VCP", "DARS1", "POLR2A", "ORAI1", "POLRMT", "CBX2", "NSD2", "TLE1", "CACNG1", "ZMYND11", "DCAF5", "MED25", "SEMA6A-AS2", "POLD3", "RFT1", "SETDB1", "POLE", "CHD4", "POLR2C", "HTT", "MED14", "PER3", "TOR3A", "EEF1A2", "IL33", "SETD5", "GNRH1", "POLR2E", "IGF1R", "IGF2R", "MTMR1", "EEF2", "RBPMS", "EIF2S3B", "NIPBL", "VCP"]

    # Generate antisense sequences
    antisense_ASE1 = generate_antisense(ASE1)
    antisense_ASE2 = generate_antisense(ASE2)
    logger.info(f"Antisense sequence ASE1: {antisense_ASE1}")
    logger.info(f"Antisense sequence ASE2: {antisense_ASE2}")

    # Fetch and store gene sequences
    gene_sequences = {}
    for gene in selected:
        sequence = fetch_gene_sequence(gene, email)
        if sequence:
            gene_sequences[gene] = sequence
        else:
            logger.warning(f"Failed to fetch sequence for {gene}")

    logger.info(f"Successfully fetched sequences for {len(gene_sequences)} out of {len(selected)} genes")

    # Save gene sequences
    with open('./output/gene_sequences.json', 'w') as f:
        json.dump(gene_sequences, f)
    logger.info("Gene sequences saved to gene_sequences.json")

    # Perform alignments
    aligned_genes_ASE1 = check_alignment(antisense_ASE1, selected, email, allowed_mismatches=3, gene_sequences=gene_sequences)
    aligned_genes_ASE2 = check_alignment(antisense_ASE2, selected, email, allowed_mismatches=1, gene_sequences=gene_sequences)

    # Print results
    logger.info("Genes with perfect alignment for ASE1:")
    for gene in aligned_genes_ASE1:
        logger.info(gene)

    logger.info("Genes with perfect alignment for ASE2:")
    for gene in aligned_genes_ASE2:
        logger.info(gene)

if __name__ == "__main__":
    main()