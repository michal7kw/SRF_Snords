def check_fasta_version(fasta_path):
    with open(fasta_path, 'r') as f:
        # Read first few lines
        header = [next(f) for _ in range(5)]
    return '\n'.join(header)
print(check_fasta_version("DATA/hg38.fa"))

import os
import hashlib

def get_fasta_info(fasta_path):
    stats = os.stat(fasta_path)
    
    # Calculate MD5 of first 1MB (for quick checking)
    md5 = hashlib.md5()
    with open(fasta_path, 'rb') as f:
        md5.update(f.read(1024*1024))
    
    return {
        'size': stats.st_size,
        'modified_date': stats.st_mtime,
        'partial_md5': md5.hexdigest()
    }
print(get_fasta_info("DATA/hg38.fa"))

from Bio import SeqIO

def get_chromosome_lengths(fasta_path):
    lengths = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        lengths[record.id] = len(record.seq)
    return lengths
lengths = get_chromosome_lengths("DATA/hg38.fa")
print(lengths)

# Known hg38 chromosome 1 length: 248,956,422
if 'chr1' in lengths:
    print(f"Is chr1 length correct for hg38? {lengths['chr1'] == 248956422}")

import requests

def verify_ucsc_version(fasta_path):
    # UCSC hg38 chromosome 1 first 1000 bases
    ucsc_url = "https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr1;start=0;end=1000"
    
    # Get local sequence
    local_seq = ""
    for record in SeqIO.parse(fasta_path, "fasta"):
        if record.id == "chr1":
            local_seq = str(record.seq[:1000])
            break
    
    # Get UCSC sequence
    response = requests.get(ucsc_url)
    if response.ok:
        ucsc_seq = response.json()['dna'].upper()
        
        # Compare
        return {
            'matches_ucsc': local_seq == ucsc_seq,
            'local_start': local_seq[:50],
            'ucsc_start': ucsc_seq[:50]
        }
    return None
print(verify_ucsc_version("DATA/hg38.fa"))