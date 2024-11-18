#!/bin/bash
#SBATCH --job-name=salmon_index
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/salmon_index.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/salmon_index.out"

source /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb/bin/activate
conda activate jupyter_nb   

set -e  # Exit on error

# Configuration
WORK_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/salmon_index"
GENCODE_VERSION="v38"
THREADS=16
GENOME="human"
LOG_FILE="${WORK_DIR}/salmon_index_build_$(date +%Y%m%d_%H%M%S).log"

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

# Function to check if program exists
check_command() {
    if ! command -v "$1" &> /dev/null; then
        log "ERROR: $1 could not be found. Please install it first."
        exit 1
    fi
}

# Function to check disk space (requires at least 50GB)
check_disk_space() {
    local required_space=50  # GB
    local available_space=$(df -BG "$WORK_DIR" | awk 'NR==2 {gsub("G",""); print $4}')
    
    if [ "$available_space" -lt "$required_space" ]; then
        log "ERROR: Insufficient disk space. Required: ${required_space}GB, Available: ${available_space}GB"
        exit 1
    fi
}

# Create work directory
mkdir -p "$WORK_DIR"

# Initialize log file
> "$LOG_FILE"

# Check required programs
log "Checking required programs..."
check_command salmon
check_command wget
check_command gunzip

# Log system information
log "System information:"
log "Salmon version: $(salmon --version)"
log "Working directory: $WORK_DIR"
log "Threads: $THREADS"

# Check disk space
log "Checking disk space..."
check_disk_space

# Set up paths
FASTA_GZ="${WORK_DIR}/gencode.${GENCODE_VERSION}.transcripts.fa.gz"
FASTA="${WORK_DIR}/gencode.${GENCODE_VERSION}.transcripts.fa"
INDEX_DIR="${WORK_DIR}/salmon_index_${GENCODE_VERSION}"

# Download transcriptome if needed
if [ ! -f "$FASTA_GZ" ]; then
    log "Downloading transcriptome..."
    BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode"
    SPECIES="Gencode_human"
    if [ "$GENOME" != "human" ]; then
        SPECIES="Gencode_mouse"
    fi
    URL="${BASE_URL}/${SPECIES}/release_${GENCODE_VERSION/v/}/${FASTA_GZ##*/}"
    
    wget -O "$FASTA_GZ" "$URL" || {
        log "ERROR: Failed to download transcriptome"
        exit 1
    }
    
    # Verify file size
    FILE_SIZE=$(stat -f%z "$FASTA_GZ" 2>/dev/null || stat -c%s "$FASTA_GZ")
    if [ "$FILE_SIZE" -lt 1000000 ]; then  # Less than 1MB
        log "ERROR: Downloaded file is too small ($FILE_SIZE bytes)"
        rm -f "$FASTA_GZ"
        exit 1
    fi
fi

# Decompress if needed
if [ ! -f "$FASTA" ]; then
    log "Decompressing FASTA file..."
    gunzip -c "$FASTA_GZ" > "$FASTA" || {
        log "ERROR: Failed to decompress FASTA"
        rm -f "$FASTA"
        exit 1
    }
fi

# Remove old index if it exists
if [ -d "$INDEX_DIR" ]; then
    log "Removing existing index directory..."
    rm -rf "$INDEX_DIR"
fi

# Build the index
log "Building Salmon index..."
log "Command: salmon index -t $FASTA -i $INDEX_DIR -p $THREADS --gencode -k 31 --keepDuplicates"

START_TIME=$(date +%s)

salmon index \
    -t "$FASTA" \
    -i "$INDEX_DIR" \
    -p "$THREADS" \
    --gencode \
    -k 31 \
    --keepDuplicates 2>&1 | tee -a "$LOG_FILE" || {
        log "ERROR: Failed to build Salmon index"
        rm -rf "$INDEX_DIR"
        exit 1
    }

END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Verify index
log "Verifying index..."
if [ ! -d "$INDEX_DIR" ]; then
    log "ERROR: Index directory not found"
    exit 1
fi

# Check for required files
REQUIRED_FILES=(
    "complete_ref_lens.bin"
    "ctable.bin"
    "ctg_offsets.bin"
    "duplicate_clusters.tsv"
    "info.json"
    "mphf.bin"
    "pos.bin"
    "pre_indexing.log"
    "rank.bin"
    "refseq.bin"
    "seq.bin"
    "versionInfo.json"
)

MISSING_FILES=0
for file in "${REQUIRED_FILES[@]}"; do
    if [ ! -f "${INDEX_DIR}/${file}" ]; then
        log "ERROR: Missing required file: $file"
        MISSING_FILES=1
    fi
done

if [ $MISSING_FILES -eq 1 ]; then
    log "ERROR: Index validation failed"
    exit 1
fi

# Calculate index size
INDEX_SIZE=$(du -sh "$INDEX_DIR" | cut -f1)

# Log completion
log "Index built successfully!"
log "Build time: ${DURATION} seconds"
log "Index size: ${INDEX_SIZE}"
log "Index location: $INDEX_DIR"
log "Log file: $LOG_FILE"

