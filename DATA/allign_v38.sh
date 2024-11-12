#!/bin/bash
#SBATCH --job-name=alignment
#SBATCH --account=kubacki.michal
#SBATCH --mem=128GB
#SBATCH --time=INFINITE
#SBATCH --exclusive
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/alignment.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/alignment.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Set paths
FASTQ_DIR="/beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/fastq/PW_RNA_seq_deep"
OUTPUT_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/v38"
STAR_INDEX="$OUTPUT_DIR/STAR_index"
GENOME_FASTA="$OUTPUT_DIR/hg38.fa"
GTF_FILE="$OUTPUT_DIR/gencode.v38.annotation.gtf"

# First, remove the existing STAR index since it seems to be corrupted
rm -rf "$STAR_INDEX"

# Create output directory
mkdir -p $OUTPUT_DIR

# Create STAR index
create_star_index() {
    echo "Creating STAR index..."
    mkdir -p $STAR_INDEX
    
    # Add error checking and verbose output
    echo "Using genome: $GENOME_FASTA"
    echo "Using GTF: $GTF_FILE"
    
    if [ ! -f "$GENOME_FASTA" ]; then
        echo "ERROR: Genome file not found: $GENOME_FASTA"
        exit 1
    fi
    
    if [ ! -f "$GTF_FILE" ]; then
        echo "ERROR: GTF file not found: $GTF_FILE"
        exit 1
    fi
    
    STAR --runMode genomeGenerate \
         --genomeDir $STAR_INDEX \
         --genomeFastaFiles $GENOME_FASTA \
         --sjdbGTFfile $GTF_FILE \
         --sjdbOverhang 100 \
         --runThreadN 8 \
         --outFileNamePrefix "$STAR_INDEX/"  # Add output prefix
    
    echo "STAR index creation complete."
}

# Check if STAR index exists, if not create it
if [ ! -d "$STAR_INDEX" ]; then
    create_star_index
else
    echo "STAR index already exists. Skipping index creation."
fi

# Function to run STAR alignment
run_star() {
    sample=$1
    
    # Check if input files exist
    if [ ! -f "$FASTQ_DIR/${sample}/${sample}_L001_R1_001.fastq.gz" ] || [ ! -f "$FASTQ_DIR/${sample}/${sample}_L001_R2_001.fastq.gz" ]; then
        echo "Error: Input files not found for sample $sample"
        return 1
    fi
    
    mkdir -p "$OUTPUT_DIR/$sample"

    STAR \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$FASTQ_DIR/${sample}/${sample}_L001_R1_001.fastq.gz" "$FASTQ_DIR/${sample}/${sample}_L001_R2_001.fastq.gz" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$OUTPUT_DIR/$sample/" \
        --outSAMtype BAM SortedByCoordinate \
        --runThreadN 32 \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 3 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000
}

# Replace the sample reading section with direct sample list since samplesheet.csv isn't found
samples=("EDO_1" "EDO_2" "EDO_3" "ND1_1" "ND1_2" "ND1_3" "PW1_1" "PW1_2" "PW1_3")
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    run_star "$sample"
done

# Modify BAM file handling to check for existence
for sample in "${samples[@]}"; do
    if [ -f "$OUTPUT_DIR/$sample/Aligned.sortedByCoord.out.bam" ]; then
        mv "$OUTPUT_DIR/$sample/Aligned.sortedByCoord.out.bam" "$OUTPUT_DIR/${sample}.bam"
        samtools index "$OUTPUT_DIR/${sample}.bam"
    else
        echo "Warning: BAM file not found for sample $sample"
    fi
done

# Run MultiQC
multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR/multiqc_report"

echo "Alignment complete. Results are in $OUTPUT_DIR"
