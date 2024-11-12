#!/bin/bash
#SBATCH --job-name=DexSeq
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --cpus-per-task=4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=0-8
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/logs/DexSeq_%A_%a.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/logs/DexSeq_%A_%a.out"

# Define array of sample names
samples=(EDO_1 EDO_2 EDO_3 ND1_1 ND1_2 ND1_3 PW1_1 PW1_2 PW1_3)
# Get current sample from array using SLURM_ARRAY_TASK_ID
sample=${samples[$SLURM_ARRAY_TASK_ID]}

echo "Processing sample: $sample"
echo "Job started" > DexSeq_${sample}.log
date >> DexSeq_${sample}.log

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# Process single sample with 4 cores
python /home/kubacki.michal/.conda/envs/jupyter_nb/lib/R/library/DEXSeq/python_scripts/dexseq_count.py \
    -p 'yes' \
    -s 'no' \
    -r pos \
    -f bam \
    -a 10 \
    -f bam \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/v31/gencode.v31.basic.annotation.dexseq.gff \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/v31/${sample}.bam \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output_v31/${sample}_2.dexeq_counts

# Clean up output file
sed 's/\"//g' /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output_v31/${sample}_2.dexeq_counts | \
    grep -v "^_" > /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output_v31/${sample}_2.fixed.dexeq_counts

date >> DexSeq_${sample}.log
echo "all done!!" >> DexSeq_${sample}.log

# # -p 'yes' \  # Paired-end reads
# # -s 'no'  \  # Strand-specific protocol (no = unstranded)
# # -r pos   \  # Order reads by position
# # -f bam   \  # Input file format
# # -a 10    \  # Minimum alignment quality score