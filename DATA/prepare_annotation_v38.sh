#!/bin/bash
#SBATCH --job-name=prepare_annotation
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --exclusive
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/prepare_annotation.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/prepare_annotation.out"

echo "Job started at $(date)" > /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/prepare_annotation.log

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

echo "Running dexseq_prepare_annotation.py..." >> /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/prepare_annotation.log

python /home/kubacki.michal/.conda/envs/jupyter_nb/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/v38/gencode.v38.annotation.gtf \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/v38/gencode.v38.annotation.dexseq.gff

echo "Job completed at $(date)" >> /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/prepare_annotation.log
