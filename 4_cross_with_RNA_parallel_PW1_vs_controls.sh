#!/bin/bash
#SBATCH --job-name=cross_with_RNA_parallel_PW1_vs_controls
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --ntasks-per-node=24
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/cross_with_RNA_parallel_PW1_vs_controls.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/cross_with_RNA_parallel_PW1_vs_controls.out"

# Set OpenBLAS and OpenMP thread limits
export OPENBLAS_NUM_THREADS=16
export OMP_NUM_THREADS=16

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb   

python /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/4_cross_with_RNA_parallel_PW1_vs_combined_controls.py

