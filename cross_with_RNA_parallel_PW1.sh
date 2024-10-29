#!/bin/bash
#SBATCH --job-name=cross_with_RNA_parallel
#SBATCH --account=kubacki.michal
#SBATCH --mem=312GB
#SBATCH --time=INFINITE
#SBATCH --exclusive
#SBATCH --ntasks=128
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/cross_with_RNA_parallel.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/cross_with_RNA_parallel.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

python /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/cross_with_RNA_parallel.py
