#!/bin/bash
#SBATCH --job-name=rna_seq_analysis3
#SBATCH --account=kubacki.michal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=256GB
#SBATCH --time=48:00:00
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/3_%x_%j.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/3_%x_%j.out"

# Set strict error handling
set -e
set -o pipefail

# Change to working directory
cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords

# Activate conda environment
source /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb/bin/activate

# Set environment variables for optimal performance
export OPENBLAS_NUM_THREADS=64
export MKL_NUM_THREADS=64
export OMP_NUM_THREADS=64
export NUMEXPR_NUM_THREADS=64
export VECLIB_MAXIMUM_THREADS=64
export NUMBA_NUM_THREADS=64

# Set R-specific environment variables
export R_LIBS_USER="/beegfs/scratch/ric.broccoli/kubacki.michal/R/library"
export R_MAX_NUM_DLLS=150
export R_ENABLE_JIT=0

# Set temporary directory to fast local storage
export TMPDIR=/tmp
mkdir -p $TMPDIR

# Memory management settings
export PYTHONMALLOC=malloc
export MALLOC_TRIM_THRESHOLD_=100000
export PYTHONGC=1

# Add performance monitoring
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Memory available: $(free -h)"
echo "CPU info: $(lscpu | grep 'Model name' | cut -f 2 -d ":")"

# Create a temporary Python script
cat << 'EOF' > run_analysis.py
import os
import sys

# Set environment variables
os.environ['MKL_NUM_THREADS'] = '64'
os.environ['OPENBLAS_NUM_THREADS'] = '64'

# Add project directory to Python path
sys.path.append('/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords')

# Import and run main function - using string import to handle filename starting with number
import importlib
module = importlib.import_module('4_cross_with_RNA')
main = module.main
main()
EOF

# Add memory monitoring
monitor_memory() {
    while true; do
        free -h >> memory_usage.log
        sleep 300  # Log every 5 minutes
    done
}
monitor_memory &
MONITOR_PID=$!

# Run the analysis script with memory limits
time python -O run_analysis.py

# Kill memory monitor
kill $MONITOR_PID

# Cleanup
rm run_analysis.py

# Cleanup and reporting
echo "Job finished at: $(date)"
echo "Final memory status: $(free -h)"

# Archive logs
LOG_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/logs/archive"
mkdir -p $LOG_DIR
cp slurm-${SLURM_JOB_ID}.out ${LOG_DIR}/
cp slurm-${SLURM_JOB_ID}.err ${LOG_DIR}/
