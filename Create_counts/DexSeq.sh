#!/bin/bash
#SBATCH --job-name=DexSeq
#SBATCH --account=kubacki.michal
#SBATCH --mem=64GB
#SBATCH --time=INFINITE
#SBATCH --exclusive
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/logs/DexSeq.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/logs/DexSeq.out"

echo "my job start now" > DexSeq.log;
date >> DexSeq.log;

source /opt/common/tools/ric.cosr/miniconda3/bin/activate
conda activate jupyter_nb

# First, create the DEXSeq-formatted GFF
# python /home/kubacki.michal/.conda/envs/jupyter_nb/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
#     /beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/gencode.v31.basic.annotation.gff \
#     /beegfs/scratch/ric.broccoli/ric.broccoli/PW_RNA_seq_deep/gencode.v31.basic.annotation.DEXSeq.gff

Then use the DEXSeq-formatted GFF for counting
Alternative parallel approach using background processes
for i in EDO_1 EDO_2 EDO_3 ND1_1 ND1_2 ND1_3 PW1_1 PW1_2 PW1_3; do
    python /home/kubacki.michal/.conda/envs/jupyter_nb/lib/R/library/DEXSeq/python_scripts/dexseq_count.py \
    -p 'yes' \
    -s 'no' \
    -r pos \
    -f bam \
    -a 10 \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.dexseq.gff \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/$i.bam \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/Create_counts/output/$i.dexeq_counts &
done
wait

date >> DexSeq.log;
echo "all done!!" >> DexSeq.log