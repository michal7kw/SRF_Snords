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

# First, create a modified GFF3 file with strand-specific gene IDs
echo "Modifying gene IDs to be strand-specific..." >> /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/prepare_annotation.log

awk -F'\t' 'BEGIN{OFS="\t"} {
    if ($3 == "gene" || $3 == "exon") {
        split($9, attrs, ";");
        for (i in attrs) {
            if (attrs[i] ~ /^gene_id/) {
                gsub(/gene_id=/, "", attrs[i]);
                attrs[i] = "gene_id=" attrs[i] ($7 == "+" ? "plus" : "minus");
            }
        }
        $9 = "";
        for (i in attrs) {
            if ($9 != "") $9 = $9 ";";
            $9 = $9 attrs[i];
        }
    }
    print
}' /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.gff3 > /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.modified.gff3

echo "Running dexseq_prepare_annotation.py..." >> /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/prepare_annotation.log

python /home/kubacki.michal/.conda/envs/jupyter_nb/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.modified.gff3 \
    /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/gencode.v31.basic.annotation.dexseq.gff

echo "Job completed at $(date)" >> /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_Snords/DATA/logs/prepare_annotation.log
