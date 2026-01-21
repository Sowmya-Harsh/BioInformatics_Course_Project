#!/bin/bash
#SBATCH --job-name=PINK1_Enhanced
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --output=pink1_enhanced_%j.out
#SBATCH --error=pink1_enhanced_%j.err
#SBATCH --constraint=skylake

echo "========================================"
echo "scRNAseq Analysis"
echo "Job ID: $SLURM_JOB_ID"
echo "Start time: $(date)"
echo "========================================"

# Load modules
source /etc/profile.d/modules.sh
module --force purge
module load env/release/2023b
module load lang/R-bundle-CRAN/2024.06-foss-2023b

# Set environment
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo ""
echo "Using $SLURM_CPUS_PER_TASK CPU cores"
echo "Allocated memory: ${SLURM_MEM_PER_NODE}MB"

# Create output directories
mkdir -p Results
mkdir -p result_images

# Run analysis
echo ""
echo "Running enhanced analysis..."
Rscript scRNAseq_Analysis.R

echo ""
echo "========================================"
echo "Analysis Complete!"
echo "End time: $(date)"
echo "========================================"
echo "Output images: result_images/ directory"
echo "Excel files: Results/ directory"
echo ""
echo "Generated files: plots and xlsx"