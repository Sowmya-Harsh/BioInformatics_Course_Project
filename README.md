# Single Cell RNA Analysis

**Bioinformatics Course Project** - HPC-based single-cell RNA-seq analysis.

## Dataset
- Provided by Professor  
- Processed on High Performance Computing cluster

## Code
All analysis scripts in `codes/`:
- `install_packages.R` - Environment setup
- `scRNAseq_Analysis.R` - Main single-cell analysis pipeline  
- `RunAnalysis.slurm` - HPC batch submission script

## Usage
Run on HPC cluster:
```bash
sbatch codes/RunAnalysis.slurm
