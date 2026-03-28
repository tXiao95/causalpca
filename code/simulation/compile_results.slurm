#!/bin/bash

#SBATCH --job-name=compile_CSDR
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=thsiao3@emory.edu
#SBATCH --time=00:05:00
#SBATCH --partition=short-cpu
#SBATCH --output=/scratch/thsiao3/output/compile_CSDR_%j.out
#SBATCH --error=/scratch/thsiao3/error/compile_CSDR_%j.err

# Load necessary modules or set any other environment variables if needed
module load singularity
# navigate to directory
cd /projects/hhchang/causal-pca/code/simulation

# Path to the Singularity image
SINGULARITY_IMAGE=~/singularity-images/geospatial_4.4.1.sif

# Path to the R script
R_SCRIPT=compile_results.R

# Launch Singularity Job
srun singularity exec --bind /projects:/projects --bind /scratch:/scratch $SINGULARITY_IMAGE Rscript $R_SCRIPT $1 $2