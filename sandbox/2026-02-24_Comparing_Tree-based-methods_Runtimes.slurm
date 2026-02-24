#!/bin/bash
#SBATCH --job-name=tree_runtime
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=thsiao3@emory.edu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=10GB
#SBATCH --time=24:00:00
#SBATCH --partition=chang,day-long-cpu
#SBATCH --output=/scratch/%u/output/tree_runtime_%j.log
#SBATCH --error=/scratch/%u/error/tree_runtime_%j.log

# Load Singularity Module
module load singularity

# Execute R script using Singularity Image
singularity exec --bind /projects:/projects --bind /scratch:/scratch ~/singularity-images/geospatial_4.4.1.sif Rscript "2026-02-24_Comparing_Tree-based-methods_Runtimes.R"
