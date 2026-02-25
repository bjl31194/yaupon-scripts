#!/bin/bash
#SBATCH --job-name=snapper
#SBATCH --partition=batch
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=3-00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/scratch/bjl31194/logs/%x_%j.out
#SBATCH --error=/scratch/bjl31194/logs/%x_%j.error
  
# Set up job environment
set -o errexit  # Exit the script on any error
set -o nounset  # Treat any unset variables as an error
module --quiet purge  # Reset the modules to the system default
  
# Load the beast2 module
module load Beast/2.7.7-GCC-12.3.0-CUDA-12.1.1

# set parameters
DATADIR="/scratch/bjl31194/yaupon/wgs/plates1-5/vcfnew/phylo"

cd $DATADIR

# Run ruby script to create input XML file
#ruby /home/bjl31194/yaupon/yaupon-scripts/snapp_prep.rb -a SNAPPER -v Ilex_redrep_filter.vcf -s Ilex_redrep_NJtree.newick -t individuals.txt -c constraints.txt -m 2000 -l 100000

# Run snapper
beast -working -threads 16 snapper.xml