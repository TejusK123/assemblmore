#!/bin/bash
#SBATCH --account=rog
#SBATCH --partition=notchpeak
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=8
#SBATCH --mem=64G
#SBATCH --mail-user=u6067145@utah.edu
#SBATCH --mail-type=ALL
#SBATCH -J medaka1
#SBATCH -o mapping_run.log
#SBATCH -e mapping_run.err


#### Print Start Date ####
date



#### Variables ####


conda activate Nanopore


#### Work ####



./fill_gaps.sh ordered_and_oriented_assembly.fasta C_briggsae_AF16.fasta


#### Print End Date ####
date
