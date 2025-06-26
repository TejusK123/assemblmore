#!/bin/bash
#SBATCH --account=rog
#SBATCH --partition=notchpeak
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=8
#SBATCH --mem=32G
#SBATCH --mail-user=u6067145@utah.edu
#SBATCH --mail-type=ALL
#SBATCH -J medaka1
#SBATCH -o medaka_run.log
#SBATCH -e medaka_run.err


#### Print Start Date ####
date



#### Variables ####


conda activate Nanopore


#### Work ####

medaka_consensus -i /uufs/chpc.utah.edu/common/home/rog-group2/Tejus/C_briggsae/C_briggsae_AF16.fastq -d /uufs/chpc.utah.edu/common/home/rog-group2/Tejus/Research/FLYE_ASSEMBLY/split_contigs/assembly_improvement_4/AF16_draft_penultimate.fasta -m r1041_e82_400bps_sup_v4.3.0 -q







#### Print End Date ####
date
