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

#minimap2 -a /uufs/chpc.utah.edu/common/home/rog-group2/Tejus/Research/FLYE_ASSEMBLY/gapfilling/C_briggsae_AF16.fasta 
#minimap2 -a -Y /uufs/chpc.utah.edu/common/home/rog-group2/Tejus/Research/FLYE_ASSEMBLY/split_contigs/contig_57_last_100k.fasta /uufs/chpc.utah.edu/common/home/rog-group2/Tejus/Research/FLYE_ASSEMBLY/gapfilling/C_briggsae_AF16.fasta > reads_mapped_to_contig_57_last_100k.sam

#minimap2 -ax map-ont ../assembly.fasta ../../../C_briggsae/C_briggsae_AF16.fastq > reads_mapped_to_full_assembly.sam

#minimap2 -ax map-ont ../split_contigs/contig_45_revcomp_last_100k.fasta ../../../C_briggsae/C_briggsae_AF16.fastq > reads_mapped_to_contig_45_revcomp_last_100k.sam

#minimap2 -ax map-ont ../split_contigs/contig_45_revcomp_len_15651416.fasta ../../../C_briggsae/C_briggsae_AF16.fastq > reads_mapped_to_contig_45_whole.sam

#minimap2 -ax map-ont ../split_contigs/contig_36_first_100k.fasta ../../../C_briggsae/C_briggsae_AF16.fastq > reads_mapped_to_contig_36_first_100k.sam


#minimap2 -ax map-ont ../split_contigs/contig_36_last_100k.fasta ../../../C_briggsae/C_briggsae_AF16.fastq > reads_mapped_to_contig_36_last_100k.sam


#for i in ../split_contigs/contig*{first,last}*.fasta;
#do ./fill_gaps.sh $i ../../../C_briggsae/C_briggsae_AF16.fasta;
#done


#./fill_gaps.sh ../split_contigs/contig_57_revcomp_last_100k.fasta ../../../C_briggsae/C_briggsae_AF16.fastq
#./fill_gaps.sh ../split_contigs/contig_57_revcomp_first_100k.fasta ../../../C_briggsae/C_briggsae_AF16.fastq
#./fill_gaps.sh ../split_contigs/contig_57_revcomp_extension.fasta ../../../C_briggsae/C_briggsae_AF16.fastq

#./fill_gaps.sh assembly1.03.fasta ../../../../C_briggsae/C_briggsae_AF16.fastq
#./fill_gaps.sh assembly1.03_concise2.fasta ../../../../C_briggsae/C_briggsae_AF16.fastq

#./fill_gaps.sh AF16_draft_penultimate.fasta ../../../../C_briggsae/C_briggsae_AF16_50kb.fastq

#./fill_gaps.sh /uufs/chpc.utah.edu/common/home/rog-group2/Tejus/Research/References/c_briggsae/C_briggsae_45S_rDNA.fasta ../../../../C_briggsae/C_briggsae_AF16.fasta

./fill_gaps.sh ordered_and_oriented_assembly.fasta C_briggsae_AF16.fasta


#### Print End Date ####
date
