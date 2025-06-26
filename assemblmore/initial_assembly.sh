#!/bin/bash


MAPPED_CONTIGS=$1
READS=$2
ASSEMBLY=$3

rem_suffix="${MAPPED_CONTIGS%.sorted.paf}"
base_name_reference="${rem_suffix##*to_}"
base_name_reference="${base_name_reference%.fasta}"
base_name_reference="${base_name_reference%.fa}"
base_name_reference="${base_name_reference%.fna}"

python contig_placements.py "$MAPPED_CONTIGS" "$READS"

> ordered_and_oriented_to_${base_name_reference}_assembly.fasta

while IFS=$'\t' read -r id _ strand _ _; do
    if [ "$strand" == "-" ]; then
        seqkit grep -w 0 -ip "$id" "$ASSEMBLY" | seqkit seq -r -p -t dna >> ordered_and_oriented_to_${base_name_reference}_assembly.fasta
    else
        seqkit grep -w 0 -ip "$id" "$ASSEMBLY" >> ordered_and_oriented_to_${base_name_reference}_assembly.fasta
    fi
done < filtered_by_${base_name_reference}_contigs.tsv

