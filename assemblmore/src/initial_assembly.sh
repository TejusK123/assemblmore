#!/bin/bash


MAPPED_CONTIGS=$1
READS=$2
ASSEMBLY=$3

if [ -z "$MAPPED_CONTIGS" ] || [ -z "$READS" ] || [ -z "$ASSEMBLY" ]; then
    echo "Usage: $0 <mapped_contigs.paf> <reads.fastq> <assembly.fasta>"
    exit 1
fi


rem_suffix="${MAPPED_CONTIGS%.sorted.paf}"
base_name_reference="${rem_suffix##*to_}"
base_name_reference="${base_name_reference%.fasta}"
base_name_reference="${base_name_reference%.fa}"
base_name_reference="${base_name_reference%.fna}"


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python "$SCRIPT_DIR/contig_placements.py" "$MAPPED_CONTIGS" "$READS"

> ordered_and_oriented_to_${base_name_reference}_assembly.fasta

while IFS=$'\t' read -r id _ strand _ _ _ _; do
    if [ "$strand" == "-" ]; then
        seqkit grep -w 0 -ip "$id" "$ASSEMBLY" | seqkit seq -v -r -p -t dna >> ordered_and_oriented_to_${base_name_reference}_assembly.fasta
    else
        seqkit grep -w 0 -ip "$id" "$ASSEMBLY" >> ordered_and_oriented_to_${base_name_reference}_assembly.fasta
    fi
done < filtered_by_${base_name_reference}_contigs.tsv

# fill_gaps.sh ordered_and_oriented_to_${base_name_reference}_assembly.fasta "$READS"