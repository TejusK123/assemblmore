#!/bin/bash


# test=$1


# t1="${test%.sorted.paf}"
# t2="${t1##*to_}"

# echo "${t2}"

ALIGNMENTS=$1
ORDERINGS=$2
CONTIGS=$3
READS=$4

if [ -z "$ALIGNMENTS" ] || [ -z "$ORDERINGS" ] || [ -z "$CONTIGS" ] || [ -z "$READS" ]; then
    echo "Usage: $0 <alignments.paf> <orderings.tsv> <contigs.fasta> <reads.fasta>"
    exit 1
fi


SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python "$SCRIPT_DIR/span_contigs.py" "$ALIGNMENTS" "$ORDERINGS" "$CONTIGS" "$READS"