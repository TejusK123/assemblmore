#!/bin/bash
set -euo pipefail

READS=${1:-}
TR=${2:-}
MAX_ROUNDS=${3:-20}

if [ -z "$READS" ] || [ -z "$TR" ]; then
    echo "Usage: $0 <reads.fastq> <transposon_repeat.fasta> [max_rounds]"
    exit 1
fi

OUTPUT_PREFIX="${TR%.fasta}_polished.fasta"

seq_total_length() {
    local fasta="$1"
    awk '/^>/ {if(seqlen){s+=seqlen; seqlen=0}} !/^>/ {seqlen+=length($0)} END{if(seqlen) s+=seqlen; print s+0}' "$fasta"
}

prev_len=-1
round=1
current_input="$TR"

while [ "$round" -le "$MAX_ROUNDS" ]; do
    outdir="polish_round_${round}"
    mkdir -p "$outdir"
    echo "[polishTR] Round $round: running medaka_consensus (draft=$current_input)"

    medaka_consensus -i "$READS" -d "$current_input" -o "$outdir" -q || { echo "medaka_consensus failed on round $round" >&2; exit 1; }

    # find consensus fasta in output dir
    if [ -f "$outdir/consensus.fasta" ]; then
        result="$outdir/consensus.fasta"
    else
        result=$(find "$outdir" -maxdepth 1 -type f -name "*.fasta" -print -quit || true)
    fi

    if [ -z "${result:-}" ] || [ ! -f "$result" ]; then
        echo "No consensus fasta found in $outdir" >&2
        exit 1
    fi

    curr_len=$(seq_total_length "$result")
    echo "[polishTR] Round $round total sequence length: $curr_len"

    if [ "$prev_len" -ge 0 ]; then
        if [ "$curr_len" -eq "$prev_len" ]; then
            echo "[polishTR] Length unchanged (difference 0). Stopping polishing."
            cp "$result" "$OUTPUT_PREFIX"
            echo "[polishTR] Final polished fasta: $OUTPUT_PREFIX"
            exit 0
        else
            diff=$(( curr_len > prev_len ? curr_len - prev_len : prev_len - curr_len ))
            echo "[polishTR] Length change from previous round: $diff"
        fi
    fi

    prev_len=$curr_len
    current_input="$result"
    round=$((round + 1))
done

echo "[polishTR] Reached max rounds ($MAX_ROUNDS). Saving last consensus."
cp "$result" "$OUTPUT_PREFIX"
echo "[polishTR] Final polished fasta: $OUTPUT_PREFIX"

