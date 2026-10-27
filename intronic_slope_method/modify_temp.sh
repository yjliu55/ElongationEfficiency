#!/bin/bash
set -euo pipefail

# === CONFIG ===
MIN_READ_SUPPORT=3
MAX_SHIFT=5
MIN_SIGNAL=5
INTRON_BED="TAIR10_pure_intron.bed"
EXPRESSED_GENES="genes.txt"
GENE_MATCH_SCRIPT="match_geneID.py"
SLOPE_SCRIPT="calculate_slope.py"

# === USAGE ===
if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <SJ.out.tab files...>"
  exit 1
fi

mkdir -p temp

# === STEP 1: Filter SJ.out.tab files by read support and collect base names ===
declare -a SJ_IDS=()

echo "Filtering splice junctions..."
for sj_file in "$@"; do
  base=$(basename "$sj_file" _SJ.out.tab)
  SJ_IDS+=("$base")
  awk -v min="$MIN_READ_SUPPORT" '($7 >= min)' "$sj_file" > "temp/selected_${base}_SJ.txt"
done

# === STEP 2: Find common introns across all SJ files ===
echo "Finding common introns across datasets..."

# Start with the first file intersected with introns
first_id=${SJ_IDS[0]}
bedtools intersect -a "$INTRON_BED" -b "temp/selected_${first_id}_SJ.txt" -wa -wb | \
  awk -v shift="$MAX_SHIFT" 'BEGIN{OFS="\t"} {d1=$8-$2; d2=$9-$3; if (d1>-shift && d1<shift && d2>-shift && d2<shift) print $1,$2,$3,$4,$5,$6}' | \
  sort -u > temp/common_introns.bed

# Iteratively intersect with remaining SJ files
for sid in "${SJ_IDS[@]:1}"; do
  bedtools intersect -a temp/common_introns.bed -b "temp/selected_${sid}_SJ.txt" -wa -wb | \
    awk -v shift="$MAX_SHIFT" 'BEGIN{OFS="\t"} {d1=$8-$2; d2=$9-$3; if (d1>-shift && d1<shift && d2>-shift && d2<shift) print $1,$2,$3,$4,$5,$6}' | \
    sort -u > temp/tmp.bed
  mv temp/tmp.bed temp/common_introns.bed
done

echo "Found $(wc -l < temp/common_introns.bed) shared introns."

# === STEP 3: Match gene IDs ===
#echo "Matching gene IDs..."
python "$GENE_MATCH_SCRIPT" "$EXPRESSED_GENES" temp/common_introns.bed temp/annotated_introns.txt

# Split by strand
awk '$6 == "+"' temp/annotated_introns.txt | sort -k1,1 -k2,2n > temp/introns_plus.bed
awk '$6 == "-"' temp/annotated_introns.txt | sort -k1,1 -k2,2n > temp/introns_minus.bed

# === STEP 3.5: Filter bedgraph files to keep signal > 0 ===
echo "Filtering bedgraph files (signal > 0 )..."
for sid in "${SJ_IDS[@]}"; do
  for strand in plus minus; do
    input_graph="pNET_${sid}_${strand}.bedgraph"
    output_graph="selected_${sid}_${strand}.bedgraph"

    if [[ -f "$input_graph" ]]; then
      awk '$4 > 0' "$input_graph" > "$output_graph"
    else
      echo "Warning: $input_graph not found. Skipping."
    fi
  done
done

# === STEP 4: Final filtering using bedgraph signal ===
echo "Filtering by expression signal..."
for sid in "${SJ_IDS[@]}"; do
  for strand in plus minus; do
    in_bed="temp/introns_${strand}.bed"
    graph_file="selected_${sid}_${strand}.bedgraph"
    temp_bed="temp/introns_${strand}_${sid}.bed"

    if [[ -f "$graph_file" ]]; then
      bedtools intersect -a "$in_bed" -b "$graph_file" -c > "$temp_bed"
    else
      echo "Warning: $graph_file not found. Skipping."
      continue
    fi
  done
done

# Keep introns with signal ≥ MIN_SIGNAL in all samples
for strand in plus minus; do
  in_bed="temp/introns_${strand}.bed"
  for sid in "${SJ_IDS[@]}"; do
    bedtools intersect -a "$in_bed" -b "selected_${sid}_${strand}.bedgraph" -c > "temp/tmp_${sid}_${strand}.bed"
    in_bed="temp/tmp_${sid}_${strand}.bed"
  done
  awk -v min="$MIN_SIGNAL" 'BEGIN{pass=1} {for (i=7;i<=NF;i++) if ($i<min) pass=0; if (pass) print $1,$2,$3,$4,$5,$6; pass=1}' \
    "$in_bed" | sort -k1,1 -k2,2n | sed 's/ /\t/g' > "Final_introns_${strand}.bed"
done

# === STEP 5: Prepare input for slope calculation ===
echo "Preparing data for slope calculation..."
for sid in "${SJ_IDS[@]}"; do
  for strand in plus minus; do
    if [[ -f "Final_introns_${strand}.bed" && -f "selected_${sid}_${strand}.bedgraph" ]]; then
      bedtools intersect -a "Final_introns_${strand}.bed" -b "selected_${sid}_${strand}.bedgraph" -wa -wb > "temp/${sid}_${strand}.txt"
    fi
  done
  cat temp/${sid}_plus.txt temp/${sid}_minus.txt | sort -k1,1 -k2,2n > "pNET_${sid}_points.txt"
  python "$SLOPE_SCRIPT" "pNET_${sid}_points.txt" > "pNET_${sid}_slope.txt"
done

echo "Done."
