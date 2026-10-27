#!/bin/bash
set -euo pipefail

# === CONFIG ===
MIN_READ_SUPPORT=3
MAX_SHIFT=5
MIN_LOCI_COVERED=5  
INTRON_BED="TAIR10_pure_intron.bed"
EXPRESSED_GENES="genes.txt"
GENE_MATCH_SCRIPT="match_geneID.py"
SLOPE_SCRIPT="calculate_slope_v2.py"

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
  base=$(basename "$base" SJ.out.tab)
  
  SJ_IDS+=("$base")
  awk -v min="$MIN_READ_SUPPORT" '($7 >= min)' "$sj_file" > "temp/selected_${base}_SJ.txt"
done

# === STEP 2: Find common introns across all SJ files ===
echo "Finding common introns across datasets..."

first_id=${SJ_IDS[0]}
bedtools intersect -a "$INTRON_BED" -b "temp/selected_${first_id}_SJ.txt" -wa -wb | \
  awk -v shift="$MAX_SHIFT" 'BEGIN{OFS="\t"} {d1=$8-$2; d2=$9-$3; if (d1>-shift && d1<shift && d2>-shift && d2<shift) print $1,$2,$3,$4,$5,$6}' | \
  sort -u > temp/common_introns.bed

for sid in "${SJ_IDS[@]:1}"; do
  bedtools intersect -a temp/common_introns.bed -b "temp/selected_${sid}_SJ.txt" -wa -wb | \
    awk -v shift="$MAX_SHIFT" 'BEGIN{OFS="\t"} {d1=$8-$2; d2=$9-$3; if (d1>-shift && d1<shift && d2>-shift && d2<shift) print $1,$2,$3,$4,$5,$6}' | \
    sort -u > temp/tmp.bed
  mv temp/tmp.bed temp/common_introns.bed
done

echo "Found $(wc -l < temp/common_introns.bed) shared introns."

# === STEP 3: Match gene IDs ===
python "$GENE_MATCH_SCRIPT" "$EXPRESSED_GENES" temp/common_introns.bed temp/annotated_introns.txt

awk '$6 == "+"' temp/annotated_introns.txt | sort -k1,1 -k2,2n > temp/introns_plus.bed
awk '$6 == "-"' temp/annotated_introns.txt | sort -k1,1 -k2,2n > temp/introns_minus.bed

# === STEP 4: Filter bedgraph files to keep signal > 0 ===
echo "Filtering normalized bedgraph files (signal > 0)..."
for sid in "${SJ_IDS[@]}"; do
  for strand in plus minus; do
    input_graph="pNET_${sid}_${strand}.bedgraph"
    output_graph="selected_${sid}_${strand}.bedgraph"

    if [[ -f "$input_graph" ]]; then
      awk '$4 > 0' "$input_graph" > "$output_graph"
    fi
  done
done

# === STEP 5: Final filtering using Coverage Breadth (Loci Count) ===
echo "Filtering by coverage breadth (Introns must have >= ${MIN_LOCI_COVERED} covered loci)..."
for sid in "${SJ_IDS[@]}"; do
  for strand in plus minus; do
    in_bed="temp/introns_${strand}.bed"
    graph_file="selected_${sid}_${strand}.bedgraph"
    temp_bed="temp/introns_${strand}_${sid}.bed"

    if [[ -f "$graph_file" ]]; then
      bedtools intersect -a "$in_bed" -b "$graph_file" -c > "$temp_bed"
    fi
  done
done

for strand in plus minus; do
  in_bed="temp/introns_${strand}.bed"
  for sid in "${SJ_IDS[@]}"; do
    bedtools intersect -a "$in_bed" -b "selected_${sid}_${strand}.bedgraph" -c > "temp/tmp_${sid}_${strand}.bed"
    in_bed="temp/tmp_${sid}_${strand}.bed"
  done
  
  awk -v min="$MIN_LOCI_COVERED" 'BEGIN{pass=1} {for (i=7;i<=NF;i++) if ($i<min) pass=0; if (pass) print $1,$2,$3,$4,$5,$6; pass=1}' \
    "$in_bed" | sort -k1,1 -k2,2n | sed 's/ /\t/g' > "temp/Final_introns_${strand}.bed"
done

# === STEP 6: Calculate trimmed coverage & filter out the bottom 5% per sample ===
echo "Generating internal Python script for 5% coverage filtering..."
cat << 'EOF' > temp/filter_cov_5percent.py
import sys
import numpy as np
from collections import defaultdict

def load_bedgraph(bg_file):
    bg_data = defaultdict(lambda: np.zeros(40000000, dtype=np.float32))
    try:
        with open(bg_file, 'r') as f:
            for line in f:
                if line.strip() == "" or line.startswith(("track", "browser")): continue
                parts = line.strip().split()
                chrom, start, end, count = parts[0], int(parts[1]), int(parts[2]), float(parts[3])
                if count > 0: bg_data[chrom][start:end] = count
    except FileNotFoundError:
        pass
    return bg_data

intron_bed, plus_bg, minus_bg, out_bed = sys.argv[1:5]
plus_data, minus_data = load_bedgraph(plus_bg), load_bedgraph(minus_bg)

introns = []
with open(intron_bed, 'r') as f:
    for line in f:
        parts = line.strip().split()
        if len(parts) < 6: continue
        chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]
        new_start, new_end = start + 5, end - 5
        N = new_end - new_start
        cov = 0.0
        if N > 1:
            if strand == '+': cov = np.sum(plus_data[chrom][new_start:new_end]) / N
            elif strand == '-': cov = np.sum(minus_data[chrom][new_start:new_end]) / N
        introns.append((line.strip(), cov))

# Sort by coverage and drop the bottom 5%
introns.sort(key=lambda x: x[1])
cutoff_idx = int(len(introns) * 0.05)
passed_introns = introns[cutoff_idx:]

with open(out_bed, 'w') as f:
    for item in passed_introns:
        f.write(item[0] + '\n')
EOF

cat temp/Final_introns_plus.bed temp/Final_introns_minus.bed | sort -k1,1 -k2,2n > temp/Final_introns_combined.bed

echo "Filtering bottom 5% lowest coverage introns (trimmed 5bp) per sample..."
for sid in "${SJ_IDS[@]}"; do
  plus_bg="selected_${sid}_plus.bedgraph"
  minus_bg="selected_${sid}_minus.bedgraph"
  if [[ -f "$plus_bg" && -f "$minus_bg" ]]; then
    python temp/filter_cov_5percent.py temp/Final_introns_combined.bed "$plus_bg" "$minus_bg" "temp/passed_introns_${sid}.bed"
  else
    echo "Warning: Bedgraph missing for ${sid}. Proceeding with unfiltered combined introns."
    cp temp/Final_introns_combined.bed "temp/passed_introns_${sid}.bed"
  fi
done

# === STEP 7: Intersect passed introns across ALL samples ===
echo "Intersecting retained introns across all samples to unify input regions..."
cp temp/Final_introns_combined.bed temp/Final_introns_intersected.bed

for sid in "${SJ_IDS[@]}"; do
  awk 'NR==FNR {a[$0]=1; next} {if(a[$0]) print $0}' "temp/passed_introns_${sid}.bed" temp/Final_introns_intersected.bed > temp/tmp_intersect.bed
  mv temp/tmp_intersect.bed temp/Final_introns_intersected.bed
done

echo "Final identical regions retained for velocity calculation: $(wc -l < temp/Final_introns_intersected.bed)"

# === STEP 8: Prepare data and calculate elongation efficiency (PEE) ===
echo "Computing emergent elongation speeds (Probability Distribution Model)..."

for sid in "${SJ_IDS[@]}"; do
  in_bed="temp/Final_introns_intersected.bed"
  bg_plus="selected_${sid}_plus.bedgraph"
  bg_minus="selected_${sid}_minus.bedgraph"
  speed_out="pNET_${sid}_speed.txt"
  
  if [[ -f "$bg_plus" && -f "$bg_minus" ]]; then
    python "$SLOPE_SCRIPT" "$in_bed" "$bg_plus" "$bg_minus" "$speed_out"
  else
    echo "Warning: Missing bedgraph files for sample ${sid}, skipping speed calculation."
  fi
done

# === STEP 9: Data Aggregation for Plotting ===
echo "Aggregating speed values into a single plotting matrix..."
PLOT_OUTPUT="Combined_Elongation_Speeds.tsv"
echo -e "Sample\tIntron_ID\tSpeed" > "$PLOT_OUTPUT"

for sid in "${SJ_IDS[@]}"; do
  speed_file="pNET_${sid}_speed.txt"
  if [[ -f "$speed_file" ]]; then
    awk -v sample="$sid" 'BEGIN{OFS="\t"} NR>1 {print sample, $4, $12}' "$speed_file" >> "$PLOT_OUTPUT"
  fi
done

echo "Done. Final aggregated data saved to '$PLOT_OUTPUT'."
