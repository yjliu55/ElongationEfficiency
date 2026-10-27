# Intron-position-weighted method to calculate Pol II Elongation Efficiency (PEE)

Pipeline to calculate RNA Polymerase II (Pol II) transcription elongation efficiency from pNET-seq data using an upgraded probability distribution model.

$$V = \sum_{i=1}^{N} \left( C_i \times \frac{k_i}{N - 1} \times \frac{100}{X} \right)$$

Where:
- $N$: Effective intron length (trimmed by 5bp on both boundaries).
- $C_i$: Signal count at position $i$.
- $k_i$: Distance index from 5' to 3' ($0, 1, \dots, N-1$).
- $X$: Total reads across the effective region ($\sum C_i$).

---

## 🛠 Prerequisites

- **Python 3** (`numpy`)
- **bedtools**

---

## 📁 File Structure

- `modify_temp_v2.sh` - Main workflow execution script.
- `calculate_slope_v2.py` - Core calculation script implementing the probability distribution model.
- `match_geneID.py` - Script for matching introns with gene annotations.
- `TAIR10_pure_intron.bed` - Reference intron annotation file.
- `genes.txt` - Target/expressed gene list.

### Example Data (`example/`)
- `WT_22_rep1_SJ.out.tab` - Splice junction file from STAR.
- `pNET_WT_22_rep1_plus.bedgraph` - Normalized pNET-seq signal file for the plus strand.
- `pNET_WT_22_rep1_minus.bedgraph` - Normalized pNET-seq signal file for the minus strand.

---

## 🚀 Usage

Run the pipeline by specifying one or more `*_SJ.out.tab` files:

```bash
bash modify_temp_v2.sh sample_rep1_SJ.out.tab
```

For multiple replicates/samples:

```bash
bash modify_temp_v2.sh sample_rep1_SJ.out.tab sample_rep2_SJ.out.tab ...
```
> **Note:** Ensure that the corresponding `pNET_<ID>_plus.bedgraph` and `pNET_<ID>_minus.bedgraph` files are located in the working directory (where `<ID>` matches the suffix prefix before `_SJ.out.tab`).

---

## ⚙️ Workflow Steps

1. **Junction Filtering:** Filters splice junctions from STAR outputs by read support (`MIN_READ_SUPPORT=3`).
2. **Common Intron Identification:** Intersects and retains introns shared across all input samples.
3. **Gene Annotation:** Matches common introns to expressed genes (`genes.txt`).
4. **Boundary Trimming:** Excludes 5 bp from both 5' and 3' boundaries of each intron to avoid spliceosome stalling artifacts.
5. **Coverage Filtering:** Drops the bottom 5% lowest-coverage introns per sample and unifies the remaining common set across samples.
6. **Velocity Calculation::** Evaluates transcription velocity ($V$) for each retained intron.
7. **Aggregation:** Combines output values across samples into a single matrix for downstream plotting.

---

## 📊 Output

### 1. The script outputs `pNET_<ID>_speed.txt` for each sample with the following tab-separated columns:
1. `chrom` - Chromosome name
2. `start` - Intron start position
3. `end` - Intron end position
4. `gene_id` - Gene identifier
5. `rank` - Intron rank/index
6. `strand` - Genomic strand (`+` / `-`)
7. `trimmed_start` - Boundary-trimmed start position (+5 bp)
8. `trimmed_end` - Boundary-trimmed end position (-5 bp)
9. `effective_length` - Effective intron length ($N$)
10. `total_reads_X` - Total signal read count within effective region ($X$)
11. `coverage` - Average signal coverage ($X / N$)
12. `velocity_V` - Calculated transcription velocity score ($V \in [0, 100]$)
