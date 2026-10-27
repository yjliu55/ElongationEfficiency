# Slope-based method to calculate Pol II elongation speed

Pipeline to calculate RNA Polymerase II (Pol II) transcription rate from pNET-seq data using the intronic slope method.

$$\text{Pol II Speed} = \log_2\left(-\frac{1}{\text{slope}}\right)$$

---

## Prerequisites

- **Python 3** (`numpy`, `scipy`)
- **bedtools**

---

## File Structure

- `modify_temp.sh` - Main workflow execution script.
- `calculate_slope.py` - Python script for calculating intronic slope and Pol II speed.
- `match_geneID.py` - Script for matching introns with gene annotations.
- `TAIR10_pure_intron.bed` - Reference intron annotation file.
- `genes.txt` - List of target/expressed genes.

### Example Data
- `WT_22_rep1_SJ.out.tab` - Splice junction file (STAR output).
- `pNET_WT_22_rep1_plus.bedgraph` - pNET-seq signal file for the plus strand.
- `pNET_WT_22_rep1_minus.bedgraph` - pNET-seq signal file for the minus strand.

---

## Usage

Run the pipeline by passing one or more `*_SJ.out.tab` files as arguments:

```bash
bash modify_temp_v2.sh sample_rep1_SJ.out.tab
```

For multiple replicates/samples:

```bash
bash modify_temp_v2.sh sample_rep1_SJ.out.tab sample_rep2_SJ.out.tab ...
```

> **Note:** Ensure that the corresponding `pNET_<ID>_plus.bedgraph` and `pNET_<ID>_minus.bedgraph` files are located in the working directory (where `<ID>` matches the suffix prefix before `_SJ.out.tab`).

---

## Output

The script outputs `pNET_<ID>_slope.txt` for each sample with the following tab-separated columns:
1. `Chrom`
2. `Start`
3. `End`
4. `Gene_ID`
5. `Score`
6. `Strand`
7. `Slope`
8. `R_squared`
9. `PolII_Speed`