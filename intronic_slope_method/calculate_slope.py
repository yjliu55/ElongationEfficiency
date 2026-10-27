#!/usr/bin/env python3
import sys
import numpy as np
from collections import defaultdict
from scipy.stats import linregress

BIN_SIZE = 10  # bp per bin

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py pNET_NH4OH_final_points.bed > output.txt")
        sys.exit(1)

    input_file = sys.argv[1]
    regions = defaultdict(list)

    # Read and parse the BED file
    with open(input_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) < 10:
                continue

            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            name, score, strand = fields[3], fields[4], fields[5]
            point_start, point_count = int(fields[7]), int(fields[9])

            intron_length = end - start
            if intron_length <= 0:
                continue

            # Position relative to 5′ end
            if strand == '+':
                rel_pos = point_start - start
            else:
                rel_pos = end - point_start

            key = (chrom, start, end, name, score, strand)
            regions[key].extend([rel_pos] * point_count)

    for region_key, positions in regions.items():
        chrom, start, end, name, score, strand = region_key
        intron_length = end - start

        if intron_length <= 0 or len(positions) < 2:
            print("\t".join([chrom, str(start), str(end), name, score, strand, "NA", "NA"]))
            continue

        total_counts = len(positions)
        num_bins = int(np.ceil(intron_length / BIN_SIZE))

        # Initialize bins
        bin_counts = np.zeros(num_bins)

        # Count reads in each bin
        for pos in positions:
            bin_idx = min(int(pos // BIN_SIZE), num_bins - 1)
            bin_counts[bin_idx] += 1

        # Cumulative counts from 5′ to 3′
        if strand == '+':
            cumulative = np.cumsum(bin_counts[::-1])[::-1]  # from 5' to 3'
        else:
            cumulative = np.cumsum(bin_counts)  # from 3' to 5'

        # Normalize cumulative counts
        y_values = cumulative / total_counts

        # Normalized x values: bin centers scaled to intron length (0 to 0.1)
        bin_centers = (np.arange(num_bins) * BIN_SIZE + BIN_SIZE / 2) / (intron_length * 0.1)
        x_values = bin_centers

        # Reverse x for negative strand
        if strand == '-':
            x_values = x_values[::-1]
            y_values = y_values

        # Linear regression
        try:
            slope_val, _, r_value, _, _ = linregress(x_values, y_values)
            slope = f"{slope_val:.4f}"
            r_squared = f"{r_value**2:.4f}"

            # Calculate Pol II speed: Log10(-1 / slope_val)
            if slope_val < 0:
                polii_speed = np.log2(-1 / slope_val)
                polii_speed = f"{polii_speed:.4f}"
            else:
                polii_speed = "NA"
        except:
            slope, r_squared, polii_speed = "NA", "NA", "NA"

        print("\t".join([chrom, str(start), str(end), name, score, strand, slope, r_squared, polii_speed]))

if __name__ == "__main__":
    main()
