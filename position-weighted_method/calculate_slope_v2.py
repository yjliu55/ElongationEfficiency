import sys
import numpy as np
from collections import defaultdict

def load_bedgraph(bedgraph_file):
    """
    Load bedgraph into memory and use numpy arrays for high-speed queries.
    """
    print(f"Loading {bedgraph_file} into memory...")
    bg_data = defaultdict(lambda: np.zeros(40000000, dtype=np.float32))
    
    with open(bedgraph_file, 'r') as f:
        for line in f:
            if line.strip() == "" or line.startswith(("track", "browser")):
                continue
            parts = line.strip().split()
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            count = float(parts[3])
            
            if count > 0:
                bg_data[chrom][start:end] = count
    return bg_data

def calculate_velocity(intron_bed, plus_bg, minus_bg, output_file):
    bg_plus_data = load_bedgraph(plus_bg)
    bg_minus_data = load_bedgraph(minus_bg)
    
    print("Calculating transcriptional velocities for perfectly intersected regions...")
    
    with open(intron_bed, 'r') as f_in, open(output_file, 'w') as f_out:
        header = ["chrom", "start", "end", "gene_id", "rank", "strand", 
                  "trimmed_start", "trimmed_end", "effective_length", "total_reads_X", "coverage", "velocity_V"]
        f_out.write("\t".join(header) + "\n")
        
        retained_count = 0
        total_count = 0
        
        for line in f_in:
            parts = line.strip().split()
            if len(parts) < 6:
                continue
                
            total_count += 1
            
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene_id = parts[3]
            rank = parts[4]
            strand = parts[5]
            
            # 1. Exclude 5bp from both left and right boundaries to avoid strong stall signal interference from the spliceosome
            new_start = start + 5
            new_end = end - 5
            N = new_end - new_start
            
            # If the length is insufficient for distribution calculation (<=1 bp) after excluding boundaries, discard it directly
            if N <= 1:
                continue
                
            # 2. Match plus and minus strands to get the base-level count array within the effective interval
            if strand == '+':
                counts = bg_plus_data[chrom][new_start:new_end]
            elif strand == '-':
                counts = bg_minus_data[chrom][new_start:new_end]
            else:
                continue
                
            # 3. Calculate total reads (X) and coverage
            X = np.sum(counts)
            coverage = X / N
            
            # 4. Generate the site step size matrix k (representing the travel distance from 5' to 3')
            if strand == '+':
                k_array = np.arange(N)
            else: 
                k_array = np.arange(N)[::-1]
                
            # 5. Substitute into the probability model formula (add denominator protection to prevent errors from all-zero signals)
            if X > 0:
                weights = (k_array / (N - 1)) * (100 / X)
                V = np.sum(counts * weights)
            else:
                V = 0.0
            
            # Limit extremely small out-of-bounds values that may be caused by floating-point operations
            V = min(max(V, 0.0), 100.0)
            retained_count += 1
            
            # Write results
            f_out.write(f"{chrom}\t{start}\t{end}\t{gene_id}\t{rank}\t{strand}\t"
                        f"{new_start}\t{new_end}\t{N}\t{X:.2f}\t{coverage:.4f}\t{V:.2f}\n")
                        
    print(f"Done! Calculated velocity for {retained_count}/{total_count} shared introns.")
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python calculate_slope_v9.py <intron_bed> <plus_bg> <minus_bg> <output_file>")
        sys.exit(1)
        
    INTRON_BED = sys.argv[1]
    PLUS_BG = sys.argv[2]
    MINUS_BG = sys.argv[3]
    OUTPUT_FILE = sys.argv[4]
    
    calculate_velocity(INTRON_BED, PLUS_BG, MINUS_BG, OUTPUT_FILE)
