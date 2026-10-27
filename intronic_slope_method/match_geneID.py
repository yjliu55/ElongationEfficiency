import sys

def main():
    if len(sys.argv) != 4:
        print("Usage: python filter_genes.py gene_id_file bed_file output_file")
        sys.exit(1)

    gene_id_file, bed_file, output_file = sys.argv[1:4]

    # Read gene IDs into a set
    gene_ids = set()
    with open(gene_id_file, 'r') as f:
        for line in f:
            gene_id = line.strip()
            if gene_id:
                gene_ids.add(gene_id)

    # Filter BED file and write output
    with open(bed_file, 'r') as bed, open(output_file, 'w') as out:
        for line in bed:
            columns = line.strip().split('\t')
            if len(columns) >= 4:
                bed_gene_id = columns[3].strip()
                if bed_gene_id in gene_ids:
                    out.write(line)

if __name__ == "__main__":
    main()
