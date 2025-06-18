import os
import glob
import sys

# ==== GET INPUT & OUTPUT PATHS ====
if len(sys.argv) != 3:
    print("Usage: python convert_quartet_to_bed.py <input_directory> <output_bed_file>")
    sys.exit(1)

input_dir = sys.argv[1]
output_file = sys.argv[2]

# ==== CHROMOSOME NAME MAPPING ====
chromosome_order = ["chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"]
chromosome_mapping = {}

# ==== PROCESS FILES ====
candidate_files = sorted(glob.glob(os.path.join(input_dir, "quarTeT.*.candidate")))
candidate_files = [file for file in candidate_files if "JAMFTS0" not in file]

with open(output_file, "w") as out_f:
    for file in candidate_files:
        print(f"Processing file: {file}")
        with open(file, "r") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 3 and '@TR_' not in parts[0] and not parts[0].startswith("#"):
                    original_chrom = parts[0]                 

                    if original_chrom not in chromosome_mapping:
                        if len(chromosome_mapping) < len(chromosome_order):
                            chromosome_mapping[original_chrom] = chromosome_order[len(chromosome_mapping)]
                    
                    mapped_chrom = chromosome_mapping[original_chrom]
                    start, end = parts[1], parts[2]
                    out_f.write(f"{mapped_chrom}\t{start}\t{end}\tquarTeT\n") 

print(f"Conversion complete! Output saved in: {output_file}")
