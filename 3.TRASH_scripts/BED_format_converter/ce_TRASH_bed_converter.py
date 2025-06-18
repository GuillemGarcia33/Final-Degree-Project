import sys
import pandas as pd

# ==== CHECK ARGUMENTS ====
if len(sys.argv) != 3:
    print("Usage: python3 convert_trash_to_bed_chrI_to_chrX.py <input_csv> <output_bed>")
    sys.exit(1)

input_csv = sys.argv[1]
output_bed = sys.argv[2]

# ==== VALID CHROMOSOMES ====
valid_chromosomes = {
    "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrX"
}

# ==== READ CSV ====
df = pd.read_csv(input_csv)

# ==== FILTER FOR VALID CHROMOSOMES ====
df = df[df['fasta.name'].isin(valid_chromosomes)]

# ==== CREATE BED DATA ====
bed_df = df[['fasta.name', 'start', 'end']].copy()
bed_df['label'] = 'TRASH'

# ==== SAVE TO BED ====
bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
print(f"Generated BED file with {len(bed_df)} entries: {output_bed}")
