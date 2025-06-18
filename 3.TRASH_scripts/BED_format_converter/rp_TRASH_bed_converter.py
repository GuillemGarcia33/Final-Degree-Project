import sys
import pandas as pd

# Get arguments
if len(sys.argv) != 3:
    print("Usage: python3 rp_TRASH_converter.py <input_csv> <output_bed>")
    sys.exit(1)

csv_file = sys.argv[1]
output_bed = sys.argv[2]

# Mapping to rename scaffolds to Chr names
chromosome_map = {
    "CM051459.1": "Chr01",
    "CM051460.1": "Chr02",
    "CM051461.1": "Chr03",
    "CM051462.1": "Chr04",
    "CM051463.1": "Chr05"
}

# Read CSV
df = pd.read_csv(csv_file)
df = df[df['fasta.name'].isin(chromosome_map.keys())]
df['fasta.name'] = df['fasta.name'].map(chromosome_map)

# Create BED DataFrame
bed_df = df[['fasta.name', 'start', 'end']].copy()
bed_df['label'] = 'TRASH'

# Save as BED
bed_df.to_csv(output_bed, sep='\t', header=False, index=False)
print(f"Generated BED file with {len(bed_df)} entries: {output_bed}")
