import sys
import pandas as pd
import re

# ==== CHECK ARGUMENTS ====
if len(sys.argv) != 3:
    print("Usage: python3 convert_trash_to_bed_bombyx.py <input_csv> <output_bed>")
    sys.exit(1)

input_csv = sys.argv[1]
output_bed = sys.argv[2]

# ==== READ CSV ====
df = pd.read_csv(input_csv)

# ==== FILTER CHROMOSOMES MATCHING 'Bomo_Chr' + NUMBER ====
pattern = re.compile(r'^Bomo_Chr(\d+)$')
df = df[df['fasta.name'].apply(lambda x: bool(pattern.match(str(x))))]

# ==== RENAME CHROMOSOMES TO 'chrX' ====
df['fasta.name'] = df['fasta.name'].apply(lambda x: 'chr' + pattern.match(x).group(1))

# ==== CREATE BED DATA ====
bed_df = df[['fasta.name', 'start', 'end']].copy()
bed_df['label'] = 'TRASH'

# ==== SAVE TO BED ====
bed_df.to_csv(output_bed, sep='\t', header=False, index=False)

print(f"Generated BED file with {len(bed_df)} entries: {output_bed}")
