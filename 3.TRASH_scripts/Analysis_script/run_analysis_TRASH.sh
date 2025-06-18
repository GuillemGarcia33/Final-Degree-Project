#!/bin/bash

# ==== CONFIGURATION ====
INPUT_DIR="/home/guillem/Descargas/Internship/Genomes/B.mori/TRASH_RESULTS/B_mori_TRASH/B_mori/Summary.of.repetitive.regions.Bomo_genome_assembly.fa.csv"          # Directory containing TRASH CSV summary files
REFERENCE_BED="/home/guillem/Descargas/Internship/Genomes/B.mori/bmori/cenptdomains_5kbremoved.bed"           # Reference BED file path
OUTPUT_DIR="/home/guillem/Descargas/Internship/Genomes/B.mori/TRASH_RESULTS/New_Results"           # Output directory path
PYTHON_SCRIPT="/home/guillem/Descargas/Internship/Scripts/bm_TRASH_converter.py"     # Python script to convert TRASH CSV to BED

# ==== EXECUTION ====
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/BED_FILES"
mkdir -p "$OUTPUT_DIR/BED_CLUSTERED_CAND_REF"
mkdir -p "$OUTPUT_DIR/METRICS"
mkdir -p "$OUTPUT_DIR/METRICS/CLUSTER"

# Step 1: Convert TRASH CSV summary into BED format
echo "Converting TRASH summary CSV into BED format..."
python3 "$PYTHON_SCRIPT" "$INPUT_DIR" "$OUTPUT_DIR/BED_FILES/trash_predictions.bed"

# Step 2: Format reference BED file
echo "Formatting reference BED file..."
awk '{print $1, $2, $3, "reference"}' OFS='\t' "$REFERENCE_BED" > "$OUTPUT_DIR/BED_FILES/reference_formatted.bed"

# Step 3: Merge and sort both BED files
echo "Merging and sorting TRASH predictions with reference BED file..."
cat "$OUTPUT_DIR/BED_FILES/reference_formatted.bed" "$OUTPUT_DIR/BED_FILES/trash_predictions.bed" | sort -k1,1 -k2,2n > "$OUTPUT_DIR/BED_FILES/merged_sorted.bed"

# Step 4: Cluster merged BED file
echo "Clustering the merged BED file..."
bedtools cluster -i "$OUTPUT_DIR/BED_FILES/merged_sorted.bed" > "$OUTPUT_DIR/BED_CLUSTERED_CAND_REF/merged_clustered.bed"

# Step 5: Compute reference-centric metrics
python3 - << EOF > "$OUTPUT_DIR/METRICS/CLUSTER/cluster_metrics_trash.txt"
import pandas as pd
import math

file_path = "$OUTPUT_DIR/BED_CLUSTERED_CAND_REF/merged_clustered.bed"
data = pd.read_csv(file_path, sep="\t", header=None, names=["chr", "start", "end", "source", "cluster"])

total_TP = 0
total_FP = 0
total_FN = 0
total_TP_star = 0
cluster_logs = []

for cluster_id, group in data.groupby("cluster"):
    references = group[group["source"] == "reference"]
    predictions = group[group["source"] == "TRASH"]
    cluster_log = [f"Cluster {cluster_id}"]

    TP = 0
    FP = 0
    FN = 0
    TP_star = 0

    if len(references) == 1 and len(predictions) == 1:
        TP_star += 1
        TP += 1
        ref = references.iloc[0]
        pred = predictions.iloc[0]
        cluster_log.append(f"  reference {ref['start']}-{ref['end']} ≈ trash {pred['start']}-{pred['end']} [TP*]")
    else:
        for _, ref in references.iterrows():
            if not predictions.empty:
                TP += 1
                cluster_log.append(f"  reference {ref['start']}-{ref['end']} [TP]")
            else:
                FN += 1
                cluster_log.append(f"  reference {ref['start']}-{ref['end']} [FN]")

        if references.empty:
            for _, pred in predictions.iterrows():
                FP += 1
                cluster_log.append(f"  trash {pred['start']}-{pred['end']} [FP]")

    total_TP += TP
    total_FN += FN
    total_FP += FP
    total_TP_star += TP_star

    labels = []
    if TP: labels.append(f"TP x {TP}")
    if TP_star: labels.append(f"TP* x {TP_star}")
    if FP: labels.append(f"FP x {FP}")
    if FN: labels.append(f"FN x {FN}")
    cluster_log.insert(1, f"  Summary: [{', '.join(labels)}]")
    cluster_logs.extend(cluster_log)

# Global metrics
sensitivity = total_TP / (total_TP + total_FN) if (total_TP + total_FN) else 0
precision = total_TP / (total_TP + total_FP) if (total_TP + total_FP) else 0
f1_score = (2 * precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) else 0
accg = math.sqrt(sensitivity * precision) if (precision and sensitivity) else 0

sensitivity_star = total_TP_star / (total_TP + total_FN) if (total_TP + total_FN) else 0
precision_star = total_TP_star / (total_TP + total_FP) if (total_TP + total_FP) else 0
f1_score_star = (2 * precision_star * sensitivity_star) / (precision_star + sensitivity_star) if (precision_star + sensitivity_star) else 0
accg_star = math.sqrt(sensitivity_star * precision_star) if (precision_star and sensitivity_star) else 0

# Output
print("Reference-centric Cluster-based Metrics for TRASH")
print(f"TP: {total_TP}")
print(f"TP* (1:1 cluster match): {total_TP_star}")
print(f"FP: {total_FP}")
print(f"FN: {total_FN}")
print(f"Sensitivity: {sensitivity:.6f}")
print(f"Precision: {precision:.6f}")
print(f"F1 Score: {f1_score:.6f}")
print(f"Geometric Accuracy (ACCG): {accg:.6f}")
print(f"Strict Sensitivity* (TP*): {sensitivity_star:.6f}")
print(f"Strict Precision* (TP*): {precision_star:.6f}")
print(f"Strict F1 Score* (TP*): {f1_score_star:.6f}")
print(f"Strict Geometric Accuracy* (ACCG*): {accg_star:.6f}")
print("\nCluster Details:")
print("\n".join(cluster_logs))
EOF

echo "Analysis complete! Results saved in $OUTPUT_DIR/METRICS/CLUSTER/"

