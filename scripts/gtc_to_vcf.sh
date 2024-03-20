#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 9 ]; then
    echo "Usage: $0 <manifest_bpm> <manifest_csv> <clusterfile> <gtcs-dir> <ref_genome> <stats_tsv> <vcf_gz_name> <n_threads> <csi_name>"
    exit 1
fi

# Assign command line arguments to variables
manifest_bpm="$1"
manifest_csv="$2"
clusterfile="$3"
gtcs_dir="$4"
ref_genome="$5"
stats_tsv="$6"
vcf_gz_name="$7"
n_threads="$8"
csi_name="$9"

bcftools +gtc2vcf \
  --no-version -Ou \
  --bpm "$manifest_bpm" \
  --csv "$manifest_csv" \
  --egt "$clusterfile" \
  --gtcs "$gtcs_dir" \
  --fasta-ref "$ref_genome" \
  --extra "$stats_tsv" | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --no-version -Oz -c x -f "$ref_genome" | \
  tee "$vcf_gz_name" | \
  bcftools index --threads "$n_threads" --force --output "$csi_name"

echo "Converting GTCs to VCF completed"