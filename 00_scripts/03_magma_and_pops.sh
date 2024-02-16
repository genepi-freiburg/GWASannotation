#!/bin/bash

# Define input and output directories
input_file="$1" #input file for magma
output_path="$2" #path for output + prefix
bfile="${3:-/data/programs/bin/gwas/PoPS/data/1000G.EUR}"  # Default bfile if not provided

# Generate MAGMA scores
echo "Generate MAGMA score"
out_magma="${output_path}_magma"
magma_command="/data/programs/bin/gwas/magma/magma-1.09b/magma"
$magma_command \
    --bfile $bfile \
    --gene-annot /data/programs/bin/gwas/PoPS/data/magma_0kb.genes.annot \
    --pval $input_file \
    ncol=N --gene-model snp-wise=mean --out $out_magma

# Check exit status
if [ $? -ne 0 ]; then
    echo "Magma command failed. Exiting."
    exit 1
fi

# Run PoPS
echo "Run pops"
pops_command="python /data/programs/bin/gwas/pops_v0.2/pops.py"
bn=$(basename "$output_path")
out="${bn}"
$pops_command \
  --gene_annot_path /data/programs/bin/gwas/pops_v0.2/example/data/utils/gene_annot_jun10.txt \
  --feature_mat_prefix /data/programs/bin/gwas/gene_features/features_munged/pops_features \
  --num_feature_chunks 9 \
  --magma_prefix $out_magma \
  --out_prefix $out \
  --verbose

# Check exit status
if [ $? -ne 0 ]; then
    echo "PoPs command failed. Exiting."
    exit 1
fi

echo "Magma and pops analysis finished"
