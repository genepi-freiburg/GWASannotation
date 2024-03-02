#!/bin/bash

# Define input and output directories
input_file="$1" #input file for magma
output_path="$2" #path for output + prefix
bfile="$3"
#bfile="${3:-/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_150k_RandomSubset_Cleaned/hg38/UKBB_14k_hg38_chr1-22}"  # Default bfile if not provided #hg38
magma_annotated_genes="$4"

if [ ! -f "$magma_annotated_genes" ]; then
    echo "Creating MAGMA annotated genes"
    /data/programs/bin/gwas/magma/magma_v1.10/magma \
    --annotate --snp-loc "$bfile".bim \
    --gene-loc /data/programs/bin/gwas/pops_v0.2/magma_files/Ensembl_hg38.gene.loc --out "${output_path}_magma_annotate"
    magma_annotated_genes="${output_path}_magma_annotate.genes.annot"
else
    echo "Using existing file: $magma_annotated_genes"
fi

# Check exit status
if [ $? -ne 0 ]; then
    echo "Create MAGMA annotated genes fails."
    exit 1
fi

# Generate MAGMA scores
echo "Generate MAGMA score"
#use rsid or Chrpos input depending on bfile

# Check if bfile contains rsIDs
if grep -q "rs" "$bfile".bim; then
    echo "using _rsid.txt"
    input_file="${input_file}_rsid.txt"
else
    echo "using _CHRPOS.txt"
    input_file="${input_file}_CHRPOS.txt"
fi

out_magma="${output_path}_magma"
magma_command="/data/programs/bin/gwas/magma/magma-1.09b/magma"
$magma_command \
    --bfile $bfile \
    --gene-annot $magma_annotated_genes \
    --pval $input_file \
    ncol=N --gene-model snp-wise=mean --out $out_magma

# Check exit status
if [ $? -ne 0 ]; then
    echo "MAGMA command failed. Exiting."
    exit 1
fi

echo "Run  PoPS"
pops_command="python /data/programs/bin/gwas/pops_v0.2/pops.py"
bn=$(basename "$output_path")
out="${bn}"
$pops_command \
  --gene_annot_path /data/programs/bin/gwas/pops_v0.2/example/data/utils/gene_annot_hg38.txt \
  --feature_mat_prefix /data/programs/bin/gwas/gene_features/features_munged/pops_features \
  --num_feature_chunks 9 \
  --magma_prefix $out_magma \
  --out_prefix $output_path \
  --verbose

# Check exit status
if [ $? -ne 0 ]; then
    echo "PoPs command failed. Exiting."
    exit 1
fi

echo "Magma and  PoPS analysis finished"
