#!/bin/bash

# Define input and output directories
input_file="$1" #input file for VEP
output=${input_file/_vep/_vepout}

# Generate MAGMA scores
echo "Run VEP"
/data/programs/bin/ngs/VEP/101/vep \
--input_file $input_file --format ensembl --no_check_variants_order \
--output_file $output \
--tab \
--force_overwrite \
--cache --dir /data/programs/bin/ngs/VEP/cache \
--canonical --biotype --symbol \
--verbose \
--af_gnomad \
--regulatory \
--offline \
--assembly GRCh38 \
--verbose \
--fork 8

# Check exit status
if [ $? -ne 0 ]; then
    echo "VEP command failed. Exiting."
    exit 1
fi

echo "VEP finished"
