#!/bin/bash

# Define input and output directories
#input_file="$1" #input file for VEP
input_file="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/model1_qnorm_tkv.bsa_final_volume/GWAS_Annotation_MRI_38k_LD/model1_qnorm_tkv.bsa_proxies_vep_2.txt"
#output=${input_file/_vep/_vepout}
output="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/model1_qnorm_tkv.bsa_final_volume/GWAS_Annotation_MRI_38k_LD/model1_qnorm_tkv.bsa_proxies_vepout2_2.txt"

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
