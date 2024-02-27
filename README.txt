#31Jan2024
GOAL: create an annotation pipeline for GWAS loci, by editing/ upgrading the ProGeM pipline.
- botom up: add coloc with eQTL and pQTLs
- top dow: add PoPS 
optparse was added to scripts/ProGeM_settings.r to facilitate using the pipeline in parallel for different tratis.

<<<<<<< HEAD
INPUT:
- sumstats RDS file in hg38 build with the following columns 

| Column | Description |
| --- | --- |
| Name | CHR:POS:A2:A1 (example: *chr1:10177:A:AC*) |
| rsID |  Reference SNP cluster ID (example *rs367896724*) |
| CHR | chromosome number (example: *1*) |
| POS | SNP position (example: *10177*) |
| A1 | effect allele (example: *AC*) |
| A2 | other allele (example: *A*) |
| BETA | beta |
| SE | se|
| nlog10P | -log10(P) |
| AF | frequency of A1 |
| N | sample size |
*script 00_process_sumstats.R can be used to create an input file from a gz and tabix regenie output
* analysis will be conducted in hg38!!! 

- output path: full output path (including the basename of the output files)

=======
Usage example:
gwas="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/20Oct/maf001/model1_qnorm_tkv.bsa_chr1-22_maf001.regenie.gz"
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
#mkdir /data/programs/pipelines/GWASannotation/02_output/test
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/01_preprocessing_$(date +'%Y%m%d').txt --job-name=preprocess_${v} --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/01_preprocessing.R \
--GWAS $gwas \
--output_path $output_path"

input_file="${output_path}_proxies_vep.txt"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/magma_and_pops_$(date +'%Y%m%d').txt --job-name=vep --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de /data/programs/pipelines/GWASannotation/00_scripts/02_vep.sh $input_file

input_file="${output_path}_input_magma.txt"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/03_magma_and_pops_$(date +'%Y%m%d').txt --job-name=magma_pops --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de /data/programs/pipelines/GWASannotation/00_scripts/03_magma_and_pops.sh $input_file $output_path

sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/ProGeM_$(date +'%Y%m%d').txt --job-name=ProGeM --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/ProGeM_settings.r \
        --sentinel_filepath /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/01_VEP/02_output/model2_qnorm_${v}.bsa_cov_eGFR/model2_qnorm_${v}.bsa_cov_eGFR_sentinel.txt \
        --proxy_filepath /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/01_VEP/02_output/model2_qnorm_${v}.bsa_cov_eGFR/model2_qnorm_${v}.bsa_cov_eGFR_proxies.txt \
        --VEP_filepath /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/01_VEP/02_output/model2_qnorm_${v}.bsa_cov_eGFR/model2_qnorm_${v}.bsa_cov_eGFR_proxies_vepout.txt \
        --PoPS_filepath /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/02_Pops/02_output/model2_qnorm_${v}.bsa_cov_eGFR/model2_qnorm_${v}.bsa_cov_eGFR.RData \
        --output_folder $output_folder"
>>>>>>> 3f93e05e8aa5039d4ef16e1e377837909f4a5f4b
