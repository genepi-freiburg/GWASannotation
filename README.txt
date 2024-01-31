#31Jan2024
GOAL: create an annotation pipeline for GWAS loci, by editing/ upgrading the ProGeM pipline.
- botom up: keep it as it is for now
- top dow: add PoPS and coloc with eQTL and pQTLs
optparse was added to scripts/ProGeM_settings.r to facilitate using the pipeline in parallel for different tratis.

Usage example:
 sbatch -p imbiPCompute --nodelist="imbip-compute-213" --output /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/03_logs/ProGeM_model2_qnorm_${v}.bsa_cov_eGFR_Jan2024.txt --job-name=ProGeM --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/00_scripts/ProGeM_settings.r \
        --sentinel_filepath /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/01_VEP/02_output/model2_qnorm_${v}.bsa_cov_eGFR/model2_qnorm_${v}.bsa_cov_eGFR_sentinel.txt \
        --proxy_filepath /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/01_VEP/02_output/model2_qnorm_${v}.bsa_cov_eGFR/model2_qnorm_${v}.bsa_cov_eGFR_proxies.txt \
        --VEP_filepath /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/01_VEP/02_output/model2_qnorm_${v}.bsa_cov_eGFR/model2_qnorm_${v}.bsa_cov_eGFR_proxies_vepout.txt \
        --PoPS_filepath /data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/02_Pops/02_output/model2_qnorm_${v}.bsa_cov_eGFR/model2_qnorm_${v}.bsa_cov_eGFR.RData \
        --output_folder $output_folder"
