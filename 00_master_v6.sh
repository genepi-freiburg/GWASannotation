
#/opt/R/R-4.6.0-shared/bin/Rscript /dsk/epidata/programs/pipelines/GWASannotation/05_internal/00_process_sumstats_metal.R

#Test with imbi cluster
input_file="/dsk/epidata/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/09_metanalysis_Visconti2024/02_output/GWASAnnotation/HYT_N5H4S1/HYT_N5H4S1_hg38_dedup.RDS"
output_path="/dsk/epidata/programs/pipelines/GWASannotation/02_output/test_12June2026_HYT_N5H4S1_METAL"
sbatch --output /dsk/epidata/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_test_imbi_HYT_N5H3S1_METAL_qnorm_$(date +'%Y%m%d').txt --job-name=HYT_N5H4S1 --mem=40G --mail-type=FAIL --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="/opt/R/R-4.6.0-shared/bin/Rscript /dsk/epidata/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
--GWAS_RDS $input_file \
--bfile "/dsk/epidata/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/05_TopMed_2022/07_Plink/GCKD-TopMed_rsq03_maf001" \
--eQTL_datasets_coloc "GTEXv8" \
--pQTL_datasets_coloc "Icelanders_pGWAS" \
--output_path $output_path "
