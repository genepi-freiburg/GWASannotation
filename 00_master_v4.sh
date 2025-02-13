#September 2024
sbatch --output  /data/programs/pipelines/GWASannotation/dev/GWASannotation/03_logs/test.txt --job-name=Nora --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/dev/GWASannotation/00_scripts/test.R"

#using module GWAS as test
mod_name="ME30"
output_path="/data/programs/pipelines/GWASannotation/dev/GWASannotation/02_output/${mod_name}"
input_file="/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/30_M30_GWASAnnotation/02_output_v2/${mod_name}/${mod_name}_liftOver_hg38_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-213" --output /data/programs/pipelines/GWASannotation/dev/GWASannotation/03_logs/GWASAnno_${mod_name}_$(date +'%Y%m%d').txt --job-name=Anno${mod_name} --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/dev/GWASannotation/00_scripts/GWASAnno_main.R \
--GWAS_RDS $input_file \
--bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/02_HRC/Plink_R2_0.3_MAF_0.01/hg38/GCKD_HRC_R2_03_MAF_001_hg38_chr1-22" \
--eQTL_datasets_coloc "Kidney_eQTL,GTEXv8,eQTLGen" \
--eQTL_tissues_interest_coloc "Adrenal_Gland,Heart_Left_Ventricle,Kidney_Cortex,Liver,Pancreas,Small_Intestine_Terminal_Ileum,Colon_Transverse,Stomach,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" \
--output_path $output_path "

#using module GWAS as test with coloc input path
mod_name="ME30"
output_path="/data/programs/pipelines/GWASannotation/dev/GWASannotation/02_output/${mod_name}"
input_file="/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/30_M30_GWASAnnotation/02_output_v2/${mod_name}/${mod_name}_liftOver_hg38_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-213" --output /data/programs/pipelines/GWASannotation/dev/GWASannotation/03_logs/GWASAnno_${mod_name}_with_colocinput_$(date +'%Y%m%d').txt --job-name=Anno${mod_name} --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/dev/GWASannotation/00_scripts/GWASAnno_main.R \
--GWAS_RDS $input_file \
--bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/02_HRC/Plink_R2_0.3_MAF_0.01/hg38/GCKD_HRC_R2_03_MAF_001_hg38_chr1-22" \
--eQTL_datasets_coloc "Kidney_eQTL,GTEXv8,eQTLGen" \
--eQTL_tissues_interest_coloc "Adrenal_Gland,Heart_Left_Ventricle,Kidney_Cortex,Liver,Pancreas,Small_Intestine_Terminal_Ileum,Colon_Transverse,Stomach,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" \
--coloc_input_path "/data/programs/pipelines/GWASannotation/dev/GWASannotation/02_output/coloc/output/" \
--output_path $output_path "


output_path="/data/programs/pipelines/GWASannotation/dev/GWASannotation/02_output/${mod_name}"
sbatch --output /data/programs/pipelines/GWASannotation/dev/GWASannotation/03_logs/06_postprocessing_$(date +'%Y%m%d').txt --job-name=post --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/dev/GWASannotation/00_scripts/06_postprocessing.R \
    --output_path $output_path "

#09Feb2025 -using /data/programs/pipelines/GWASannotation_v2 for dev
v="medulla"
mkdir /data/programs/pipelines/GWASannotation/02_output/model1_qnorm_${v}.bsa_final_volume/
output_path="/data/programs/pipelines/GWASannotation/02_output/model1_qnorm_${v}.bsa_final_volume/model1_qnorm_${v}.bsa"
input_file="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/10Feb2024_final_volumes/maf001/model1_qnorm_${v}.bsa_chr1-22_maf001_liftOver_hg38_dedup.RDS"
sbatch -p TempCompute --nodelist="imbi12" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_MRI_38k_LD_main_window_1000kb_model1_qnorm_${v}.bsa_dev_$(date +'%Y%m%d').txt --job-name=Anno${v} --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="/scratch/global/martins/R-4.4.1/bin/Rscript /data/programs/pipelines/GWASannotation_dev/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --interval_window_kb 1000 \
    --bfile "/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_MRI_38k_MAF001/plink/chrpos/hg38/UKBB_MRI_38k_MAF001_hg38_allchr_v2" \
    --eQTL_datasets_coloc "GTEXv8" \
    --eQTL_tissues_interest_coloc "Kidney_Cortex,Liver,Whole_Blood" \
    --output_path $output_path "


