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
