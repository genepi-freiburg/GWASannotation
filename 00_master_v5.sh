#TESTS and Examples

#test 1: run GWASAnnotation with a set of --eQTL_tissues_interest_coloc and then run 06_posprocessing with a different tissue selection
output_path="/data/programs/pipelines/GWASannotation/02_output/test_Feb2025/test"
#mkdir /data/programs/pipelines/GWASannotation/02_output/test_Feb2025/
input_file="/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/30_M30_GWASAnnotation/02_output/ME30/ME30_liftOver_hg38_dedup.RDS"
sbatch -p TempCompute --nodelist="imbi12" --output  /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_test_Feb2025.txt --job-name=test1 --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="/scratch/global/martins/R-4.4.1/bin/Rscript /data/programs/pipelines/GWASannotation_dev/00_scripts/GWASAnno_main.R  \
--GWAS_RDS $input_file \
--bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/02_HRC/Plink_R2_0.3_MAF_0.01/hg38/GCKD_HRC_R2_03_MAF_001_hg38_chr1-22" \
--eQTL_datasets_coloc "Kidney_eQTL,GTEXv8,eQTLGen" \
--eQTL_tissues_interest_coloc "Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" \
--coloc_input_path "/data/programs/pipelines/GWASannotation/02_output/test_Feb2025/coloc/output/" \
--output_path $output_path "

output_path="/data/programs/pipelines/GWASannotation/02_output/test_Feb2025/test"
#mkdir /data/programs/pipelines/GWASannotation/02_output/test_Feb2025/
sbatch -p TempCompute --nodelist="imbi12" --output  /data/programs/pipelines/GWASannotation/03_logs/06_postprocessing_test_Feb2025.txt --job-name=Anno --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="/scratch/global/martins/R-4.4.1/bin/Rscript /data/programs/pipelines/GWASannotation_dev/00_scripts/06_postprocessing.R  \
    --output_path $output_path \
    --output_file_name "GWASAnno_summary_v2.txt" \
    --eQTL_tissues_interest_coloc "Adrenal_Gland,Heart_Left_Ventricle,Kidney_Cortex,Liver,Pancreas,Small_Intestine_Terminal_Ileum,Colon_Transverse,Stomach,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" "

#test 2, run GWASAnnotation without --eQTL_tissues_interest_coloc and add afterwards
output_path="/data/programs/pipelines/GWASannotation/02_output/test2_Feb2025/test"
mkdir /data/programs/pipelines/GWASannotation/02_output/test2_Feb2025/
input_file="/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/30_M30_GWASAnnotation/02_output/ME30/ME30_liftOver_hg38_dedup.RDS"
sbatch -p TempCompute --nodelist="imbi12" --output  /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_test2_Feb2025.txt --job-name=Anno --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="/scratch/global/martins/R-4.4.1/bin/Rscript /data/programs/pipelines/GWASannotation_dev/00_scripts/GWASAnno_main.R  \
--GWAS_RDS $input_file \
--bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/02_HRC/Plink_R2_0.3_MAF_0.01/hg38/GCKD_HRC_R2_03_MAF_001_hg38_chr1-22" \
--eQTL_datasets_coloc "Kidney_eQTL,GTEXv8,eQTLGen" \
--coloc_input_path "/data/programs/pipelines/GWASannotation/02_output/test_Feb2025/coloc/output/" \
--output_path $output_path "

output_path="/data/programs/pipelines/GWASannotation/02_output/test2_Feb2025/test"
#mkdir /data/programs/pipelines/GWASannotation/02_output/test_Feb2025/
sbatch -p TempCompute --nodelist="imbi12" --output  /data/programs/pipelines/GWASannotation/03_logs/06_postprocessing_test2_Feb2025.txt --job-name=Anno --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="/scratch/global/martins/R-4.4.1/bin/Rscript /data/programs/pipelines/GWASannotation_dev/00_scripts/06_postprocessing.R  \
    --output_path $output_path \
    --output_file_name "GWASAnno_summary_v2.txt" \
    --eQTL_tissues_interest_coloc "Adrenal_Gland,Heart_Left_Ventricle,Kidney_Cortex,Liver,Pancreas,Small_Intestine_Terminal_Ileum,Colon_Transverse,Stomach,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" "

#Test Sahar sumstats after merge with main folder - use only 1 region on chr 4
#Summary statistics: /data/meta_analyses/00_CKDGen/01_analyses/09_GWAS-r5/urate-Sahar/urate-multi/urate_20250103.multi.tbl.gz
#LD file: /data/studies/06_UKBB/02_Projects/19_LD_RefPanel/bed_ABC/ALL/LD_Panel_ALL_combined
#Genome build: GRCh38

sbatch -p TempCompute --nodelist="imbi12" --output  /data/programs/pipelines/GWASannotation/03_logs/00_process_sumstats_metal_Sahar_test.txt --job-name=pp --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/05_internal/00_process_sumstats_metal_Sahar_test.R"



mkdir /data/programs/pipelines/GWASannotation/02_output/test_Sahar_2025/
output_path="/data/programs/pipelines/GWASannotation/02_output/test_Sahar_2025/test_chr4"
input_file="/data/programs/pipelines/GWASannotation/02_output/test_Sahar_2025/test_chr4_hg38_dedup.RDS"
sbatch -p TempCompute --nodelist="imbi12" --output  /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_test_Sahar_chr4_Feb2025.txt --job-name=Anno --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="/scratch/global/martins/R-4.4.1/bin/Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R  \
--GWAS_RDS $input_file \
--bfile "/data/studies/06_UKBB/02_Projects/19_LD_RefPanel/bed_ABC/ALL/LD_Panel_ALL_combined" \
--eQTL_datasets_coloc "GTEXv8,eQTLGen" \
--output_path $output_path "
