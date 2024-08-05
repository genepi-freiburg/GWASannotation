## 10June 2024
#Test with everything from scratch and no options
#Test a problematic conditional summary stats from hilus
F="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/model1_qnorm_hilus.bsa_final_volume/gcta_cojo/cojo_conditional_pvav5e-8_col09/snps_model1_qnorm_hilus.bsa_4_18751406_19751406_chr4:18938095:T:C.cma.cojo"
BN=`basename ${F}`
#echo $BN
name=`echo "$BN" | sed -n 's/.*chr\([^\.]*\)\.cma\.cojo/\1/p'`

output_path="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/model1_qnorm_hilus.bsa_final_volume/GWAS_Annotation_MRI_38k_LD/cond_stats/${BN%.cma.cojo}/${name}"
input_file="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/model1_qnorm_hilus.bsa_final_volume/gcta_cojo/cojo_conditional_pvav5e-8_col09/${BN%.cma.cojo}.RDS"

sbatch -p imbiPCompute --nodelist="imbip-compute-221" --output /data/studies/06_UKBB/02_Projects/14_MRI-kidney/03_logs/GWASAnno_main_${BN%.cma.cojo}_38kLD_$(date +'%Y%m%d').txt --job-name=Anno${BN%.cma.cojo} --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/GWASannotation_dev/00_scripts/GWASAnno_main.R \
--GWAS_RDS $input_file \
--bfile "/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_MRI_38k_MAF001/plink/chrpos/hg38/UKBB_MRI_38k_MAF001_hg38_allchr" \
--eQTL_datasets_coloc "Kidney_eQTL,GTEXv8,eQTLGen" \
--eQTL_tissues_interest_coloc "Kidney_Cortex,Spleen,Whole_Blood,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs,Kidney_eQTL_Meta_S686_Significant.q0.01,CXTubsigeQTLs,CXGlomsigeQTLs" \
--coloc_input_path "/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/model1_qnorm_hilus.bsa_final_volume/GWAS_Annotation_MRI_38k_LD/cond_stats/${BN%.cma.cojo}/coloc/output/" \
--output_path $output_path "

#Test Nora
#first preprocess
output_path="/data/programs/pipelines/GWASannotation/02_output/C100000961/test"
input_file="/dsk/data1/studies/04_ARIC/01_analyses/03_paired_mGWAS/01_input/serum/EA/ARIC_EA_TopMed_2023-05-10_C100000961_271_C100000961.regenie.gz"
sbatch -p imbiPCompute --nodelist="imbip-compute-221" --output /data/programs/pipelines/GWASannotation/03_logs/00_process_C100000961_tes_$(date +'%Y%m%d')_raw.txt --job-name=prepro --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/GWASannotation_dev/04_utils/00_process_sumstats_REGENIE_Name_by_position.R \
    --GWAS $input_file \
    --output_path $output_path "
    
    
output_path="/data/programs/pipelines/GWASannotation/02_output/C100000961/test"
input_file="/data/programs/pipelines/GWASannotation/02_output/C100000961/test_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_C100000961_test_$(date +'%Y%m%d').txt --job-name=GeneAnno --begin=now --wrap="Rscript /data/programs/pipelines/GWASannotation/GWASannotation_dev/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --bfile "/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_MRI_38k_MAF001/plink/chrpos/hg38/UKBB_MRI_38k_MAF001_hg38_allchr" \
    --output_path $output_path "



output_path="/data/programs/pipelines/GWASannotation/02_output/C100000961/test"
mkdir /data/programs/pipelines/GWASannotation/02_output/C100000961/
input_file="/dsk/data1/studies/04_ARIC/01_analyses/03_paired_mGWAS/05_gene_annotation/04_test/01_input/C100000961_hg38.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_C100000961_test_$(date +'%Y%m%d').txt --job-name=GeneAnno --begin=now --wrap="Rscript /data/programs/pipelines/GWASannotation/GWASannotation_dev/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --bfile "/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_MRI_38k_MAF001/plink/chrpos/hg38/UKBB_MRI_38k_MAF001_hg38_allchr" \
    --output_path $output_path "

######test tkv

v="tkv"
output_path="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/model1_qnorm_${v}.bsa_final_volume/GWAS_Annotation_MRI_38k_LD/model1_qnorm_${v}.bsa"
input_file="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/10Feb2024_final_volumes/maf001/model1_qnorm_${v}.bsa_chr1-22_maf001_liftOver_hg38_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-213" --output /data/studies/06_UKBB/02_Projects/14_MRI-kidney/03_logs/GWASAnno_MRI_38k_LD_main_model1_qnorm_${v}.bsa_$(date +'%Y%m%d').txt --job-name=Anno${v} --begin=now+2400 --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/GWASannotation_dev/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --eQTL_datasets_coloc "Kidney_eQTL,GTEXv8,eQTLGen" \
    --eQTL_tissues_interest_coloc "Kidney_Cortex,Spleen,Whole_Blood,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs,Kidney_eQTL_Meta_S686_Significant.q0.01,CXTubsigeQTLs,CXGlomsigeQTLs" \
    --output_path $output_path "
