#Test with everything from scratch and no options
output_path="/data/programs/pipelines/GWASannotation/02_output/test_raw/test"
mkdir /data/programs/pipelines/GWASannotation/02_output/test_raw
input_file="/data/programs/pipelines/GWASannotation/02_output/test2/test2_liftOver_hg38.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_$(date +'%Y%m%d')_raw.txt --job-name=Anno_raw --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --output_path $output_path "

#Test giving coloc folder and some other options
output_path="/data/programs/pipelines/GWASannotation/02_output/test_with_coloc_path/test"
mkdir /data/programs/pipelines/GWASannotation/02_output/test_with_coloc_path
input_file="/data/programs/pipelines/GWASannotation/02_output/test2/test2_liftOver_hg38.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_$(date +'%Y%m%d')_with_coloc_path.txt --job-name=Anno_coloc --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --magma_annotated_genes "${output_path}_magma_annotate.genes.annot" \
    --eQTL_datasets_coloc "Kidney_eQTL,GTEXv8" \
    --coloc_input_path "/data/programs/pipelines/GWASannotation/02_output/test2/coloc/output/" \
    --eQTL_tissues_interest_coloc "Kidney_Cortex,Spleen,Whole_Blood,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" \
    --output_path $output_path "


#TEST Uromodulin
output_path="/data/programs/pipelines/GWASannotation/02_output/Uromodulin/"
mkdir $output_path
input_file="/data/meta_analyses/19_UMOD_2023/meta_UKB50k_Regenie/GWASannotation_run/01_input/Uromodulin_liftOver_hg38.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_Uromodulin_$(date +'%Y%m%d')_raw.txt --job-name=Anno_raw --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --eQTL_datasets_coloc "Kidney_eQTL,GTEXv8" \
    --eQTL_tissues_interest_coloc "Kidney_Cortex,Spleen,Whole_Blood,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" \
    --output_path $output_path "
    
#TEST Uromodulin with coloc results from Oleg
output_path="/data/programs/pipelines/GWASannotation/02_output/Uromodulin"
mkdir $output_path
input_file="/data/meta_analyses/19_UMOD_2023/meta_UKB50k_Regenie/GWASannotation_run/01_input/Uromodulin_liftOver_hg38.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_Uromodulin_$(date +'%Y%m%d')_with_coloc.txt --job-name=Anno_raw --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --eQTL_datasets_coloc "Kidney_eQTL,GTEXv8" \
    --coloc_input_path "/data/meta_analyses/19_UMOD_2023/meta_UKB50k_Regenie/colocalization" \
    --eQTL_tissues_interest_coloc "Kidney_Cortex,Spleen,Whole_Blood,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" \
    --output_path $output_path "

output_path="/data/programs/pipelines/GWASannotation/02_output/Uromodulin/Uromodulin"
file_out="GWAS_anno_scoring_Uromodulin_test.txt"
sbatch --output /data/programs/pipelines/GWASannotation/03_logs/08_postprocessing_Uromodulin_$(date +'%Y%m%d').txt --job-name=score --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/08_postprocessing_ProGeM.R \
    --output_path $output_path \
    --eQTL_tissues_interest_coloc "Kidney_Cortex,Spleen,Whole_Blood,Kidney_eQTL.TubsigeQTLs,Kidney_eQTL.GlomsigeQTLs" \
    --nearest 1 \
    --second_nearest  0.8  \
    --third_nearest 0.6 \
    --output_file_name $file_out "

#TEST Olink
# 3 proteins to test
#HDGF - OID21455 - cis is not the nearest
#ANGPTL3 - OID20407 - cis is not the nearest and wasn't found with Yong GWAS annotation
#FCGF2A - OID20391 - cis is the nearest but wasn't found with Yong GWAS annotation
#"OID31493" (PSCA gene)
#"OID20678" (IFNGR1)
#HDGF - OID21455 - cis is not the nearest
output_path="/data/programs/pipelines/GWASannotation/02_output/HDGF_OID21455/HDGF_OID21455"
mkdir "/data/programs/pipelines/GWASannotation/02_output/HDGF_OID21455/"
input_file="/dsk/ge_netssd/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/22_Olink/19_meta_analysis/02_output/oids/OID21455/output/Olink_2batches_METAL1.tbl.gwas.gz"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/00_process_sumstats_metal_HDGF_$(date +'%Y%m%d')_raw.txt --job-name=prepro --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/04_utils/00_process_sumstats_metal.R \
    --GWAS $input_file \
    --output_path $output_path "

output_path="/data/programs/pipelines/GWASannotation/02_output/HDGF_OID21455/HDGF_OID21455"
input_file="/data/programs/pipelines/GWASannotation/02_output/HDGF_OID21455/HDGF_OID21455_hg38_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_HDGF_$(date +'%Y%m%d')_raw.txt --job-name=HDGF --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/05_TopMed_2022/07_Plink/GCKD-TopMed_rsq03_maf001" \
    --eQTL_datasets_coloc "GTEXv8,eQTLGen" \
    --eQTL_tissues_interest_coloc "Whole_Blood" \
    --output_path $output_path "

 #run with coloc finished
output_path="/data/programs/pipelines/GWASannotation/02_output/HDGF_OID21455/HDGF_OID21455"
input_file="/data/programs/pipelines/GWASannotation/02_output/HDGF_OID21455/HDGF_OID21455_hg38_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_HDGF_$(date +'%Y%m%d')_with_coloc.txt --job-name=HDGF --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main_test.R \
    --GWAS_RDS $input_file \
    --bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/05_TopMed_2022/07_Plink/GCKD-TopMed_rsq03_maf001" \
    --eQTL_datasets_coloc "GTEXv8,eQTLGen" \
    --eQTL_tissues_interest_coloc "Whole_Blood" \
    --coloc_input_path "/data/programs/pipelines/GWASannotation/02_output/HDGF_OID21455/coloc/output/" \
    --output_path $output_path "

#test only script08
output_path="/data/programs/pipelines/GWASannotation/02_output/HDGF_OID21455/HDGF_OID21455"
file_out="GWAS_anno_scoring_test.txt"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/08_postprocessing_HDGF_$(date +'%Y%m%d').txt --job-name=HDGF --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/08_postprocessing_ProGeM.R \
    --output_path $output_path \
    --output_file_name $file_out "
    
#ANGPTL3 - OID20407 - cis is not the nearest and wasn't found with Yong GWAS annotation
output_path="/data/programs/pipelines/GWASannotation/02_output/ANGPTL3_OID20407/ANGPTL3_OID20407"
mkdir "/data/programs/pipelines/GWASannotation/02_output/ANGPTL3_OID20407/"
input_file="/dsk/ge_netssd/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/22_Olink/19_meta_analysis/02_output/oids/OID20407/output/Olink_2batches_METAL1.tbl.gwas.gz"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/00_process_sumstats_metal_ANGPTL3_$(date +'%Y%m%d').txt --job-name=prepro --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/04_utils/00_process_sumstats_metal.R \
    --GWAS $input_file \
    --output_path $output_path "

output_path="/data/programs/pipelines/GWASannotation/02_output/ANGPTL3_OID20407/ANGPTL3_OID20407"
input_file="/data/programs/pipelines/GWASannotation/02_output/ANGPTL3_OID20407/ANGPTL3_OID20407_hg38_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_ANGPTL3_$(date +'%Y%m%d')_raw.txt --job-name=ANGPTL3 --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/05_TopMed_2022/07_Plink/GCKD-TopMed_rsq03_maf001" \
    --eQTL_datasets_coloc "GTEXv8,eQTLGen" \
    --eQTL_tissues_interest_coloc "Whole_Blood" \
    --output_path $output_path "
    
#FCGF2A - OID20391 - cis is the nearest but wasn't found with Yong GWAS annotation
output_path="/data/programs/pipelines/GWASannotation/02_output/FCGF2A_OID20391/FCGF2A_OID20391"
mkdir "/data/programs/pipelines/GWASannotation/02_output/FCGF2A_OID20391/"
input_file="/dsk/ge_netssd/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/22_Olink/19_meta_analysis/02_output/oids/OID20391/output/Olink_2batches_METAL1.tbl.gwas.gz"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/00_process_sumstats_metal_FCGF2A_$(date +'%Y%m%d').txt --job-name=prepro --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/04_utils/00_process_sumstats_metal.R \
    --GWAS $input_file \
    --output_path $output_path "

output_path="/data/programs/pipelines/GWASannotation/02_output/FCGF2A_OID20391/FCGF2A_OID20391"
input_file="/data/programs/pipelines/GWASannotation/02_output/FCGF2A_OID20391/FCGF2A_OID20391_hg38_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_FCGF2A_$(date +'%Y%m%d')_raw.txt --job-name=FCGF2A --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/05_TopMed_2022/07_Plink/GCKD-TopMed_rsq03_maf001" \
    --eQTL_datasets_coloc "GTEXv8,eQTLGen" \
    --eQTL_tissues_interest_coloc "Whole_Blood" \
    --output_path $output_path "

#PSCA_OID31493
output_path="/data/programs/pipelines/GWASannotation/02_output/PSCA_OID31493/PSCA_OID31493"
mkdir "/data/programs/pipelines/GWASannotation/02_output/PSCA_OID31493/"
input_file="/dsk/ge_netssd/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/22_Olink/19_meta_analysis/02_output/oids/OID31493/output/Olink_2batches_METAL1.tbl.gwas.gz"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/00_process_sumstats_metal_PSCA_$(date +'%Y%m%d').txt --job-name=prepro --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/04_utils/00_process_sumstats_metal.R \
    --GWAS $input_file \
    --output_path $output_path "

output_path="/data/programs/pipelines/GWASannotation/02_output/PSCA_OID31493/PSCA_OID31493"
input_file="/data/programs/pipelines/GWASannotation/02_output/PSCA_OID31493/PSCA_OID31493_hg38_dedup.RDS"
sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/programs/pipelines/GWASannotation/03_logs/GWASAnno_main_PSCA_$(date +'%Y%m%d')_raw.txt --job-name=PSCA --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/GWASAnno_main.R \
    --GWAS_RDS $input_file \
    --bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/05_TopMed_2022/07_Plink/GCKD-TopMed_rsq03_maf001" \
    --eQTL_datasets_coloc "GTEXv8,eQTLGen" \
    --eQTL_tissues_interest_coloc "Whole_Blood" \
    --output_path $output_path "

#check RAP
sbatch --output /data/programs/pipelines/GWASannotation/03_logs/07_regplot_locuszoom_HDGF_$(date +'%Y%m%d').txt --job-name=locuszoom --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/studies/06_UKBB/02_Projects/14_MRI-kidney/00_scripts/07_regplot_locuszoom_test_GWASAnno_pgwas.R "


