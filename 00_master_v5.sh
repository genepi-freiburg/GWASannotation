mkdir /data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/02_GWASAnnotation/02_output/glycopeptides/v4_allData/
traits=("HYT_N3H3S2_qnorm")
for trait in "${traits[@]}"; do
    input_file="/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/01_regenie/02_output/step2/glycopeptides/v4_allData/model1/${trait}/maf001/${trait}_hg38_dedup.RDS"
    output_path="/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/02_GWASAnnotation/02_output/glycopeptides/v4_allData/${trait}/${trait}"
    mkdir /data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/02_GWASAnnotation/02_output/glycopeptides/v4_allData/${trait}
    log_file="/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/05_regenie/03_logs/GWASAnno_main_v3_allData_${trait}_20250410.txt"

    if [[ $(zcat /data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/01_regenie/02_output/step2/glycopeptides/v4_allData/model1/${trait}/maf001/${trait}_chr1-22_X_maf001.regenie.gz | awk 'NR>1 && $13 > 7.3' | wc -l) -gt 0 ]]; then
            echo ${trait}
            mkdir /data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/01_regenie/step2/glycopeptides/v4_allData/${trait}
            if ! grep -q "Pipeline complete!" $log_file; then
                sbatch -p imbiPCompute --nodelist="imbip-compute-214" --output /data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/02_GWASAnnotation/03_logs/GWASAnno_main_allData_${trait}_$(date +'%Y%m%d').txt --job-name=A_${trait} --begin=now --mem=40G --mail-type=FAIL --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="/scratch/global/martins/R-4.4.1/bin/Rscript /data/programs/pipelines/GWASannotation_dev/00_scripts/GWASAnno_main.R \
                --GWAS_RDS $input_file \
                --bfile "/data/studies/00_GCKD/00_data/01_genotypes/03_imputed_data/05_TopMed_2022/07_Plink/GCKD-TopMed_rsq03_maf001" \
                --eQTL_datasets_coloc "GTEXv8" \
                --pQTL_datasets_coloc "Icelanders_pGWAS" \
                --coloc_input_path "/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/32_IgA_N-glycome/02_GWASAnnotation/02_output/glycopeptides/v4_allData/${trait}/coloc/output/" \
                --output_path $output_path "
            else
                echo "Pipeline already complete for ${trait}, skipping sbatch."
            fi
        else
            echo ${trait} "HAS NO SIGNIFICANT FINDINGS"
        fi
done
