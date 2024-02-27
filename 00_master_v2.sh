#Prepare input for pipeline
#ouput:
##_liftOver_hg38.tx
##_liftOver_hg38.RDS
gwas="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/20Oct/maf001/model1_qnorm_tkv.bsa_chr1-22_maf001.regenie.gz"
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
mkdir /data/programs/pipelines/GWASannotation/02_output/test2
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/00_process_sumstats_test2_$(date +'%Y%m%d').txt --job-name=preprocess --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/00_process_sumstats.R \
--GWAS $gwas \
--genome_build hg37 \
--output_path $output_path"

#Preprocess files (for VEP, magma, coloc and ProGeM pipline
#output:
##"_get_coloc_regions_log.txt" #not needed
#for coloc (and ProGeM):
##"_subset.RDS")) #not needed?
##"_subset.txt.gz" and tbi
##"_coloc_regions.RDS"
#for ProGeM:
##"_sentinel.txt"
##"_proxies.txt"
#for VEP:
#"_proxies_vep.txt"
#for magma:
#"_input_magma_rsid.txt" and "_input_magma_CHRPOS.txt" - filtered on script 03 depending if the bfile has rsids or chrpos
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
input_file="${output_path}_liftOver_hg38.RDS"
#mkdir /data/programs/pipelines/GWASannotation/02_output/test
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/01_preprocessing_test2_$(date +'%Y%m%d').txt --job-name=preprocess --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/01_preprocessing.R \
--GWAS_RDS $input_file \
--output_path $output_path"

#run VEP
#output
##_proxies_vepout.txt
##_proxies_vepout.txt_summary.html
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
input_file="${output_path}_proxies_vep.txt"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/02_vep_test2_$(date +'%Y%m%d').txt --job-name=vep --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de /data/programs/pipelines/GWASannotation/00_scripts/02_vep.sh $input_file

#run magma and Pops
#input_file="$1" #input file for magma (from scripr 01) - it will be selected _rsid.txt or _CHRPOS.txt depending on the bfile
#bfile="${2:-/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_150k_RandomSubset_Cleaned/hg38/UKBB_14k_hg38_chr1-22}"  # Default bfile if not provided #hg38
#magma_annotated_genes="$3" - generates if not given
#output_path="$4" #path for output + prefix - same as previous output path

#output_
##_input_magma_CHRPOS.txt
##_input_magma_rsid.txt
##_magma_annotate.genes.annot
##_magma_annotate.log
##_magma.genes.out
##_magma.genes.raw
##_magma.log
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
input_file="${output_path}_input_magma" # dont add the ending
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/03_magma_and_pops_test2_$(date +'%Y%m%d').txt --job-name=magma_pops --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de /data/programs/pipelines/GWASannotation/00_scripts/03_magma_and_pops.sh $input_file $output_path "" ""

#make Granges file for PoPS to use on ProGeM
#output:
##_PoPS.RData
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
input_file="${output_path}.preds"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/04_Pops_make_GRanges_test2_$(date +'%Y%m%d').txt --job-name=magma_pops --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/04_Pops_make_GRanges.R \
--PoPS_results $input_file \
--output_path $output_path"

#run coloc for eQTL and pQTL
#output in coloc/output/
##ARIC_pGWAS.RDS
##ARIC_pGWAS.xlsx
##GTEXv8.RDS
##GTEXv8.xlsx
##Icelanders_pGWAS.RDS
##Icelanders_pGWAS.xlsx
##Kidney_eQTL.RDS
##Kidney_eQTL.xlsx
##UKB_PPP_EUR.RDS
##UKB_PPP_EUR.xlsx
##summary.RDS
##summary.xlsx
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/05_coloc_test2_$(date +'%Y%m%d').txt --job-name=coloc --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/05_coloc.R \
--input_path $output_path"

#Prepare coloc QTL results for ProGEM
#output in coloc/output/
##input_Progem_eQTL.txt
##input_Progem_pQTL.txt
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/06_prepare_QTL_test2_$(date +'%Y%m%d').txt --job-name=prepcoloc --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/06_prepare_QTL.R \
--input_path $output_path"

#run ProGeM
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/ProGeM_test2_$(date +'%Y%m%d').txt --job-name=ProGeM --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/ProGeM_settings.r \
        --sentinel_filepath ${output_path}_sentinel.txt \
        --proxy_filepath ${output_path}_proxies.txt \
        --VEP_filepath ${output_path}_proxies_vepout.txt \
        --PoPS_filepath ${output_path}_PoPS.RData \
        --coloc_dir /data/programs/pipelines/GWASannotation/02_output/test2/coloc/output/ \
        --output_path $output_path"

#Postprocessing ProGEM
output_path="/data/programs/pipelines/GWASannotation/02_output/test2/test2"
ProGeM_dir="/data/programs/pipelines/GWASannotation/02_output/test2/ProGeM/"
sbatch --output /data/programs/pipelines/GWASannotation/03_logs/07_postprocessing_ProGeM_test2_$(date +'%Y%m%d').txt --job-name=posprocess --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/07_postprocessing_ProGeM.R \
    --ProGeM_dir $ProGeM_dir \
    --pops.file ${output_path}_PoPS.RData \
    --sentinel.file ${output_path}_sentinel.txt"
