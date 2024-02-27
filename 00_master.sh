gwas="/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/20Oct/maf001/model1_qnorm_tkv.bsa_chr1-22_maf001.regenie.gz"
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
#mkdir /data/programs/pipelines/GWASannotation/02_output/test
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/00_process_sumstats_$(date +'%Y%m%d').txt --job-name=preprocess --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/00_process_sumstats.R \
--GWAS $gwas \
--genome_build hg37 \
--output_path $output_path"

#Preprocess files (for VEP, magma, coloc and ProGeM pipline
#output:
#for coloc:
##"_get_coloc_regions_log.txt"
##"_subset.RDS")) #not needed?
##"_subset.txt.gz" and tbi
##"_coloc_regions.RDS"
#for ProGeM:
##"_sentinel.txt"
##"_proxies_vep.txt"
#for VEP:
#"_proxies_vep.txt"
#for magma:
#"_input_magma_rsid.txt" and "_input_magma_CHRPOS.txt" - filtered on script 03 depending if the bfile has rsids or chrpos

output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
input_file="${output_path}_liftOver_hg38.RDS"
#mkdir /data/programs/pipelines/GWASannotation/02_output/test
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/01_preprocessing_$(date +'%Y%m%d').txt --job-name=preprocess --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/01_preprocessing.R \
--GWAS_RDS $input_file \
--output_path $output_path"

#run VEP
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
input_file="${output_path}_proxies_vep.txt"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/02_vep_$(date +'%Y%m%d').txt --job-name=vep --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de /data/programs/pipelines/GWASannotation/00_scripts/02_vep.sh $input_file

#run magma and Pops
#input_file="$1" #input file for magma (from scripr 01) - it will be selected _rsid.txt or _CHRPOS.txt depending on the bfile
#bfile="${2:-/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_150k_RandomSubset_Cleaned/hg38/UKBB_14k_hg38_chr1-22}"  # Default bfile if not provided #hg38
#magma_annotated_genes="$3" - generates if not given
#output_path="$4" #path for output + prefix - same as previous output path
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
input_file="${output_path}_input_magma" # dont add the ending
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/03_magma_and_pops_$(date +'%Y%m%d').txt --job-name=magma_pops --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de /data/programs/pipelines/GWASannotation/00_scripts/03_magma_and_pops.sh $input_file $output_path "" ""

#make Granges file for PoPS to use on ProGeM
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
input_file="${output_path}.preds"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/04_Pops_make_GRanges_$(date +'%Y%m%d').txt --job-name=magma_pops --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/04_Pops_make_GRanges.R \
--PoPS_results $input_file \
--output_path $output_path"

#run coloc for eQTL and pQTL
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/05_coloc_$(date +'%Y%m%d').txt --job-name=coloc --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/05_coloc.R \
--input_path $output_path"

#Prepare coloc QTL results for ProGEM
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/06_prepare_QTL_$(date +'%Y%m%d').txt --job-name=prepcoloc --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/06_prepare_QTL.R \
--input_path $output_path"

#run ProGeM
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
sbatch  --output /data/programs/pipelines/GWASannotation/03_logs/ProGeM_$(date +'%Y%m%d').txt --job-name=ProGeM --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/ProGeM_settings.r \
        --sentinel_filepath ${output_path}_sentinel.txt \
        --proxy_filepath ${output_path}_proxies.txt \
        --VEP_filepath ${output_path}_proxies_vepout.txt \
        --PoPS_filepath ${output_path}_PoPS.RData \
        --coloc_dir /data/programs/pipelines/GWASannotation/02_output/test/coloc/output/ \
        --output_path $output_path"

#Postprocessing ProGEM
output_path="/data/programs/pipelines/GWASannotation/02_output/test/test"
ProGeM_dir="/data/programs/pipelines/GWASannotation/02_output/test/ProGeM/"
sbatch --output /data/programs/pipelines/GWASannotation/03_logs/07_postprocessing_ProGeM_$(date +'%Y%m%d').txt --job-name=posprocess --begin=now --mail-type=END --mail-user=sara.monteiro.martins@uniklinik-freiburg.de --wrap="Rscript /data/programs/pipelines/GWASannotation/00_scripts/07_postprocessing_ProGeM.R \
    --ProGeM_dir $ProGeM_dir \
    --pops.file ${output_path}_PoPS.RData \
    --sentinel.file ${output_path}_sentinel.txt"
