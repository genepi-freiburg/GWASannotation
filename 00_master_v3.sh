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

    
