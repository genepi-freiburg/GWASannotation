#devtools::load_all("/data/programs/pipelines/genepicoloc/genepicoloc_package")

# Set parameters
cat("Colocalization analysis \n")
file <- paste0(output_path, "_subset.tsv.gz")
#files <- "/data/programs/pipelines/GWASannotation/02_output/test/test_subset.txt.gz"
regions <- paste0(output_path, "_coloc_regions.RDS")

folder_path <- dirname(output_path)
folder_path <- paste0(folder_path, "/coloc")
dir.create(folder_path)
setwd(folder_path)

regions <- readRDS(regions)

#Colocalization module
args_list_selected <- args_list_wrapper(region_args=regions$coloc_regions_PASS,
                                        sumstats_1_file=file,
                                        sumstats_1_type=sumstats_type,
                                        sumstats_1_sdY=NA,
                                        selected_studies=datasets_coloc)
print("Datasets that will be used for colocalization analysis: ")
print(datasets_coloc)

coloc_out_all <- coloc_out_wrapper(args_list_selected, output_folder = "output", mc_cores=10)
# change number of cores with mc_cores
# use debug_mode=T for debugging

coloc_out_all_annot <- coloc_out_annotate(coloc_out_all, output_folder = "output")
coloc_out_summary(coloc_out_all_annot, output_folder = "output")

