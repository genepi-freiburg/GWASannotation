# Set parameters
cat("Colocalization analysis \n")
file <- paste0(output_path, "_subset.tsv.gz")
#files <- "/data/programs/pipelines/GWASannotation/02_output/test/test_subset.txt.gz"
regions <- paste0(output_path, "_coloc_regions.RDS")
sumstats_1_type="quant"

folder_path <- dirname(output_path)
folder_path <- paste0(folder_path, "/coloc")
dir.create(folder_path)
setwd(folder_path)

# Run colocs
#cat("check R environment \n")
#print(Sys.getenv())

####################################
regions <- readRDS(regions)

sumstats_1_args <- data.frame(sumstats_1_file = file,
                              sumstats_1_function = "query_sumstats_1",
                              sumstats_1_type = sumstats_1_type,
                              sumstats_1_sdY = NA)
sumstats_1_args <- data.frame(sumstats_1_args, regions[["coloc_regions_PASS"]])

## ---- echo=T, eval = T--------------------------------------------------------
#selected_studies <- c("ARIC_pGWAS", "GTEXv8", "Icelanders_pGWAS", "Kidney_eQTL", "UKB_PPP_EUR")
selected_studies <- datasets_coloc
#cat("coloc datasets: ", selected_studies, "\n")
list_to_create_args <- list_to_create_args_list[selected_studies]
list_of_args <- lapply(list_to_create_args, function(x) {
  do.call(create_coloc_params_df,
          c(x, list(sumstats_1_args = sumstats_1_args)))
})
str(list_of_args, 1)

##############
coloc_out <- Map(parallel_wrapper, list_of_args, dry_run = F, N_cpus_per_node = 10) #change back to 10

summarize_coloc(selected_studies=selected_studies,
                output_folder = "output",
                remove_dirname = F)

#
