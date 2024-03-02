# Set parameters
file <- paste0(output_path, "_subset.txt.gz")
#files <- "/data/programs/pipelines/GWASannotation/02_output/test/test_subset.txt.gz"
regions <- paste0(output_path, "_coloc_regions.RDS")
sumstats_1_type="quant"

folder_path <- dirname(output_path)
folder_path <- paste0(folder_path, "/coloc")
dir.create(folder_path)
setwd(folder_path)
# Run colocs


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
list_to_create_args <- list_to_create_args_list[selected_studies]
list_of_args <- lapply(list_to_create_args, function(x) {
  do.call(create_coloc_params_df,
          c(x, list(sumstats_1_args = sumstats_1_args)))
})
str(list_of_args, 1)

parallel_wrapper <- function(args_df, N_cpus_per_node = 10, output_folder="output",
                             do_rbind = T, do_annotate = NULL, do_annotate_sumstats_1 = NULL,
                             save_RDS = T, save_RDS_no_annotation = F,
                             dry_run = F, debug_mode = F,
                             min_nlog10P = -log10(1e-5),
                             run_slurm = NULL) {
  # Setup
  if (!is.null(run_slurm)) {
    stop("Slurm functionality is (temporarily) disabled.",
         " Please remove 'run_slurm' parameter and the function will execute the standard run.",
         " Wrap in slurm separately, if necessary (e.g., for each sumstats_2).")
  }
  if (!"parallel" %in% rownames(installed.packages())) {
    stop("'parallel' is currently required to run this function'")
  }
  system(paste0("mkdir -p ", output_folder))
  # Arguments
  EXPERIMENT <- args_df$EXPERIMENT
  print(EXPERIMENT)
  extra_args <- args_df$extra_args
  if (!is.null(extra_args)) {
    extra_args <- list(extra_args, min_nlog10P = min_nlog10P)
  } else {
    extra_args <- list(min_nlog10P = min_nlog10P)
  }
  params_df <- args_df$params_df
  if (is.null(do_annotate)) {
    if (!is.null(args_df$annotate)) {
      do_annotate <- args_df$annotate$do_annotate
      annotation_function <- args_df$annotate$annotation_function
      annotation_function_args <- args_df$annotate$annotation_function_args
    } else {
      do_annotate <- F
    }
  }
  if (is.null(do_annotate_sumstats_1)) {
    if (!is.null(args_df$annotate_sumstats_1)) {
      do_annotate_sumstats_1 <- args_df$annotate_sumstats_1$do_annotate_sumstats_1
      annotation_function_sumstats_1 <- args_df$annotate_sumstats_1$annotation_function_sumstats_1
      annotation_function_args_sumstats_1 <- args_df$annotate_sumstats_1$annotation_function_args_sumstats_1
    } else {
      do_annotate_sumstats_1 <- F
    }
  }
  if (do_annotate | do_annotate_sumstats_1) {  if (!do_rbind) { stop("do_annotate=T or do_annotate_sumstats_1=T cannot be used with do_rbind=F")}  }
  if (dry_run & !(debug_mode)) {params_df <- params_df[1:2,]; EXPERIMENT <- paste0(EXPERIMENT, "_dryrun")}
  if (debug_mode) {
    coloc_out <- lapply(1:nrow(params_df),
                        function(i) {
                          print(i); do.call(coloc_wrapper, c(params_df[i,], extra_args))
                        })
  } else {
    print("coloc_out")
    
    print(head(params_df))
    coloc_out <- parallel::mclapply(1:nrow(params_df),
                                    function(i) {
                                      do.call(coloc_wrapper, c(params_df[i,], extra_args))
                                    }, mc.cores = N_cpus_per_node)
                                    
    print(head(coloc_out))
  }
  if (do_rbind) {
    print("do_rbind")
    print(head(coloc_out))
    coloc_out <- do.call(rbind, lapply(coloc_out, function(x) {do.call(rbind, x)}))
  }
  if (do_annotate) {
    coloc_out <- do.call(annotation_function, c(annotation_function_args, list(coloc_out = coloc_out)))
  }
  if (do_annotate_sumstats_1) {
    coloc_out <- do.call(annotation_function_sumstats_1, c(annotation_function_args_sumstats_1, list(coloc_out = coloc_out)))
  }
  if (save_RDS_no_annotation | save_RDS) {
    
    EXPERIMENT <- paste0(output_folder, "/", EXPERIMENT)
  }
  if (save_RDS_no_annotation) {
    saveRDS(coloc_out, paste0(EXPERIMENT, "_no_annotation.RDS"))
  }
  if (save_RDS) {
    print("save_RDS")
    saveRDS(coloc_out, paste0(EXPERIMENT, ".RDS"))
  }
  return(coloc_out)
}

##########################################
query_sumstats_1 <- function(sumstats_file,
                             CHR_var, BP_START_var, BP_STOP_var,
                             ...,
                             read_mode = "tabix") {
                                
  if (read_mode == "tabix") {
    sumstats <- read.table(text=system(paste0("tabix -h ", sumstats_file, " ",
                                              CHR_var, ":", BP_START_var, "-",
                                              BP_STOP_var), intern = T), sep = "\t", header = T)
    colnames(sumstats) =c("Name", "rsID", "CHR", "POS", "A1", "A2", "BETA", "SE", "nlog10P", "AF", "N")
  }
  if (read_mode == "RDS") {
    sumstats <- readRDS(sumstats_file)
  }
  if (read_mode == "read.csv") {
    sumstats <- read.csv(sumstats_file)
  }
  if (read_mode == "get") {
    sumstats <- get(sumstats_file)
  }
  if (nrow(sumstats) == 0) { return(sumstats) }
  sumstats <- subset(sumstats, CHR == CHR_var & POS >= BP_START_var & POS <= BP_STOP_var)
  return(sumstats)
}
##########################################
coloc_wrapper <- function(CHR_var, BP_START_var, BP_STOP_var,
                          sumstats_1_file, sumstats_1_function,
                          ..., sumstats_1_type, sumstats_1_sdY,
                          sumstats_2_file, sumstats_2_function,
                          sumstats_2_type, sumstats_2_sdY,
                          hyprcoloc_mode = F,
                          do_process_wrapper = T, min_nlog10P = -log10(1e-5)) {
  # Declare nested function
  process_sumstats_2_df <- function(sumstats_1_df, sumstats_2_df) {
    # if sumstats_1/2 queries return 0 rows or there is no SNP intersect
    no_intersect <- all(!(sumstats_1_df$Name %in% sumstats_2_df$Name))
    if (nrow(sumstats_1_df) == 0 | nrow(sumstats_2_df) == 0 | no_intersect) {
      if (do_process_wrapper == F) {stop("Output is not consistent when nrow=0 and do_process_wrapper=F")}
      coloc_output <- out_template(CHR_var, BP_START_var, BP_STOP_var,
                                   sumstats_1_file, sumstats_1_max_nlog10P=NA,
                                   sumstats_2_file, sumstats_2_max_nlog10P=NA,
                                   nsnps=0)
    } else {
      # calculate max nlog10 or min P
      if (is.data.frame(sumstats_1_df)) {
        sumstats_1_max_nlog10P <- max(sumstats_1_df[["nlog10P"]], na.rm=T)
      } else { sumstats_1_max_nlog10P <- NULL }
      if (is.data.frame(sumstats_2_df)) {
        sumstats_2_max_nlog10P <- max(-log10(sumstats_2_df[["P"]]), na.rm=T)
      } else { sumstats_2_max_nlog10P <- NULL }
      # Do not run coloc if there are no significant SNP
      if (sumstats_1_max_nlog10P < min_nlog10P | sumstats_2_max_nlog10P < min_nlog10P) {
        coloc_output <- out_template(CHR_var, BP_START_var, BP_STOP_var,
                                     sumstats_1_file, sumstats_1_max_nlog10P,
                                     sumstats_2_file, sumstats_2_max_nlog10P)
      } else {
        # run coloc if both sumstats have significant SNPs
        coloc_output <- run_coloc(sumstats_1_df = sumstats_1_df,
                                  sumstats_1_type = sumstats_1_type,
                                  sumstats_1_sdY = sumstats_1_sdY,
                                  sumstats_2_df = sumstats_2_df,
                                  sumstats_2_type = sumstats_2_type,
                                  sumstats_2_sdY = sumstats_2_sdY)
        coloc_output$region <- data.frame(CHR_var = CHR_var,
                                          BP_START_var = BP_START_var,
                                          BP_STOP_var = BP_STOP_var,
                                          sumstats_1_file = sumstats_1_file,
                                          sumstats_1_max_nlog10P = sumstats_1_max_nlog10P,
                                          sumstats_2_file = sumstats_2_file,
                                          sumstats_2_max_nlog10P = sumstats_2_max_nlog10P)

        if (do_process_wrapper) {
          coloc_output <- process_wrapper(coloc_output)
        }
      }
    }
    return(coloc_output)
  }
  # run code
  extra_args <- list(...)
  args_list <- list(CHR_var = CHR_var, BP_START_var = BP_START_var,
                    BP_STOP_var = BP_STOP_var)
  if (length(extra_args) > 0) {
    args_list <- c(args_list, extra_args)
  }
  sumstats_1_df <- do.call(sumstats_1_function,
                           c(list(sumstats_file = sumstats_1_file), args_list))
  #head(sumstats_1_df)
  if (hyprcoloc_mode) {
    return(sumstats_1_df)
  }
  # handle character P in case of underflow - already done at preprocessing
  # sumstats_1_df[["P"]] <- as.numeric(sumstats_1_df[["P"]])
  sumstats_2_obj <- do.call(sumstats_2_function,
                            c(list(sumstats_file = sumstats_2_file), args_list))
  #head(sumstats_2_obj)

  # sumstats_2_obj can be either a list of data.frames or a data.frame
  # next block will process sumstats_2_obj as a list, so convert first if needed
  if (is.data.frame(sumstats_2_obj)) {sumstats_2_obj <- list(sumstats_2_obj)}
  # process_sumstats_2_df
  coloc_output <- lapply(sumstats_2_obj, function(sumstats_2_df) {
    # handle character P in case of underflow - already done at preprocessing
    # sumstats_2_df[["P"]] <- as.numeric(sumstats_2_df[["P"]])
    df_out <- process_sumstats_2_df(sumstats_1_df = sumstats_1_df,
                                    sumstats_2_df = sumstats_2_df)
    if ("Phenotype" %in% colnames(sumstats_2_df)) {
      df_out[["sumstats_2_file"]] <- paste0(df_out[["sumstats_2_file"]], "_", unique(sumstats_2_df[["Phenotype"]]))
    }
    return(df_out)
  })
  return(coloc_output)
}
##############
coloc_out <- Map(parallel_wrapper, list_of_args)

summarize_coloc(selected_studies=selected_studies,
                output_folder = "output",
                remove_dirname = F)
