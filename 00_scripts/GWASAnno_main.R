### IMPORT RELEVANT PACKAGES
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(parallel)
library(coloc)
library(data.table)
library(readxl)
library(writexl)
library(dplyr)
library(biomaRt)

devtools::load_all("/data/programs/pipelines/genepicoloc/genepicoloc_package")
sapply(list.files("/data/programs/pipelines/genepicoloc/custom_scripts/source", full.names = T), source)

suppressMessages(library(optparse))
cat("\nImported required packages.\n")
### Functions
source(file = file.path("/data/programs/pipelines/GWASannotation/04_utils/ProGeM_functions.R"))
start_time_main <- Sys.time()

#create a temporary symlink with my python path
system("TEMP_SYMLINK=\"$HOME/bin/python\"; if [ ! -f $TEMP_SYMLINK ]; then ln -s /scratch/global/martins/anaconda3/bin/python3.9 $TEMP_SYMLINK; fi; which python")
#system("ln -s /scratch/global/martins/anaconda3/bin/python3.9 $TEMP_SYMLINK") #use my python path
## check python version and dependencies
system_status <- system("/data/programs/pipelines/GWASannotation/04_utils/check_environment.sh")

# Check the return status
if (system_status != 0) {
  cat("Environment check failed\n")
  quit(save = "no", status = 1)
}


option_list = list(
    make_option("--GWAS_RDS", action="store", default=NA, type='character', help="GWAS summary stats .RSD [required]"),
    make_option("--output_path", action="store", default=NA, type='character', help="output folder path [required]"),
    make_option("--bfile", action="store", default="/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_150k_RandomSubset_Cleaned/hg38/UKBB_14k_hg38_chr1-22", type='character', help="bfile to use for selecting proxies and for MAGMA (if not provided UKB_14K_hg38 will be used)"),
    make_option("--magma_annotated_genes", action="store", default="", type='character', help="magma_annotated_genes to use for MAGMA (if not provided it will be created)"),
    make_option("--eQTL_datasets_coloc", action="store", default="GTEXv8", type='character', help="comma separated eQTL datasets to use for coloc (default=GTEXv8"),
    make_option("--eQTL_tissues_interest_coloc", action="store", default=NA, type='character', help="comma separated tissues of interest to be selected from the eQTL datasets -check 04_utils/tissues_eQTL.txt"),
    make_option("--pQTL_datasets_coloc", action="store", default=c("Icelanders_pGWAS", "UKB_PPP_EUR"), type='character', help="comma separated pQTL datasets to use for coloc (default=c('Icelanders_pGWAS','UKB_PPP_EUR')"),
    make_option("--coloc_input_path", action="store", default=NA, type='character', help="path with coloc results for the selected databases (default=NA and runs coloc analysis)"),
    make_option("--r2_thresh", action="store", default=0.8, type='numeric', help="Threshold for r2 values (default=0.8)"),
    make_option("--interval_window_kb", action="store", default=500, type='numeric', help="Genomic interval (in kb) either side of the sentinel SNP (default=500), used to define genes that nearby the index SNP"),
    make_option("--LD_region_overhang_kb", action="store", default=5, type='numeric', help="Genomic interval (in kb) of the overhang of the left-most and right-most proxy/index SNP at each locus (default=5), used to identify genes with variants in LD with index or proxy SNPs"),
    make_option("--sumstats_2_max_nlog10P_thresh", action="store", default=-log10(5e-8), type='numeric', help="sumstats_2_max_nlog10P will be filtered for > this threshold for the eQTL and pQTL datasets used for colocalization (default=-log10(5e-8))"),
    make_option("--PP.H4.abf_thresh", action="store", default=0.8, type='numeric', help="PP.H4.abf (coloc) will be filtered for > this threshold (default=0.8)"))

 

opt = parse_args(OptionParser(option_list=option_list))

#############################################################################
## 1. DIRECTORIES AND FILES
folder_path <- dirname(opt$output_path)
output_dir <- paste0(folder_path, "/GWASAnno/")
dir.create(output_dir)

# Directory for GTEx eQTL data set:
eQTLdata_dir <- "/data/public_resources/GTeX/eQTLs/V8/GTEx_Analysis_v8_eQTL"

# File containing reference genes:
gene_model_filename <- "/data/programs/bin/gwas/PRoGeM/GRCh38_genes.RData"                            # default provided.

#MAGMA files
bfile=opt$bfile
magma_annotated_genes=opt$magma_annotated_genes

#Coloc datasets


possible_eQTL_dataset <- c("GTEXv8", "Kidney_eQTL", "eQTLGen")
eQTL_datasets_coloc <- unlist(strsplit(opt$eQTL_datasets_coloc, ","))
if (!all(eQTL_datasets_coloc %in% possible_eQTL_dataset)) {
  cat("--eQTL_datasets_coloc (", eQTL_datasets_coloc,  ") doesn't correspond to an available eQTL dataset \n")
  cat("Available eQTL datasets are:", paste(possible_eQTL_dataset,collapse = ", "),  "\n")
  stop("Aborting pipeline due to invalid eQTL dataset.")
}
possible_pQTL_dataset <- c("ARIC_pGWAS", "Icelanders_pGWAS","UKB_PPP_EUR")
pQTL_datasets_coloc<- unlist(strsplit(opt$pQTL_datasets_coloc, ","))
if (!all(pQTL_datasets_coloc %in% possible_pQTL_dataset)) {
  cat("--pQTL_datasets_coloc (", pQTL_datasets_coloc,  ") doesn't correspond to an available pQTL dataset \n")
  cat("Available pQTL datasets are:", paste(possible_pQTL_dataset,collapse = ", "),  "\n")
  stop("Aborting pipeline due to invalid pQTL dataset.")
}

datasets_coloc=c(eQTL_datasets_coloc, pQTL_datasets_coloc)

#tissues of interest - needs a couple of CHECKS!!!!
tissues_interest <- unlist(strsplit(opt$eQTL_tissues_interest_coloc, ","))


if(is.na(opt$coloc_input_path)){
    coloc_path <- paste0(folder_path, "/coloc/output/")   ################ ADD CHECK IF coloc_input_path has the datasets selected
} else {
    files <- list.files(opt$coloc_input_path, pattern = "\\.RDS$")
    files <- files[!grepl("summary\\.RDS$", files)]

    if (all(paste0(datasets_coloc, ".RDS") %in% files)) {
      cat("Coloc folder provided contains all selected QTL datasets \n")
    } else {
      cat("Coloc folder provided does not contain the selected QTL datasets \n")
      stop("Aborting pipeline")
    }
}


#------------------------------------------------------------------------------------------------------
## 2. PARAMETERS FOR TOP-DOWN APPROACH

# Set filter for proxy variants based on r2 values:
# filtering_required <- TRUE
filtering_required <- FALSE
		# TRUE or FALSE.

# If TRUE, set threshold for r2 values:
#r2_thresh <- 0.8										# default is 0.8.
r2_thresh <- opt$r2_thresh

# Genomic interval (in kb) either side of the sentinel SNP:
#interval_kb <- 1000										# default is 1000kb.
interval_kb <- opt$interval_window_kb
#------------------------------------------------------------------------------------------------------
## 3. PARAMETERS FOR BOTTOM-UP APPROACH

# Genomic interval (in kb) of the overhang of the left-most and right-most proxy/sentinel variant at 
# each locus:
LD_region_overhang_kb <- opt$LD_region_overhang_kb	# default is 5kb.

# Number of nearest genes that reside nearest to the sentinel variant:
number_of_nearest <- 3		# default is 3.

# Biotype(s) of candidate genes:
# further information: http://vega.archive.ensembl.org/info/about/gene_and_transcript_types.html
biotype_of_interest <- "protein_coding"								# default is "protein_coding".

# Column indices in the VEP output that contain the following information:
sentinel_rsID_column <- 13
proxy_rsID_column <- 1
ensembl_gene_id_column <- 4
IMPACT_column <- 14
r2_column <- NULL    # OPTIONAL: only required if your VEP output needs to be filtered based on r2_thresh

# Tissue(s) of interest for GTEx eQTL data set:
#GTEx_tissues <- dir(eQTLdata_dir)[grep("signifpairs", dir(eQTLdata_dir))]     # variable contains filenames
GTEx_tissues <- dir(eQTLdata_dir)[grep("signif_variant_gene_pairs", dir(eQTLdata_dir))]

tissues_of_interest <- GTEx_tissues[2]     # default is all tissues.
# tissues_of_interest <- GTEx_tissues[48]
                                          # alternatively the user can select specific tissues by
                                          # providing the appropriate indices in the square brackets.

# Column indices and threshold in the COLOC eQTL file
coloc_eqtl_sentinel_rsID_col <- 1
# Index of Ensembl gene ID column
coloc_eqtl_ensembl_gene_id_col <- 7
# Index of the tissue column
eqtl_tissue_col <- 6

# Only these tissues will be considered
# Reformat the tissues_of_interest from GTEx eQTL lookup above and include
tissues_of_interest.short <-sub(pattern=".v8.signif_variant_gene_pairs.txt.gz",
                                replacement="", x=GTEx_tissues, fixed=TRUE)
eqtl_tissues_of_interest <- c(tissues_of_interest.short,
                              "Kidney_eQTL.GlomsigeQTLs", "Kidney_eQTL_Meta_S686_Significant.q0.01",
                              "Kidney_eQTL.TubsigeQTLs")

eqtl_sumstats_2_max_nlog10P_col <- 4
# If set, sumstats_2_max_nlog10P will be filtered for > this threshold
eqtl_sumstats_2_max_nlog10P_thresh <- opt$sumstats_2_max_nlog10P_thresh
eqtl_PP.H4.abf_col <- 5
# If set, PP.H4.abf will be filtered for > this threshold
eqtl_PP.H4.abf_thresh <- opt$PP.H4.abf_thresh
# Index of the gene type column
coloc_eqtl_gene_type_col <- 8
eqtl_cis_trans_col <- 10
eqtl_cis_trans_sel <- c("cis")

# Column indices and threshold in the COLOC pQTL file
coloc_pqtl_sentinel_rsID_col <- 1
coloc_pqtl_ensembl_gene_id_col <- 6
pqtl_sumstats_2_max_nlog10P_col <- 4
pqtl_sumstats_2_max_nlog10P_thresh <- opt$sumstats_2_max_nlog10P_thresh
pqtl_PP.H4.abf_col <- 5
pqtl_PP.H4.abf_thresh <- opt$PP.H4.abf_thresh
pqtl_cis_trans_col <- 7
pqtl_cis_trans_sel <- c("cis")

#------------------------------------------------------------------------------------------------------
## 4. EXECUTE ALL

cat("##################################################\n######## STEP1: Performing pre-processing ########\n ################################################## \n")
start_time1 <- Sys.time()
output_path = opt$output_path
GWAS_RDS = opt$GWAS_RDS
source("/data/programs/pipelines/GWASannotation/00_scripts/01_preprocessing.R")
end_time1 <- Sys.time()

execution_time_seconds <- as.numeric(difftime(end_time1, start_time1, units = "secs"))
# Convert seconds to minutes
execution_time_minutes <- execution_time_seconds / 60
# Print the execution time
cat(sprintf("STEP1 done; took %.2f minutes\n", execution_time_minutes), "\n\n")
rm(execution_time_seconds)
rm(execution_time_minutes)

if (file.exists(paste0(output_path,"_subset.tsv.gz")) &
  file.exists(paste0(output_path,"_coloc_regions.RDS")) &
  file.exists(paste0(output_path,"_sentinel.txt")) &
  file.exists(paste0(output_path,"_proxies.txt")) &
  file.exists(paste0(output_path,"_proxies_vep.txt")) &
  file.exists(paste0(output_path,"_input_magma.txt"))) {
  cat("##################################################\n################# STEP2: run VEP ################# \n################################################## \n")
  input_vep=paste0(output_path,"_proxies_vep.txt")
  start_time2 <- Sys.time()
  system(paste0("/data/programs/pipelines/GWASannotation/00_scripts/02_vep.sh ",input_vep))
  end_time2 <- Sys.time()

  execution_time_seconds <- as.numeric(difftime(end_time2, start_time2, units = "secs"))
  # Convert seconds to minutes
  execution_time_minutes <- execution_time_seconds / 60
  # Print the execution time
  cat(sprintf("STEP2 done; took %.2f minutes\n", execution_time_minutes), "\n\n")
  rm(execution_time_seconds)
  rm(execution_time_minutes)
} else {
  # If output file doesn't exist, print an error message and abort the pipeline
  stop("Output file of 01_preprocessing.R does not exist. Aborting the pipeline.")
}
if (is.na(opt$coloc_input_path)){
    cat("##################################################\n########### STEP3: run MAGMA and PoPS ############\n########### AND ############\n\n########### ###### STEP5: Performing eQTL and pQTL colocalization ######\n################################################## \n")
    input_magma=paste0(output_path,"_input_magma.txt")
    start_time3 <- Sys.time()
    system(paste0("/data/programs/pipelines/GWASannotation/00_scripts/03_magma_and_pops.sh ", input_magma, " ", output_path, " ", bfile, " ", magma_annotated_genes, " > ", output_path, "_STEP3_log.txt  2>&1"), wait = FALSE)

    source("/data/programs/pipelines/GWASannotation/00_scripts/05_coloc.R")
    
    # Wait for the process started by 03_magma_and_pops.sh to finish
    while(TRUE) {
      check_process <- system("ps aux | grep '/data/programs/pipelines/GWASannotation/00_scripts/03_magma_and_pops.sh' | grep -v grep | wc -l", intern = TRUE)
      #f=file.exists(paste0(output_path,".preds"))
      #if (f == TRUE) {
      if (as.integer(check_process) == 0) {
        break  # Exit the loop if the process has finished
      }
      
      Sys.sleep(10)  # Wait for 10 seconds before checking again
    }
    cat("STEP3 has finished.")
    cat("STEP3 log file saved in ", output_path, "_STEP3_log.txt \n\n")

    end_time3 <- Sys.time()
    execution_time_seconds <- as.numeric(difftime(end_time3, start_time3, units = "secs"))
    execution_time_minutes <- execution_time_seconds / 60
    cat(sprintf("STEP3 and STEP5 done; took %.2f minutes\n", execution_time_minutes), "\n\n")
    rm(execution_time_seconds)
    rm(execution_time_minutes)
} else {
    cat("##################################################\n########### STEP3: run MAGMA and PoPS ############\n################################################## \n")
    input_magma=paste0(output_path,"_input_magma.txt")
    start_time3 <- Sys.time()
    system(paste("/data/programs/pipelines/GWASannotation/00_scripts/03_magma_and_pops.sh ", input_magma, output_path, bfile, magma_annotated_genes, sep=" "))
    end_time3 <- Sys.time()
    execution_time_seconds <- as.numeric(difftime(end_time3, start_time3, units = "secs"))
    execution_time_minutes <- execution_time_seconds / 60
    cat(sprintf("STEP3 done; took %.2f minutes\n", execution_time_minutes), "\n\n")
    rm(execution_time_seconds)
    rm(execution_time_minutes)
    coloc_path <- opt$coloc_input_path
}
if (file.exists(paste0(output_path,".preds"))) {
  cat("##################################################\n############ STEP4: Pops_make_GRanges ############ \n################################################## \n")
  PoPS_results = paste0(output_path, ".preds")
  start_time4 <- Sys.time()
  source("/data/programs/pipelines/GWASannotation/00_scripts/04_Pops_make_GRanges.R")
  end_time4 <- Sys.time()
  execution_time_seconds <- as.numeric(difftime(end_time4, start_time4, units = "secs"))
  execution_time_minutes <- execution_time_seconds / 60
  cat(sprintf("STEP4 done; took %.2f minutes\n", execution_time_minutes), "\n\n")
  rm(execution_time_seconds)
  rm(execution_time_minutes)
  
} else {
  stop("Output file of 03_magma_and_pops.sh does not exist. Aborting the pipeline.")
}

if (file.exists(paste0(coloc_path,"/",datasets_coloc[1], ".RDS"))) {
    cat("########################################################\n###### STEP6: Prepare coloc QTL results for ProGEM ###### \n######################################################## \n")
    start_time6 <- Sys.time()
    source("/data/programs/pipelines/GWASannotation/00_scripts/06_prepare_QTL.R")
    cat("Step6 done\n\n")
    end_time6 <- Sys.time()
    execution_time_seconds <- as.numeric(difftime(end_time6, start_time6, units = "secs"))
    # Convert seconds to minutes
    execution_time_minutes <- execution_time_seconds / 60
    # Print the execution time
    cat(sprintf("STEP6 took %.2f minutes\n", execution_time_minutes), "\n\n")
    rm(execution_time_seconds)
    rm(execution_time_minutes)
    } else {
    stop("Output file of 05_coloc.R does not exist. Aborting the pipeline.")
}

if (file.exists(paste0(coloc_path,"input_Progem_eQTL.txt")) &
file.exists(paste0(coloc_path,"input_Progem_pQTL.txt"))) {
    #Files needed to run ProGeM
    #Coloc
    COLOC_EQTL_filename=paste0(coloc_path,"input_Progem_eQTL.txt")
    COLOC_PQTL_filename=paste0(coloc_path,"input_Progem_pQTL.txt")

    # Files containing GWAS sentinel and proxy variants:
    sentinel_filename <- paste0(output_path,"_sentinel.txt")
    proxy_filename <- paste0(output_path,"_proxies.txt")

    # File containing VEP annotation:
    VEP_filename <- paste0(output_path,"_proxies_vepout.txt")

    # File containing a PoPS results:#full path
    filename_PoPS <- paste0(output_path,"_PoPS.RData")

    cat("##################################################\n############### STEP7: Run ProGEM ################ \n################################################## \n")
    start_time7 <- Sys.time()
    source("/data/programs/pipelines/GWASannotation/00_scripts/07_run_ProGeM.R")
    end_time7 <- Sys.time()
    execution_time_seconds <- as.numeric(difftime(end_time7, start_time7, units = "secs"))
    execution_time_minutes <- execution_time_seconds / 60
    cat(sprintf("STEP7 took %.2f minutes\n", execution_time_minutes), "\n\n")
    rm(execution_time_seconds)
    rm(execution_time_minutes)
    } else {
    stop("Output file of 06_prepare_QTL.R does not exist. Aborting the pipeline.")
}

if (file.exists(paste0(output_dir,"OUTPUT_bottom_up_summary.txt")) &
file.exists(paste0(output_dir,"OUTPUT_top_down_scores.txt"))) {
  cat("##################################################\n###### STEP8: Pos-processing ProGEM results (include scoring) ###### \n################################################## \n")
  start_time8 <- Sys.time()
  rm(opt)
  
  source("/data/programs/pipelines/GWASannotation/00_scripts/08_postprocessing_SMM.R")
  
  end_time8 <- Sys.time()
  execution_time_seconds <- as.numeric(difftime(end_time8, start_time8, units = "secs"))
  # Convert seconds to minutes
  execution_time_minutes <- execution_time_seconds / 60
  # Print the execution time
  cat(sprintf("STEP8 took %.2f minutes\n", execution_time_minutes), "\n\n")
  rm(execution_time_seconds)
  rm(execution_time_minutes)
  } else {
  stop("Output file of 07_run_ProGeM.R does not exist. Aborting the pipeline.")
}

end_time_main <- Sys.time()
execution_time_seconds <- as.numeric(difftime(end_time_main, start_time_main, units = "secs"))
# Convert seconds to minutes
execution_time_minutes <- execution_time_seconds / 60
# Print the execution time
cat(sprintf("Pipeline took %.2f minutes\n", execution_time_minutes), "\n\n")
rm(execution_time_seconds)
rm(execution_time_minutes)
system("rm $TEMP_SYMLINK")
cat("##################################################\n############### Pipeline complete! ############### \n################################################## \n")
