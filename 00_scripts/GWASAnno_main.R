#update July 2025 - add chr x

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

devtools::load_all("/data/programs/pipelines/genepicoloc/source_package/genepicolocGWASAnnotation")

suppressMessages(library(optparse))
cat("\nImported required packages.\n")

### Functions
#Get scripts location
# Extract the script path from commandArgs
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_path <- normalizePath(script_path)
print(paste("The script path is:", script_path))
script_path <- dirname(script_path)

source(file = file.path(paste0(dirname(script_path),"/04_utils/Anno_functions.R")))
start_time_main <- Sys.time()

option_list = list(
    make_option("--GWAS_RDS", action="store", default=NA, type='character', help="GWAS summary stats .RSD [required]"),
    make_option("--output_path", action="store", default=NA, type='character', help="output folder path [required]"),
    make_option("--sumstats_type", action="store", default="quant", type='character', help="GWAS summary type (default=quant)"),
    make_option("--GWAS_max_nlog10P_thresh", action="store", default=-log10(5e-8), type='numeric', help="GWAS summary statistics will be filtered for > this threshold  (default=-log10(5e-8))"),
    make_option("--bfile", action="store", default="/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_150k_RandomSubset_Cleaned/hg38/UKBB_14k_hg38_chr1-22", type='character', help="bfile to use for selecting proxies (if not provided UKB_14K_hg38 will be used)"),
    make_option("--r2_thresh", action="store", default=0.8, type='numeric', help="Threshold for r2 values - used for proxies identification (default=0.8)"),
    make_option("--eQTL_datasets_coloc", action="store", default="GTEXv8", type='character', help="comma separated eQTL datasets to use for coloc (default=GTEXv8"),
    make_option("--eQTL_tissues_interest_coloc", action="store", default=NA, type='character', help="comma separated tissues of interest to be selected from the eQTL datasets -check 04_utils/tissues_eQTL.txt"),
    make_option("--pQTL_datasets_coloc", action="store", default=c("Icelanders_pGWAS", "UKB_PPP_EUR"), type='character', help="comma separated pQTL datasets to use for coloc (default=c('Icelanders_pGWAS','UKB_PPP_EUR')"),
    make_option("--coloc_input_path", action="store", default=NA, type='character', help="path with coloc results for the selected databases (default=NA and runs coloc analysis)"),
    make_option("--QTL_coloc_max_nlog10P_thresh", action="store", default=-log10(5e-8), type='numeric', help="eQTL and pQTL datasets used for colocalization will be filtered for > this threshold  (default=-log10(5e-8))"),
    make_option("--PP.H4.abf_thresh", action="store", default=0.8, type='numeric', help="PP.H4.abf (coloc) will be filtered for > this threshold (default=0.8)"),
    make_option("--interval_window_kb", action="store", default=500, type='numeric', help="Genomic interval (in kb) either side of the lead SNP (default=500), used to define genes that are nearby the lead SNP"),
    make_option("--LD_region_overhang_kb", action="store", default=5, type='numeric', help="Genomic interval (in kb) of the overhang of the left-most and right-most proxy/lead SNP at each locus (default=5), used to identify genes with variants in LD with lead or proxy SNPs")
)

 

opt = parse_args(OptionParser(option_list=option_list))

#############################################################################
## 1. DIRECTORIES AND FILES

folder_path <- dirname(opt$output_path)
#folder_path <- dirname("/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/model1_qnorm_hilus.bsa_final_volume/GWAS_Annotation_MRI_38k_LD/cond_stats/snps_model1_qnorm_hilus.bsa_4_18751406_19751406_chr4:18938095:T:C/4:18938095:T:C")
output_dir <- paste0(folder_path, "/GWASAnno/")
dir.create(output_dir)

# Directory for GTEx eQTL data set:
eQTLdata_dir <- "/data/public_resources/GTeX/eQTLs/V8/GTEx_Analysis_v8_eQTL"

# File containing reference genes:
gene_model_filename <- "/data/programs/bin/gwas/PRoGeM/GRCh38_genes.RData"

#Used for proxies identification
bfile=opt$bfile

#Coloc datasets
possible_eQTL_dataset <- c("GTEXv8", "Kidney_eQTL", "eQTLGen")
eQTL_datasets_coloc <- unlist(strsplit(opt$eQTL_datasets_coloc, ","))
#eQTL_datasets_coloc <- c("GTEXv8", "Kidney_eQTL", "eQTLGen") ###################################
if (!all(eQTL_datasets_coloc %in% possible_eQTL_dataset)) {
  cat("--eQTL_datasets_coloc (", eQTL_datasets_coloc,  ") doesn't correspond to an available eQTL dataset \n")
  cat("Available eQTL datasets are:", paste(possible_eQTL_dataset,collapse = ", "),  "\n")
  stop("Aborting pipeline due to invalid eQTL dataset.")
}
possible_pQTL_dataset <- c("ARIC_pGWAS", "Icelanders_pGWAS","UKB_PPP_EUR")
pQTL_datasets_coloc <- unlist(strsplit(opt$pQTL_datasets_coloc, ","))
#pQTL_datasets_coloc<- c("Icelanders_pGWAS", "UKB_PPP_EUR") #####################
if (!all(pQTL_datasets_coloc %in% possible_pQTL_dataset)) {
  cat("--pQTL_datasets_coloc (", pQTL_datasets_coloc,  ") doesn't correspond to an available pQTL dataset \n")
  cat("Available pQTL datasets are:", paste(possible_pQTL_dataset,collapse = ", "),  "\n")
  stop("Aborting pipeline due to invalid pQTL dataset.")
}

datasets_coloc=c(eQTL_datasets_coloc, pQTL_datasets_coloc)
print("datsets used for coloc")
print(datasets_coloc)
#tissues of interest - needs a couple of CHECKS!!!!
#eQTL_tissues_interest_coloc <- c("Kidney_Cortex", "Liver", "Whole_Blood", "Kidney_eQTL.TubsigeQTLs", "Kidney_eQTL.GlomsigeQTLs", "Kidney_eQTL_Meta_S686_Significant.q0.01", "CXTubsigeQTLs", "CXGlomsigeQTLs") ###################
tissues_interest <- unlist(strsplit(opt$eQTL_tissues_interest_coloc, ","))

if(is.na(opt$coloc_input_path)){
    coloc_path <- paste0(folder_path, "/coloc/output/")   
} else {
    files <- list.files(opt$coloc_input_path, pattern = "\\.RDS$")
    files <- files[!grepl("summary\\.RDS$", files)]

    if (all(paste0(datasets_coloc, "_annot_unfilt.RDS") %in% files)) {
      cat("Coloc folder provided contains all selected QTL datasets \n")
    } else {
      cat("Coloc folder provided does not contain the selected QTL datasets \n")
      stop("Aborting pipeline")
    }
}


#------------------------------------------------------------------------------------------------------
## 2. PARAMETERS
#Significance threshold
GWAS_max_nlog10P_thresh <- opt$GWAS_max_nlog10P_thresh

#Sumstats type
sumstats_type  <- opt$sumstats_type

#threshold for r2 values:
# default is 0.8.
r2_thresh <- opt$r2_thresh

# Genomic interval (in kb) either side of the lead SNP:
# default is 500kb.
interval_kb <- opt$interval_window_kb

# Genomic interval (in kb) of the overhang of the left-most and right-most proxy/lead variant at
# each locus:
LD_region_overhang_kb <- opt$LD_region_overhang_kb	# default is 5kb.

# Number of nearest genes that reside nearest to the lead variant:
number_of_nearest <- 3

# Biotype(s) of candidate genes:
biotype_of_interest <- "protein_coding"

# Column indices in the VEP output that contain the following information:
lead_rsID_column <- 40
proxy_rsID_column <- 1
ensembl_gene_id_column <- 4
IMPACT_column <- 14

# Tissue(s) of interest for GTEx eQTL data set:
GTEx_tissues <- dir(eQTLdata_dir)[grep("signif_variant_gene_pairs", dir(eQTLdata_dir))]
tissues_eQTL_associations <- GTEx_tissues     # all tissues.

# Column indices and threshold in the COLOC eQTL file
coloc_eqtl_lead_rsID_col <- 1
# lead of Ensembl gene ID column
coloc_eqtl_ensembl_gene_id_col <- 7
# lead of the tissue column
eqtl_tissue_col <- 6

# Only these tissues will be considered for eQTL_tissues_interest_coloc
# Reformat the tissues_eQTL_associations from GTEx eQTL lookup above and include
tissues_interest.short <-sub(pattern=".v8.signif_variant_gene_pairs.txt.gz",
                                replacement="", x=GTEx_tissues, fixed=TRUE)
eqtl_tissues_of_interest <- c(tissues_interest.short,
                              "Kidney_eQTL.GlomsigeQTLs", "Kidney_eQTL_Meta_S686_Significant.q0.01",
                              "Kidney_eQTL.TubsigeQTLs", "CXTubsigeQTLs", "CXGlomsigeQTLs","eQTLGen_blood")

#Check if eQTL_tissues_interest_coloc selected exist in the database
if (all(is.na(tissues_interest))) {
    cat("No tissues selected (NA detected)\n")
} else if (all(tissues_interest %in% eqtl_tissues_of_interest)) {
    cat("Valid tissue selection:\n")
    print(tissues_interest)
} else {
    cat("--eQTL_tissues_interest_coloc selection is not valid \n")
    cat("Available eQTL_tissues_interest_coloc are:", paste(eqtl_tissues_of_interest, collapse = ", "), "\n")
    stop("Aborting pipeline due to eQTL_tissues_interest_coloc selection")
}


eqtl_sumstats_2_max_nlog10P_col <- 4
# If set, sumstats_2_max_nlog10P will be filtered for > this threshold
eqtl_sumstats_2_max_nlog10P_thresh <- opt$QTL_coloc_max_nlog10P_thresh
eqtl_PP.H4.abf_col <- 5
# If set, PP.H4.abf will be filtered for > this threshold
eqtl_PP.H4.abf_thresh <- opt$PP.H4.abf_thresh
# lead of the gene type column
coloc_eqtl_gene_type_col <- 8
eqtl_cis_trans_col <- 10
eqtl_cis_trans_sel <- c("cis")

# Column indices and threshold in the COLOC pQTL file
coloc_pqtl_lead_rsID_col <- 1
coloc_pqtl_ensembl_gene_id_col <- 6
pqtl_sumstats_2_max_nlog10P_col <- 4
pqtl_sumstats_2_max_nlog10P_thresh <- opt$QTL_coloc_max_nlog10P_thresh
pqtl_PP.H4.abf_col <- 5
pqtl_PP.H4.abf_thresh <- opt$PP.H4.abf_thresh
pqtl_cis_trans_col <- 7
pqtl_cis_trans_sel <- c("cis")

#------------------------------------------------------------------------------------------------------
# 4. EXECUTE ALL

cat("##################################################\n######## STEP1: Performing pre-processing ########\n ################################################## \n")
start_time1 <- Sys.time()
output_path <- opt$output_path
GWAS_RDS <- opt$GWAS_RDS

#source("/data/programs/pipelines/GWASannotation/00_scripts/01_preprocessing.R")
source(paste0(script_path,"/01_preprocessing.R"))
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
  file.exists(paste0(output_path,"_lead.txt")) &
  file.exists(paste0(output_path,"_proxies.txt")) &
  file.exists(paste0(output_path,"_proxies_vep.txt"))) {
  cat("##################################################\n################# STEP2: run VEP ################# \n################################################## \n")
  input_vep=paste0(output_path,"_proxies_vep.txt")
  start_time2 <- Sys.time()
  #system(paste0("/data/programs/pipelines/GWASannotation/00_scripts/02_vep.sh ",input_vep))
  system(paste0(script_path,"/02_vep.sh ",input_vep))
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
    cat("########################################################\n###### STEP3: Performing eQTL and pQTL colocalization ######\n################################################## \n")
    start_time3 <- Sys.time()
    source(paste0(script_path,"/03_coloc.R"))

    end_time3 <- Sys.time()
    execution_time_seconds <- as.numeric(difftime(end_time3, start_time3, units = "secs"))
    execution_time_minutes <- execution_time_seconds / 60
    cat(sprintf("STEP3 done; took %.2f minutes\n", execution_time_minutes), "\n\n")
    rm(execution_time_seconds)
    rm(execution_time_minutes)
} else {
    cat("coloc loaded from ", opt$coloc_input_path, "\n")
    coloc_path <- opt$coloc_input_path
}

if (file.exists(paste0(coloc_path,"/",datasets_coloc[1], "_annot_unfilt.RDS"))) {
    cat("########################################################\n###### STEP4: Prepare coloc QTL results for Annotation ###### \n######################################################## \n")
    start_time6 <- Sys.time()
    source(paste0(script_path,"/04_prepare_QTL.R"))
    cat("Step6 done\n\n")
    end_time6 <- Sys.time()
    execution_time_seconds <- as.numeric(difftime(end_time6, start_time6, units = "secs"))
    # Convert seconds to minutes
    execution_time_minutes <- execution_time_seconds / 60
    # Print the execution time
    cat(sprintf("STEP4 took %.2f minutes\n", execution_time_minutes), "\n\n")
    rm(execution_time_seconds)
    rm(execution_time_minutes)
    } else {
    stop("Output file of 03_coloc.R does not exist. Aborting the pipeline.")
}

if (file.exists(paste0(coloc_path,"input_Anno_eQTL.txt")) &
file.exists(paste0(coloc_path,"input_Anno_pQTL.txt"))) {
    #Files needed to run annotation
    #Coloc
    COLOC_EQTL_filename=paste0(coloc_path,"input_Anno_eQTL.txt")
    COLOC_PQTL_filename=paste0(coloc_path,"input_Anno_pQTL.txt")

    # Files containing GWAS lead and proxy variants:
    lead_filename <- paste0(output_path,"_lead.txt")
    proxy_filename <- paste0(output_path,"_proxies.txt")

    # File containing VEP annotation:
    VEP_filename <- paste0(output_path,"_proxies_vepout.txt")

    # File containing a PoPS results:#full path
    filename_PoPS <- paste0(output_path,"_PoPS.RData")

    cat("##################################################\n############### STEP5: Run Annotation ################ \n################################################## \n")
    start_time5 <- Sys.time()
    cat("source 05_run_Anno.R \n")
    source(paste0(script_path,"/05_run_Anno.R"))
    end_time5 <- Sys.time()
    execution_time_seconds <- as.numeric(difftime(end_time5, start_time5, units = "secs"))
    execution_time_minutes <- execution_time_seconds / 60
    cat(sprintf("STEP5 took %.2f minutes\n", execution_time_minutes), "\n\n")
    rm(execution_time_seconds)
    rm(execution_time_minutes)
    } else {
    stop("Output file of 04_prepare_QTL.R does not exist. Aborting the pipeline.")
}

if (file.exists(paste0(output_dir,"OUTPUT_anno_summary.txt"))) {
  cat("##################################################\n###### STEP6: Pos-processing Annotation results (include scoring) ###### \n################################################## \n")
  start_time6 <- Sys.time()
  rm(opt)

  source(paste0(script_path,("/06_postprocessing.R")))

  end_time6 <- Sys.time()
  execution_time_seconds <- as.numeric(difftime(end_time6, start_time6, units = "secs"))
  # Convert seconds to minutes
  execution_time_minutes <- execution_time_seconds / 60
  # Print the execution time
  cat(sprintf("STEP8 took %.2f minutes\n", execution_time_minutes), "\n\n")
  rm(execution_time_seconds)
  rm(execution_time_minutes)
  } else {
  stop("Output file of 05_run_Anno.R does not exist. Aborting the pipeline.")
}

end_time_main <- Sys.time()
execution_time_seconds <- as.numeric(difftime(end_time_main, start_time_main, units = "secs"))
# Convert seconds to minutes
execution_time_minutes <- execution_time_seconds / 60
# Print the execution time
cat(sprintf("Pipeline took %.2f minutes\n", execution_time_minutes), "\n\n")
#rm(execution_time_seconds)
#rm(execution_time_minutes)
#system("rm $TEMP_SYMLINK")

cat("##################################################\n############### Pipeline complete! ############### \n################################################## \n")
