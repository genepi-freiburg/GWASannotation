suppressMessages(library(optparse))

option_list = list(
  make_option("--sentinel_filepath", action="store", default=NA, type='character', help="sentinel_filepath [required]"),
  make_option("--proxy_filepath", action="store", default=NA, type='character', help="proxy_filepath [required]"),
  make_option("--VEP_filepath", action="store", default=NA, type='character', help="VEP_filepath [required]"),
  make_option("--PoPS_filepath", action="store", default=NA, type='character', help="PoPS_filepath [required]"),
  make_option("--coloc_dir", action="store", default=NA, type='character', help="Colocalization results directory path [required]"),
  make_option("--output_path", action="store", default=NA, type='character', help="output_path [required]")

)

opt = parse_args(OptionParser(option_list=option_list))

######################################## USER-DEFINED SETTINGS ########################################
## 1. DIRECTORIES AND FILES

# Working directory:
#base_dir <- "/data/meta_analyses/13_UMOD/04_meta/01_UMOD_overall/transethic/meta_plus_UKB50k/PRoGem_run/"

# Output directory: 

folder_path <- dirname(opt$output_path)
output_dir <- paste0(folder_path, "/ProGeM/")
dir.create(output_dir)

# Directory for GTEx eQTL data set:
eQTLdata_dir <- "/data/public_resources/GTeX/eQTLs/V8/GTEx_Analysis_v8_eQTL"

# File containing a PoPS results:
#full path
filename_PoPS=opt$PoPS_filepath

#Coloc
coloc_dir=opt$coloc_dir
COLOC_EQTL_filename=paste0(coloc_dir,"input_Progem_eQTL.txt")
COLOC_PQTL_filename=paste0(coloc_dir,"input_Progem_pQTL.txt")

# File containing reference genes:
gene_model_filename <- "/data/programs/bin/gwas/PRoGeM/GRCh38_genes.RData"							# default provided.

# Files containing GWAS sentinel and proxy variants:
# sentinel_filename #full path
sentinel_filename <- opt$sentinel_filepath

# proxy_filename #full path
proxy_filename <- opt$proxy_filepath

# File containing VEP annotation:
# VEP_filename #full path
VEP_filename <- opt$VEP_filepath

#------------------------------------------------------------------------------------------------------
## 2. PARAMETERS FOR TOP-DOWN APPROACH

# Set filter for proxy variants based on r2 values:
# filtering_required <- TRUE
filtering_required <- FALSE
		# TRUE or FALSE.

# If TRUE, set threshold for r2 values:
r2_thresh <- 0.8										# default is 0.8.

# Genomic interval (in kb) either side of the sentinel SNP:
interval_kb <- 500										# default is 500kb.

#------------------------------------------------------------------------------------------------------
## 3. PARAMETERS FOR BOTTOM-UP APPROACH

# Genomic interval (in kb) of the overhang of the left-most and right-most proxy/sentinel variant at 
# each locus:
LD_region_overhang_kb <- 5									# default is 5kb.

# Number of nearest genes that reside nearest to the sentinel variant:
number_of_nearest <- 3										# default is 3.

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

tissues_of_interest <- GTEx_tissues[]     # default is all tissues.
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
eqtl_sumstats_2_max_nlog10P_thresh <- -log10(5e-8)
eqtl_PP.H4.abf_col <- 5
# If set, PP.H4.abf will be filtered for > this threshold
eqtl_PP.H4.abf_thresh <- 0.8
# Index of the gene type column
coloc_eqtl_gene_type_col <- 8
eqtl_cis_trans_col <- 10
eqtl_cis_trans_sel <- c("cis")

# Column indices and threshold in the COLOC pQTL file
coloc_pqtl_sentinel_rsID_col <- 1
coloc_pqtl_ensembl_gene_id_col <- 6
pqtl_sumstats_2_max_nlog10P_col <- 4
pqtl_sumstats_2_max_nlog10P_thresh <- -log10(5e-8)
pqtl_PP.H4.abf_col <- 5
pqtl_PP.H4.abf_thresh <- 0.8
pqtl_cis_trans_col <- 7
pqtl_cis_trans_sel <- c("cis")

#------------------------------------------------------------------------------------------------------
## 4. EXECUTE ANNOTATION

source(file = file.path("/data/programs/pipelines/GWASannotation/00_scripts/ProGeM_functions.R"))			# provided.
source(file = file.path("/data/programs/pipelines/GWASannotation/00_scripts/ProGeM_commands.R"))	    # provided.

#######################################################################################################
