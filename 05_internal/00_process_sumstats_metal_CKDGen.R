###############################################################
#
# Preparing input RDS file for GWASannotation pipeline
# Supposed to work for Metal output files from CKDGen R5
# 
# Need to adapt: in multi file the N column is n_total_sum_sum
# otherwise it is n_total_sum
#
# There might be other places to adapt as well
#
###############################################################

suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
library(Rmpfr)
library(parallel)
devtools::load_all("/data/programs/pipelines/genepicoloc/genepicoloc_package")

option_list = list(
  make_option("--GWAS", action="store", default=NA, type='character', help="GWAS summary stats - Metal output (gz and tabix) [required]"),
  make_option("--genome_build", action="store", default="hg38", type='character', help="genome build of the GWAS summary stats (‘hg37‘ or ’hg38‘) [default=hg38]"),
  make_option("--output_path", action="store", default=NA, type='character', help="output folder path [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

# opt$GWAS <- "../../uacr_int_20241014_sumstat_output/uacr_int_20241014.EUR.tbl.gz"
# opt$genome_build <- "hg38"
# opt$output_path <- "uacr_int_20241014_EUR"
# opt$GWAS <- "/data/meta_analyses/00_CKDGen/01_analyses/09_GWAS-r5/urate-Sahar/urate-multi/urate_20250103.multi.tbl.gz"

#####################################
# load input file
sumstats <- data.table::fread(opt$GWAS)
cat("number of SNPs orginal file", nrow(sumstats), "\n")

cat("head gwas\n")
print(head(sumstats))

sumstats$nlog10_p <- -log10(as.numeric(sumstats$`P-value`))
# Recalculate -log10 for small p-values (< 1e-300) using mpfr in parallel
threshold <- 1e-300
small_pvals <- which(as.numeric(sumstats$`P-value`) < threshold)
print(length(small_pvals))

sumstats$pval=as.numeric(sumstats$`P-value`)
cat("number of SNPs with pval zero", nrow(sumstats[sumstats$pval==0, ]), "\n")

num_cores <- detectCores() - 1
# Perform mpfr-based -log10 calculation for small p-values
mpfr_results <- unlist(mclapply(sumstats$`P-value`[small_pvals], function(x) {
  p_val_mpfr <- mpfr(x, precBits = 300)
  log_val <- -log10(p_val_mpfr)
  as.numeric(log_val)  # Convert mpfr result to numeric!!!!
}, mc.cores = num_cores))
sumstats$nlog10_p[small_pvals] <- mpfr_results

print("summary nlog10_p")
print(summary(sumstats$nlog10_p))

# Change chr 23 to X, to align with the LD ref Plink files
# sumstats$NEWID = sub(pattern="^23:", replacement="X:", x=sumstats$NEWID)
sumstats$MarkerName =sub(pattern="^23:",replacement="X:",x=sumstats$MarkerName)
sumstats$Chr[sumstats$Chr==23] = "X"

genome_build=opt$genome_build  #add as parameter
cat("genome build: ", genome_build, "\n")

# Need to adapt here:
# In multi file the N column is n_total_sum_sum, otherwise n_total_sum

sumstats_1 <- read_sumstats(sumstats,
                        Name = "MarkerName",
                        rsID = "MarkerName",
                        CHR = "Chr",
                        POS = "Pos_b38",
                        A1 = "Allele1", # re-check
                        A2 = "Allele2", # re-check
                        BETA = "Effect",
                        SE = "StdErr",
                        nlog10p_value = "nlog10_p",
                        AF = "Freq1",
                        N = "n_total_sum_sum") # sometimes its n_total_sum

# print(colnames(sumstats_1))
# sumstats_1$Name <- sub(pattern="_", replacement=":", x=sumstats_1$Name,
#     fixed=TRUE)
# sumstats_1$Name <- sub(pattern="/", replacement=":", x=sumstats_1$Name,
#     fixed=TRUE)
# sumstats_1$Name <- paste0("chr", sumstats_1$Name)
# print(head(sumstats_1))

print(head(sumstats_1))
colnames(sumstats_1)[2] <- "rsID"
print(colnames(sumstats_1))

# filter for one chr - testing
# sumstats_1 <- sumstats_1[sumstats_1$CHR==19, ]
print(dim(sumstats_1))

# Use Name_by_position function to get the Name_hg38 in correct format
# Note: variants not in dbSNP will be lost. This is okay

sumstats_name_1 <- Name_by_position(sumstats=sumstats_1,
                                    CHR_name="CHR",
                                    POS_name="POS",
                                    A1_name="A1",
                                    A2_name="A2",
                                    tabix_bin="tabix",
                                    dbSNP_file="/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz",
                                    mc_cores = 16)

print(dim(sumstats_name_1))
print(head(sumstats_name_1))

# use Name_hg38 for columns Name, it is in format CHR:BP:REF:ALT
# rsID stays MarkerName from input
sumstats_name_1$Name <- sumstats_name_1$Name_hg38
sumstats_name_1$unique_ID <- sumstats_name_1$Name_hg38 <-  NULL

# re-order columns
sumstats_name_1 <- sumstats_name_1[, c("Name", "rsID", "CHR", "POS", "A1",
      "A2", "BETA", "SE", "nlog10P", "AF", "N", "rs")]

# Duplicate Name is not allowed in the coloc, need to remove them
# There should not be any by doing the above steps

if (anyDuplicated(sumstats_name_1$Name)) {
  unique_rows <- sumstats_name_1 %>%
    group_by(Name) %>%
    slice(which.max(nlog10P)) %>%
    ungroup()
} else {
  unique_rows <- sumstats_name_1
}

cat("number of SNPs after removing duplicate 'Name'", nrow(unique_rows),"\n")
saveRDS(unique_rows, paste0(opt$output_path, "_hg38_dedup.RDS"))
