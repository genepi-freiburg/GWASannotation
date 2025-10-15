###############################################################
# THIS IS NOT PART OF GWAS ANNOTATION PIPELINE You need to create the input file for GWASAnnotation
# You need to create the input file for GWASAnnotation before running the pipeline. Use this code as help only
##################################################################

suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
library(parallel)
library(Rmpfr)

devtools::load_all("/data/programs/pipelines/genepicoloc/genepicoloc_package")


option_list = list(
  make_option("--GWAS", action="store", default=NA, type='character', help="GWAS summary stats - snptest .gwas"),
  make_option("--genome_build", action="store", default="hg38", type='character', help="genome build of the GWAS summary stats (‘hg37‘ or ’hg38‘) [default=hg38]"),
  make_option("--output_path", action="store", default=NA, type='character', help="output folder path [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

#####################################
# load input file
# load input file
sumstats <- read.table(opt$GWAS, header=T,colClasses = "character")
#sumstats <- read.table("/data/studies/00_GCKD/01_analyses/gwas/multi_biomarker_GWAS/common_chip/21_VAE/02_output/modules/ME30_ext_residuals/unadjusted/GCKD_Common_Clean_dedup.gwas", header=T)
sumstats$nlog10 <- -log10(as.numeric(sumstats$pval))

#Recalculate -log10 for small p-values (< 1e-300) using mpfr in parallel
threshold <- 1e-300
small_pvals <- which(as.numeric(sumstats$pval) < threshold)
length(small_pvals)

num_cores <- detectCores() - 1
# Perform mpfr-based -log10 calculation for small p-values
mpfr_results <- unlist(mclapply(sumstats$pval[small_pvals], function(x) {
  p_val_mpfr <- mpfr(x, precBits = 300)
  log_val <- -log10(p_val_mpfr)
  as.numeric(log_val)  # Convert mpfr result to numeric!!!!
}, mc.cores = num_cores))
sumstats$nlog10[small_pvals] <- mpfr_results
sumstats$nlog10 <- round(sumstats$nlog10, 3) #Round results

sumstats[which.max(sumstats$nlog10), ]


cat("number of SNPs orginal file", nrow(sumstats))
print("head gwas")
head(sumstats)
genome_build=opt$genome_build  #add as parameter
cat("genome build: ", genome_build,  "\n")
if(genome_build=="hg37"){
    cat("\n## LiftOver sumstats to hg38 (using genepi_liftOver function) ##  \n")
    sumstats_liftOver <- genepi_liftOver(sumstats, CHR_name = "chr", POS_name ="position", A1_name= "coded_all", A2_name= "noncoded_all",
                                liftOver_bin = "/scratch/global/martins/liftover/liftOver", liftOver_chain_hg19ToHg38 = "/scratch/global/martins/liftover/hg19ToHg38.over.chain.gz", dbSNP_file="/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz", tabix_bin="tabix",
                                mc_cores=4, keep_lower=F, do_soring=T, rm_tmp_liftOver=T)
   print(head(sumstats_liftOver))
   sumstats_liftOver$beta <- as.numeric(sumstats_liftOver$beta)
   sumstats_liftOver$SE <- as.numeric(sumstats_liftOver$SE)  
   sumstats_liftOver$nlog10 <- as.numeric(sumstats_liftOver$nlog10)
   sumstats_liftOver$POS_hg38 <- as.numeric(sumstats_liftOver$POS_hg38)
   
   sumstats_1 <- read_sumstats(sumstats_liftOver,
                               Name="Name_hg38",
                               rsID = "SNP",
                               CHR = "CHR_hg38",
                               POS = "POS_hg38",
                               A1 = "coded_all", #re-check
                               A2 = "noncoded_all", #re-check
                               BETA = "beta",
                               SE = "SE",
                               nlog10p_value = "nlog10",
                               AF = "AF_coded_all")
  sumstats_1$N=unique(sumstats$n_total) #N is needed for magma input
  #### save hg38 file
  cat("number of SNPs in hg38", nrow(sumstats_1))
  
  # Find duplicated names
  unique_rows <- sumstats_1 %>%
    group_by(Name) %>%
    slice(which.max(nlog10P)) %>%
    ungroup()
  cat("number of SNPs after removing duplicated 'Name'", nrow(unique_rows))
  
  write.table(sumstats_1, paste0(opt$output_path,"_liftOver_hg38.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  saveRDS(sumstats_1, paste0(opt$output_path, "_liftOver_hg38.RDS"))
  saveRDS(unique_rows, paste0(opt$output_path, "_liftOver_hg38_dedup.RDS"))

}else{
    #TO DO:
}
