###############################################################
# process the input (regneie.gz for now) and creates input files for VEP, ProGEM script and for MAGMA
# input: GWAS.regenie.gz
# outputs:
#  - loci regions (function bored from coloc pipeline) - intermediate file or will we use it in the end?
#  - sentinel file: txt file with the index SNP (top SNP) for each loci
#  - proxy file: txt file with the snps for each loci (+/- 500kb around the index snp) ***NOTE: add parameter to decide the window
#  - magma input
#  - filtered regions file for coloc!
# **** take genome version into consideration????
##################################################################

suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
devtools::load_all("/data/programs/pipelines/genepicoloc/genepicoloc_package")
sapply(list.files("/data/programs/pipelines/genepicoloc/custom_scripts/source", full.names = T), source)


option_list = list(
  make_option("--GWAS", action="store", default=NA, type='character', help="GWAS summary stats - regenie output (gz and tabix) [required]"),
  make_option("--genome_build", action="store", default="hg38", type='character', help="genome build of the GWAS summary stats (‘hg37‘ or ’hg38‘) [default=hg38]"),
  make_option("--output_path", action="store", default=NA, type='character', help="output folder path [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

#####################################
# load input file
sumstats <- data.table::fread(opt$GWAS)
cat("number of SNPs orginal file", nrow(sumstats))
print("head gwas")
head(sumstats)
genome_build=opt$genome_build  #add as parameter
cat("genome build: ", genome_build,  "\n")
if(genome_build=="hg37"){
    cat("\n## LiftOver sumstats to hg38 (using genepi_liftOver function) ##  \n")
    sumstats_liftOver <- genepi_liftOver(sumstats, CHR_name = "CHROM", POS_name ="GENPOS", A1_name= "ALLELE0", A2_name= "ALLELE1",
                                liftOver_bin = "/scratch/global/martins/liftover/liftOver", liftOver_chain_hg19ToHg38 = "/scratch/global/martins/liftover/hg19ToHg38.over.chain.gz", dbSNP_file="/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz", tabix_bin="tabix",
                                unique_ID_name="unique_ID",
                                mc_cores=4, keep_lower=F, do_soring=T, rm_tmp_liftOver=T)
   print(head(sumstats_liftOver))
   
   sumstats_1 <- read_sumstats(sumstats_liftOver,
                               Name="Name_hg38",
                               rsID = "ID",
                               CHR = "CHR_hg38",
                               POS = "POS_hg38",
                               A1 = "ALLELE0", #re-check
                               A2 = "ALLELE1", #re-check
                               BETA = "BETA",
                               SE = "SE",
                               nlog10p_value = "LOG10P",
                               AF = "A1FREQ")
  sumstats_1$N=unique(sumstats$N) #N is needed for magma input
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
    sumstats$Name <- paste(sumstats$CHROM, sumstats$GENPOS, sumstats$ALLELE1, sumstats$ALLELE0, sep = ":")
    sumstats_1 <- read_sumstats(sumstats,
                                Name="Name",
                                rsID = "ID",
                                CHR = "CHROM",
                                POS = "GENPOS",
                                A1 = "ALLELE0", #re-check
                                A2 = "ALLELE1", #re-check
                                BETA = "BETA",
                                SE = "SE",
                                nlog10p_value = "LOG10P",
                                AF = "A1FREQ")
    sumstats_1$N=sumstats$N #N is needed for magma input
    saveRDS(sumstats_1, paste0(opt$output_path, "_hg38.RDS"))
}
