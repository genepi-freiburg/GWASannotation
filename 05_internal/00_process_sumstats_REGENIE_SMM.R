###############################################################
# THIS IS NOT PART OF GWAS ANNOTATION PIPELINE You need to create the input file for GWASAnnotation
# You need to create the input file for GWASAnnotation before running the pipeline. Use this code as help only
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
#A1 using alternative allele (ALLELE1) and #A2 using reference allele (ALLELE0)
if(genome_build=="hg37"){
    cat("\n## LiftOver sumstats to hg38 (using genepi_liftOver function) ##  \n")
    sumstats_liftOver <- genepi_liftOver(sumstats, CHR_name = "CHROM", POS_name ="GENPOS", A1_name= "ALLELE1", A2_name= "ALLELE0",
                                liftOver_bin = "/scratch/global/martins/liftover/liftOver", liftOver_chain = "/scratch/global/martins/liftover/hg19ToHg38.over.chain.gz", dbSNP_file="/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz", tabix_bin="/usr/bin/tabix",
                                mc_cores=8, keep_lower=F, do_soring=T, rm_tmp_liftOver=T, do_Name_by_position=F)

   sumstats_liftOver$Name_hg38 <- paste0("chr", sumstats_liftOver$CHR_hg38, ":", sumstats_liftOver$POS_hg38,":", sumstats_liftOver$A2_hg38, ":",sumstats_liftOver$A1_hg38)
   #print(head(sumstats_liftOver))
   
   sumstats_1 <- read_sumstats(sumstats_liftOver,
                               Name="Name_hg38",
                               rsID = "ID",
                               CHR = "CHR_hg38",
                               POS = "POS_hg38",
                               A1 = "A1_hg38", #alternative ALLELE1
                               A2 = "A2_hg38", #reference ALLELE0
                               BETA = "BETA",
                               SE = "SE",
                               nlog10p_value = "LOG10P",
                               AF = "A1FREQ") #its ALLELE1 freq, no A2 FREQ now
  sumstats_1$N=unique(sumstats$N) #N is needed for magma input
 
 #keep only chr 1-22
 cat("keeping only chr 1-22 \n ")
 sumstats_1 <- sumstats_1[sumstats_1$CHR %in% 1:22, ]

 
  #### save hg38 file
  cat("number of SNPs in hg38", nrow(sumstats_1),"\n ")
  #Name is for coloc - use inverted name
  #rsID is for eveything else - use normal name that fits ld ref
  sumstats_1$rsID <- sumstats_1$Name
  sumstats_1$Name <- paste0("chr", sumstats_1$CHR, ":", sumstats_1$POS,":", sumstats_1$A1, ":",sumstats_1$A2)
  print(head(sumstats_1))
  
  # Find duplicated names
  unique_rows <- sumstats_1 %>%
    group_by(Name) %>%
    slice(which.max(nlog10P)) %>%
    ungroup()
  cat("number of SNPs after removing duplicated 'Name'", nrow(unique_rows))
  
  write.table(sumstats_1, paste0(opt$output_path,"_liftOver_hg38.txt"), row.names=F, col.names=T, sep="\t", quote=F)
  saveRDS(sumstats_1, paste0(opt$output_path, "_liftOver_hg38.RDS"))
  cat("RDS saved in ", paste0(opt$output_path, "_liftOver_hg38_dedup.RDS"))
  saveRDS(unique_rows, paste0(opt$output_path, "_liftOver_hg38_dedup.RDS"))

}else{
    sumstats$Name_inv <- paste0("chr", sumstats$CHROM, ":", sumstats$GENPOS, ":", sumstats$ALLELE1, ":", sumstats$ALLELE0)
    sumstats_1 <- read_sumstats(sumstats,
                                Name="ID",
                                rsID = "Name_inv",
                                CHR = "CHROM",
                                POS = "GENPOS",
                                A1 = "ALLELE1", #alternative
                                A2 = "ALLELE0", #reference
                                BETA = "BETA",
                                SE = "SE",
                                nlog10p_value = "LOG10P",
                                AF = "A1FREQ")
    sumstats_1$N=sumstats$N #N is needed for magma input
    
    #
    # Find duplicated names
    unique_rows <- sumstats_1 %>%
      group_by(Name) %>%
      slice(which.max(nlog10P)) %>%
      ungroup()
    cat("number of SNPs after removing duplicated 'Name'", nrow(unique_rows))
    
    
    write.table(sumstats_1, paste0(opt$output_path,"_hg38.txt"), row.names=F, col.names=T, sep="\t", quote=F)
    saveRDS(sumstats_1, paste0(opt$output_path, "_hg38.RDS"))
    saveRDS(unique_rows, paste0(opt$output_path, "_dedup.RDS"))
}
