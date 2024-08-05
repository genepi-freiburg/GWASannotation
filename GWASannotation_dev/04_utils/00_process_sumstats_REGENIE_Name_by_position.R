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
#A1 using alternative allele (ALLELE1) and #A2 using reference allele (ALLELE0)
if(genome_build=="hg37"){
    cat("\n## LiftOver sumstats to hg38 (using genepi_liftOver function) ##  \n")
    sumstats_liftOver <- genepi_liftOver(sumstats, CHR_name = "CHROM", POS_name ="GENPOS", A1_name= "ALLELE1", A2_name= "ALLELE0",
                                liftOver_bin = "/scratch/global/martins/liftover/liftOver", liftOver_chain_hg19ToHg38 = "/scratch/global/martins/liftover/hg19ToHg38.over.chain.gz", dbSNP_file="/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz", tabix_bin="tabix",
                                unique_ID_name="unique_ID",
                                mc_cores=8, keep_lower=F, do_soring=T, rm_tmp_liftOver=T, do_Name_by_position=F)

   sumstats_liftOver$Name_hg38 <- paste0("chr", sumstats_liftOver$CHR_hg38, ":", sumstats_liftOver$POS_hg38,":", sumstats_liftOver$A2_hg38, ":",sumstats_liftOver$A1_hg38)
   print(head(sumstats_liftOver))
   
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
                                A1 = "ALLELE1", #alternative
                                A2 = "ALLELE0", #reference
                                BETA = "BETA",
                                SE = "SE",
                                nlog10p_value = "LOG10P",
                                AF = "A1FREQ")
    sumstats_1$N=sumstats$N #N is needed for magma input
    
    
    #sumstats_2 <- Name_by_position(sumstats_1, tmp_name=NULL,
    #                             CHR_name="CHR", POS_name="POS",
    #                             A1_name="A1", A2_name="A2",
    #                             Name_out="Name_new", rs_name="ID",
    #                             unique_ID_name="unique_ID",
    #                             tabix_bin="tabix", dbSNP_file="/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz",
    #                             do_soring=T, mc_cores=8)
    #sumstats_2$N=sumstats$N
    print(head(sumstats_1)
    # Find duplicated names
    unique_rows <- sumstats_1 %>%
      group_by(Name) %>%
      slice(which.max(nlog10P)) %>%
      ungroup()
    cat("number of SNPs after removing duplicated 'Name'", nrow(unique_rows))
    
    write.table(sumstats_1, paste0(opt$output_path,".txt"), row.names=F, col.names=T, sep="\t", quote=F)
    saveRDS(sumstats_1, paste0(opt$output_path, ".RDS"))
    saveRDS(unique_rows, paste0(opt$output_path, "_dedup.RDS"))
}



#EXTRA
bfile="/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_MRI_38k_MAF001/plink/chrpos/hg38/UKBB_MRI_38k_MAF001_hg38_allchr.bim"
b=read.table(bfile)
c <- data.table::fread("/data/public_resources/CKDGen/preprocessing/Wuttke2019/CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt_hg38.gz")

sumstats <- data.table::fread("/dsk/data1/studies/04_ARIC/01_analyses/03_paired_mGWAS/01_input/serum/EA/ARIC_EA_TopMed_2023-05-10_C100000961_271_C100000961.regenie.gz")
length(unique(sumstats$ID))

#test bfile and coloc
#test_b=intersect(sumstats$ID, b$V2)
t_b=sumstats[which(sumstats$ID %in% b$V2),]
dim(t_b)
nrow(t_b)/length(unique(sumstats$ID))*100
#0.13
sumstats$Name_inv <- paste0("chr",sumstats$CHROM,":", sumstats$GENPOS,":", sumstats$ALLELE1,":", sumstats$ALLELE0)
t_b2=sumstats[which(sumstats$Name_inv %in% b$V2),]
dim(t_b2)
nrow(t_b2)/length(unique(sumstats$ID))*100

#coloc test
test_c=intersect(sumstats$ID, c$Name_hg38)
t_c=sumstats[which(sumstats$ID %in%  c$Name_hg38),]
dim(t_c)
nrow(t_c)/length(unique(sumstats$ID))*100
length(unique(sumstats$ID))/nrow(t_c)*100
length(unique(c$Name_hg38))
#66%
t_c2=sumstats[which(sumstats$Name_inv %in%  c$Name_hg38),]
dim(t_c2)
nrow(t_c2)/length(unique(sumstats$ID))*100
#0.014%

#kidney volume
## before liftover
sumstats <- data.table::fread("/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/10Feb2024_final_volumes/maf001/model1_qnorm_tkv.bsa_chr1-22_maf001.regenie.gz")
sumstats_2 =sumstats[which(sumstats$CHROM ==2),]
sumstats_2$Name=paste0("chr",sumstats_2$CHROM,":", sumstats_2$GENPOS,":", sumstats_2$ALLELE0,":", sumstats_2$ALLELE1)
bim=read.table("/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_MRI_38k_MAF001/plink/UKBB_MRI_38k_MAF001_chr2.bim")
bim$Name = paste0("chr",bim$V1,":", bim$V4,":",bim$V5,":", bim$V6)
t_b=sumstats_2[which(sumstats_2$Name %in% bim$Name),]
nrow(t_b)/length(unique(sumstats_2$Name))*100
#95.47%

#chrpos folder (before liftover)
rm(bim)
bim=read.table("/data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_MRI_38k_MAF001/plink/chrpos/UKBB_MRI_38k_MAF001_chr2.bim")
t_b=sumstats_2[which(sumstats_2$Name %in% bim$V2),]
nrow(t_b)/length(unique(sumstats_2$Name))*100
#95.47% -- PROBLEM is here

#after liftover
sumstats <- readRDS("/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/10Feb2024_final_volumes/maf001/model1_qnorm_tkv.bsa_chr1-22_maf001_liftOver_hg38_dedup.RDS")
sumstats$Name_inv <- paste0("chr",sumstats$CHR,":", sumstats$POS,":", sumstats$A1,":", sumstats$A2)
length(unique(sumstats$Name))

#test bfile and coloc
#test_b=intersect(sumstats$ID, b$V2)
t_b=sumstats[which(sumstats$Name %in% b$V2),]
nrow(t_b)/length(unique(sumstats$Name))*100
#19%
t_b2=sumstats[which(sumstats$Name_inv %in% b$V2),]
nrow(t_b2)/length(unique(sumstats$Name))*100
#78%

#coloc test
t_c=sumstats[which(sumstats$Name %in%  c$Name_hg38),]
dim(t_c)
nrow(t_c)/length(unique(sumstats$Name))*100
#85%
t_c2=sumstats[which(sumstats$Name_inv %in%  c$Name_hg38),]
dim(t_c2)
nrow(t_c2)/length(unique(sumstats$Name))*100
