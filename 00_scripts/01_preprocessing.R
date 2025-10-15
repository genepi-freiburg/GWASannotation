###############################################################
# process the input (regneie.gz for now) and creates input files for VEP, ProGEM script and for MAGMA
# input: .RDS
# outputs:
#  - loci regions (function bored from coloc pipeline) - intermediate file or will we use it in the end?
#  - lead file: txt file with the lead SNP (top SNP) for each loci
#  - proxy file: txt file with the snps for each loci (+/- 500kb around the lead snp) ***NOTE: add parameter to decide the window
#  - filtered regions file for coloc!
# **** take genome version into consideration????
#################################################################

#####################################
# load input file
sumstats <- readRDS(GWAS_RDS)

print("head gwas")
head(sumstats)
table(sumstats$Name == sumstats$rsID)

#To make sure that name match the coloc datasets names
sumstats_name <- Name_by_position(sumstats=sumstats,
                                   CHR_name="CHR",
                                   POS_name="POS",
                                   A1_name="A1",
                                   A2_name="A2",
                                   tabix_bin="tabix",
                                   dbSNP_file="/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz")


table(sumstats_name$Name == sumstats_name$Name_hg38)
length(unique(sumstats_name$Name_hg38))

dim(sumstats_name[duplicated(sumstats_name$Name_hg38)])
setDT(sumstats_name)
sumstats[sumstats_name, on = "rsID", Name := i.Name_hg38]
dim(sumstats[duplicated(sumstats$Name)])
sumstats <- sumstats %>%
  group_by(Name) %>%
  slice(which.max(nlog10P)) %>%
  ungroup()
print("Name == rsID")
table(sumstats$Name == sumstats$rsID)

#####################################
# loci regions - using function from coloc
cat("\n## Defining significant loci (regions around 500kb of the lead SNP) ##  \n")
regions_list <- get_coloc_regions(sumstats, nlogP_threshold = GWAS_max_nlog10P_thresh, halfwindow = 500000)

regions <- regions_list$coloc_regions
regions <- regions[which(regions$comment=="PASS"),]
cat(paste0(nrow(regions)), " loci identified \n")

regions_log <- regions_list$regions_log
sumstats_filt <- regions_list$sumstats_filt
#head(sumstats_filt)

# Save log and filtered summary statistics
writeLines(regions_log, con = paste0(output_path, "_get_coloc_regions_log.txt"))
saveRDS(sumstats_filt, paste0(output_path, "_subset.RDS"))

#Following code substitute:
#save_coloc_regions(regions_list, output_path, sumstats_1_type=sumstats_type)
write.table(regions, paste0(output_path, "_coloc_regions_PASS.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

write.table(sumstats_filt, paste0(output_path, "_subset.tsv"),
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

system(paste0("bgzip -f ", output_path, "_subset.tsv"))
system(paste0("tabix -f -s3 -b4 -e4 ", output_path, "_subset.tsv.gz -c Name"))
##

saveRDS(regions_list, paste0(output_path, "_coloc_regions.RDS"))

#####################################
# lead file
cat("\n## Creating lead file ## \n")
#3.1. A tab-separated .txt file containing rsIDs, chromosomes, and GRCh38 coordinates (both start and end) of your lead variants of interest across four columns with the below column names. In cases where there is no rsID for a lead variant then the notation "chr:start" should be used (i.e., 1:11856378).
tophit=regions

tophit$strand <- "+"

# If theres no rsID, create "chr:start" annotation
tophit<- tophit %>%
    mutate(ID = ifelse(startsWith(rsID, "rs"), rsID, paste0(CHR_var, ":", POS)))

tophit <- tophit[c("rsID","CHR_var", "POS", "POS")]
colnames(tophit) <- c("rsID", "CHR", "START", "END")

write.table(x=tophit, file = paste0(output_path, "_lead.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#####################################
# proxy file - and vep input file
cat("\n## Creating proxie file ## \n")
#3.2. A tab-separated .txt file containing rsIDs, chromosomes, and GRCh38 coordinates (both start and end) of proxies for your lead variants; along with the corresponding lead rsIDs and r2 value across six columns in total with the below column names. Your proxies can be derived either directly from your sample or, for example, from 1000 genomes data. If proxies have not been pre-filtered prior to input, then ProGeMM will filter based on a user-defined r2 threshold (default r2 is 0.8) if desired.
#PROXY_rsID    PROXY_CHR    PROXY_START    PROXY_END    LEAD_rsID    r2
#rs602950    1    20915531    20915531    rs532545    0.991
#rs1253904    1    20913519    20913519    rs532545    0.963
#rs589942    1    20916080    20916080    rs532545    0.931


#vep -input
#Default VEP input
#The default format is a simple whitespace-separated format (columns may be separated by space or tab characters), containing five required columns plus an optional identifier column:
#chromosome - just the name or number, with no 'chr' prefix
#start
#end
#allele - pair of alleles separated by a '/', with the reference allele first (or structural variant type)
#strand - defined as + (forward) or - (reverse). The strand will only be used for VEP to know which alleles to use.
#identifier - this identifier will be used in VEP's output. If not provided, VEP will construct an identifier from the given coordinates and alleles.

#! REGENIE: reference allele (allele 0), alternative allele (allele 1)

#Getting proxy SNPs based on lead.txt file previously created
#SnIPA uses hg37 - using plink for hg38
#plink2 doesnt support --ld-snp-list, need to use plink1
#if no bfile provided, will use /data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_150k_RandomSubset_Cleaned/hg38/UKBB_14k_hg38_chr1-22
system(paste0("/data/programs/bin/gwas/plink/plink-1.90_beta6.20/plink --bfile ", bfile, " --r2 with-freqs --ld-snp-list ", output_path, "_lead.txt --ld-window 1000000 --ld-window-kb 1000 --out ", output_path, "_ld_results"))
ld=read.table(paste0(output_path, "_ld_results.ld"),header=T)

r2.cutoff <- r2_thresh 

ld_results=ld[which(ld$R2>r2.cutoff),]
#sumstats_filt=readRDS(paste0(output_path, "_subset.RDS"))
sumstats_filt=as.data.frame(sumstats_filt)
#ld_results=merge(ld_results, sumstats_filt[,c("rsID", "A1","A2")], by.x="SNP_A", by.y="rsID")
#I want the ref and alt from the proxies! if the r2==1 proxie and lead is the same snp
ld_results=merge(ld_results, sumstats_filt[,c("rsID", "A1","A2")], by.x="SNP_B", by.y="rsID")

ld_results$strand <- "+"

#ld_results$ref <- sapply(ld_results$SNP_B, function(x) unlist(strsplit(x, ":"))[3])
#ld_results$alt <- sapply(ld_results$SNP_B, function(x) unlist(strsplit(x, ":"))[4])
#ld_results$allele <- paste(ld_results$ref, ld_results$alt, sep = "/") #A2 is reference
#ld_results$END =as.numeric(ld_results$BP_B) +  nchar(as.character(ld_results$ref)) - 1 #need to addapt to indels


ld_results$allele <- paste(ld_results$A2, ld_results$A1, sep = "/") #A2 is reference
ld_results$END =as.numeric(ld_results$BP_B) +  nchar(as.character(ld_results$A2)) - 1 #need to addapt to indels


proxy_data <- data.frame(
  PROXY_rsID = ld_results$SNP_B,
  PROXY_CHR = ld_results$CHR_B,
  PROXY_START = ld_results$BP_B,
  PROXY_END = ld_results$END,
  LEAD_rsID = ld_results$SNP_A,
  r2 = ld_results$R2
)
print("Dimension proxies data:")
print(dim(proxy_data))


cat("head proxie file \n")
print(head(proxy_data))
write.table(proxy_data, file= paste0(output_path, "_proxies.txt"), quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


#FOR VEP ANALYSIS NEED TO CHECK ALLELE POSITION BASED ON dbSNP - use
#Name_by_position(sumstats, tmp_name=NULL,CHR_name="CHR_hg38", POS_name="POS_hg38",A1_name="A1_hg38", A2_name="A2_hg38",Name_out="Name_hg38", rs_name="rs",unique_ID_name="unique_ID",tabix_bin, dbSNP_file,do_soring=T, mc_cores=4)

#devtools::load_all("/data/programs/pipelines/genepicoloc/genepicoloc_package")
#library(data.table)
#setDT(ld_results)
#ld_results2 <- Name_by_position(ld_results, tmp_name=NULL,CHR_name="CHR_B", POS_name="BP_B",A1_name="alt", A2_name="ref",Name_out="Name_new",unique_ID_name="unique_ID",tabix_bin="tabix", dbSNP_file="/data/public_resources/Ensembl_human_variation_b38_v109/dbSNP_v156_b38p14_rsid.vcf.gz",do_soring=T)
#ld_results2 <- separate(ld_results2, Name_new, into = c("chr_new", "pos_new", "ref_new", "alt_new"), sep = ":")

#ld_results2$allele <- paste(ld_results2$ref, ld_results2$alt, sep = "/") #A2 is reference

vep_data <- ld_results[,c("CHR_B", "BP_B", "END", "allele", "strand", "SNP_B")]
##!VEP used "X" instead of "23"
vep_data$CHR_B[vep_data$CHR_B == "23"] <- "X"

cat("head vep input file \n")
print(head(vep_data))
#vep <- separate(vep, Allele, into = c("ref", "alt"), sep = "/")

write.table(vep_data, file = paste0(output_path, "_proxies_vep.txt"), quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = FALSE)


#cat("\n## Pre-processing finished ##\n")
