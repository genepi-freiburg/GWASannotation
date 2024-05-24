###############################################################
# process the input (regneie.gz for now) and creates input files for VEP, ProGEM script and for MAGMA
# input: .RDS
# outputs:
#  - loci regions (function bored from coloc pipeline) - intermediate file or will we use it in the end?
#  - sentinel file: txt file with the index SNP (top SNP) for each loci
#  - proxy file: txt file with the snps for each loci (+/- 500kb around the index snp) ***NOTE: add parameter to decide the window
#  - magma input
#  - filtered regions file for coloc!
# **** take genome version into consideration????
#################################################################

#####################################
# load input file
sumstats <- readRDS(GWAS_RDS)

print("head gwas")
head(sumstats)


#####################################
# loci regions - using function from coloc
cat("\n## Defining significant loci (regions around 500kb of the index SNP) ##  \n")

regions_list <- get_coloc_regions(sumstats, nlogP_threshold = -log10(5e-8), halfwindow = 500000)

regions <- regions_list$coloc_regions
regions <- regions[which(regions$comment=="PASS"),]
cat(paste0(nrow(regions)), " loci identified \n")

regions_log <- regions_list$regions_log
sumstats_filt <- regions_list$sumstats_filt
head(sumstats_filt)
#sumstats_filt <- subset_sumstats(sumstats, regions)
# Save
writeLines(regions_log, con = paste0(output_path, "_get_coloc_regions_log.txt"))
saveRDS(sumstats_filt, paste0(output_path, "_subset.RDS")) #not needed?
#write.table(sumstats_filt, file = paste0(output_path, "_subset.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
#system(paste0("bgzip ", output_path, "_subset.txt"))
#system(paste0("tabix -b 4 -e 4 -S 1 -s 3 ", output_path, "_subset.txt.gz"))
save_coloc_regions(regions_list, output_path)
saveRDS(regions_list, paste0(output_path, "_coloc_regions.RDS"))

#####################################
# sentinel file
cat("\n## Creating sentinel file ## \n")
#3.1. A tab-separated .txt file containing rsIDs, chromosomes, and GRCh37 coordinates (both start and end) of your sentinel variants of interest across four columns with the below column names. In cases where there is no rsID for a sentinel variant then the notation "chr:start" should be used (i.e., 1:11856378).
tophit=regions

tophit$strand <- "+"
#tophit$rsID <- paste(tophit$CHROM, tophit$POS, sep = ":")

#tophit$allele <- paste(tophit$A1, tophit$A2, sep = "/")

# If theres no rsID, create "chr:start" annotation
tophit<- tophit %>%
    mutate(ID = ifelse(startsWith(rsID, "rs"), rsID, paste0(CHR_var, ":", POS)))

tophit <- tophit[c("rsID","CHR_var", "POS", "POS")]
colnames(tophit) <- c("rsID", "CHR", "START", "END")

write.table(x=tophit, file = paste0(output_path, "_sentinel.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#####################################
# proxy file - and vep input file
cat("\n## Creating proxie file ## \n")
#3.2. A tab-separated .txt file containing rsIDs, chromosomes, and GRCh38 coordinates (both start and end) of proxies for your sentinel variants; along with the corresponding sentinel rsIDs and r2 value across six columns in total with the below column names. Your proxies can be derived either directly from your sample or, for example, from 1000 genomes data. If proxies have not been pre-filtered prior to input, then ProGeMM will filter based on a user-defined r2 threshold (default r2 is 0.8) if desired.
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

#Getting proxy SNPs based on sentinel.txt file previously created
#SnIPA uses hg37 - using plink for hg38
#plink2 doesnt support --ld-snp-list, need to use plink1
#if no bfile provided, will use /data/studies/06_UKBB/01_Data/02_Genetic_Data/UKBB_150k_RandomSubset_Cleaned/hg38/UKBB_14k_hg38_chr1-22
system(paste0("/data/programs/bin/gwas/plink/plink-1.90_beta6.20/plink --bfile ", bfile, " --r2 with-freqs --ld-snp-list ", output_path, "_sentinel.txt --ld-window 1000000 --ld-window-kb 1000 --out ", output_path, "_ld_results"))
ld=read.table(paste0(output_path, "_ld_results.ld"),header=T)

r2.cutoff <- r2_thresh 

ld_results=ld[which(ld$R2>r2.cutoff),]
#sumstats_filt=readRDS(paste0(output_path, "_subset.RDS"))
sumstats_filt=as.data.frame(sumstats_filt)
ld_results=merge(ld_results, sumstats_filt[,c("rsID", "A1","A2")], by.x="SNP_B", by.y="rsID")
ld_results$END =as.numeric(ld_results$BP_B) +  nchar(as.character(ld_results$A1)) - 1
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



ld_results$strand <- "+"
ld_results$allele <- paste(ld_results$A2, ld_results$A1, sep = "/") #A2 is reference
vep_data <- ld_results[,c("CHR_B", "BP_B", "END", "allele",
"strand", "SNP_B")]
cat("head vep input file \n")
head(vep_data)
#dim(res)
#dim(res[!duplicated(res),])

write.table(vep_data, file = paste0(output_path, "_proxies_vep.txt"), quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = FALSE)


#####################################
# prepare magma input file
#output should look like:
#SNP P N
#rs367896724 0.810602 37941
cat("\n## Prepare magma input file ## \n")
input_magma=sumstats

input_magma$P=10^(-input_magma$nlog10P)

input_magma_rsID=input_magma[,c("rsID", "P", "N")]
colnames(input_magma_rsID)=c("SNP", "P", "N")
print(head("input magma rsID"))
head(input_magma_rsID)
write.table(input_magma_rsID,paste0(output_path,"_input_magma.txt"), row.names=F, col.names=T, sep="\t", quote=F)

#input_magma_CHRPOS=input_magma[,c("Name", "P", "N")]
#colnames(input_magma_CHRPOS)=c("SNP", "P", "N")
#print(head("input magma CHRPOS"))
#head(input_magma_CHRPOS)
#write.table(input_magma_CHRPOS,paste0(output_path,"_input_magma_CHRPOS.txt"), row.names=F, col.names=T, sep="\t", quote=F)

cat("\n## Pre-processing finished ##\n")
