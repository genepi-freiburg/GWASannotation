###############################################################
# process the input (regneie.gz for now) and creates input files for VEP, ProGEM script and for MAGMA
# input: GWAS.regenie.gz
# outputs:
#  - loci regions (function bored from coloc pipeline) - intermediate file or will we use it in the end?
#  - sentinel file: txt file with the index SNP (top SNP) for each loci
#  - proxy file: txt file with the snps for each loci (+/- 500kb around the index snp) ***NOTE: add parameter to decide the window
#  - magma input
# **** also preprocess file for coloc!
# **** take genome version into consideration????
##################################################################


suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(optparse))
suppressMessages(library(data.table))

option_list = list(
  make_option("--GWAS", action="store", default=NA, type='character', help="GWAS summary stats - regenie output (gz and tabix) [required]"),
  make_option("--output_folder", action="store", default=NA, type='character', help="output folder path [required]")
)

opt = parse_args(OptionParser(option_list=option_list))

#####################################
# load input file
input=fread(opt$GWAS)
input=fread("/data/studies/06_UKBB/02_Projects/14_MRI-kidney/02_output/regenie_output/step2/20Oct/maf001/model2_qnorm_tkv.bsa_cov_eGFR_chr1-22_maf001.regenie.gz")



#create a p-value collumn P
input$P=10^(-input$LOG10P)


print("head gwas")
head(input)



#####################################
# loci regions - using function from coloc
cat("Defining significant loci (regions aoung 500kb of the index SNP \n")
source("/data/programs/pipelines/genepicoloc/genepicoloc_package/R/02_process_sumstats.R") #*** when adding coloc preprocessing, just load all...

regions_list <- get_coloc_regions(sumstats = input, CHR_name = "CHROM", POS_name = "GENPOS", p_value_name = "P",
                                        p_threshold = 5e-8,
                                        halfwindow = 500000)
regions <- regions_list$coloc_regions_PASS
cat(paste0(nrow(regions)), " loci identified \n")


# Save
saveRDS(regions_list, paste0(opt$output,"loci_regions.RDS"))

#####################################
# sentinel file
cat("Creating sentinel file \n")
#3.1. A tab-separated .txt file containing rsIDs, chromosomes, and GRCh37 coordinates (both start and end) of your sentinel variants of interest across four columns with the below column names. In cases where there is no rsID for a sentinel variant then the notation "chr:start" should be used (i.e., 1:11856378).
tophit=regions

tophit$strand <- "+"
#tophit$rsID <- paste(tophit$CHROM, tophit$GENPOS, sep = ":")

tophit$allele <- paste(tophit$ALLELE0, tophit$ALLELE1,
    sep = "/")

tophit<- tophit %>%
    mutate(ID = ifelse(startsWith(ID, "rs"), ID, paste0(CHROM, ":", GENPOS)))

tophit <- tophit[c("ID","CHROM", "GENPOS", "GENPOS")]
colnames(tophit) <- c("rsID", "CHR", "START", "END")

write.table(x=tophit, file = paste0(opt$output_folder, "/sentinel.txt"),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

#####################################
# proxy file - and vep input file
cat("Creating proxie file \n")
#3.2. A tab-separated .txt file containing rsIDs, chromosomes, and GRCh37 coordinates (both start and end) of proxies for your sentinel variants; along with the corresponding sentinel rsIDs and r2 value across six columns in total with the below column names. Your proxies can be derived either directly from your sample or, for example, from 1000 genomes data. If proxies have not been pre-filtered prior to input, then ProGeMM will filter based on a user-defined r2 threshold (default r2 is 0.8) if desired.
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

source("/data/programs/scripts/SNiPA-annotation-new/snipa_anno_ext.R")
r2.cutoff <- 0.8 #*** ADD parameter with cutoff
rm(tophit)
tophit=regions
#In cases where there is no rsID for a sentinel variant then the notation "chr:start" should be used (i.e., 1:11856378).
tophit <- tophit %>%
    mutate(ID = ifelse(startsWith(ID, "rs"), ID, paste0(CHROM, ":", GENPOS)))



#res <- data.frame()
res <- data.frame(
  CHR = numeric(),
  POS1 = numeric(),
  POS2 = numeric(),
  R2 = numeric(),
  D = numeric(),
  DPRIME = numeric(),
  RSID = character(),
  RSALIAS = character(),
  MINOR = character(),
  MAF = numeric(),
  MAJOR = character(),
  CMMB = numeric(),
  CM = numeric(),
  stringsAsFactors = FALSE)
for (i in 1:nrow(tophit)) {
    chr <- tophit[i, "CHROM"]
    pos <- tophit[i, "GENPOS"]
    cat(paste0("leadSNP: ", tophit[i, "ID"], "\n"))
    ld.df <- ld.basic(snp_chr=chr, snp_pos1=pos, snp_pos2=pos)
    #print(colnames(ld.df))
    if (is.null(ld.df)) {
        print("ld.df is null")
        ld.df <- c(paste0("chr",chr), pos, pos, 1, NA, NA, tophit[i, "ID"], NA,
             tophit[i, "ALLELE0"],  tophit[i, "A1FREQ"],tophit[i, "ALLELE1"],NA,NA,tophit[i, "ID"])
          
        res[nrow(res) + 1,] <- ld.df
        next
    }
    print("nrow ld.df")
    print(nrow(ld.df))
    ld.df$LEAD_rsID <- tophit[i, "ID"]
    res <- rbind(res, ld.df)
}

res <- res[res$R2 > r2.cutoff, ]

#remove "chr" from chromosome name
res$CHR <- sub(pattern="chr", replacement="", x=res$CHR, fixed=TRUE)


res <- res[, c("RSID", "CHR", "POS2", "LEAD_rsID","MINOR","MAJOR" "R2")]

print("Dimension proxies data:")
dim(res) #*** this data also include the indexSNPs! should it be removed?

names(res)[names(res)=="RSID"] <- "PROXY_rsID"
names(res)[names(res)=="CHR"] <- "PROXY_CHR"
names(res)[names(res)=="POS2"] <- "PROXY_START"
names(res)[names(res)=="R2"] <- "r2"
res$PROXY_END <- as.numeric(res$PROXY_START) + nchar(as.character(res$MINOR)) - 1

cat("head proxie file \n")
head(res[c("PROXY_rsID", "PROXY_CHR", "PROXY_START", "PROXY_END","LEAD_rsID", "r2")])
write.table(x=res[c("PROXY_rsID", "PROXY_CHR", "PROXY_START", "PROXY_END",
    "LEAD_rsID", "r2")], file= paste0(opt$output_folder, "/proxies.txt"), quote = FALSE, sep = "\t",
    row.names = FALSE, col.names = TRUE)


res$strand <- "+"
res$allele <- paste(res$MINOR, res$MAJOR, sep = "/")

cat("head vep input file \n")
head(res[c("PROXY_CHR", "PROXY_START", "PROXY_END", "allele",
"strand", "PROXY_rsID")])
#dim(res)
#dim(res[!duplicated(res),])

write.table(x=res[c("PROXY_CHR", "PROXY_START", "PROXY_END", "allele",
    "strand", "PROXY_rsID")], file = paste0(opt$output_folder, "/proxies_vep.txt"), quote = FALSE, sep = "\t",
        row.names = FALSE, col.names = FALSE)


#####################################
# prepare magma input file
#output should look like:
#SNP P N
#rs367896724 0.810602 37941
cat("Prepare magma input file \n")
input_magma=input

input_magma$P=10^(-input_magma$LOG10P)

input_magma=input_magma[,c("ID", "P", "N")]
colnames(input_magma)=c("SNP", "P", "N")
print(head("input magma"))
head(input_magma)
write.table(input_magma,paste0(opt$output_folder,"/input_magma.txt"), row.names=F, col.names=T, sep="\t", quote=F)

cat("Prepprocessing finished \n")
