###############################################################
# THIS IS NOT PART OF GWAS ANNOTATION PIPELINE You need to create the input file for GWASAnnotation
# You need to create the input file for GWASAnnotation before running the pipeline. Use this code as help only
##################################################################

suppressMessages(library(readxl))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
library(Rmpfr)
library(parallel)
devtools::load_all("/data/programs/pipelines/genepicoloc/genepicoloc_package")

#Sahar's example
#Summary statistics: /data/meta_analyses/00_CKDGen/01_analyses/09_GWAS-r5/urate-Sahar/urate-multi/urate_20250103.multi.tbl.gz
#LD file: /data/studies/06_UKBB/02_Projects/19_LD_RefPanel/bed_ABC/ALL/LD_Panel_ALL_combined
#Genome build: GRCh38

#####################################
# load input file
sumstats <- data.table::fread("/data/meta_analyses/00_CKDGen/01_analyses/09_GWAS-r5/urate-Sahar/urate-multi/urate_20250103.multi.tbl.gz")

#using only chr 4 region 4  85527498    89566483 for testing
sumstats_chr4 <- sumstats[which(sumstats$Chr == 4 & sumstats$Pos_b38 > 85527498 & sumstats$Pos_b38 < 89566483),]
head(sumstats)
rm(sumstats)

sumstats_chr4$nlog10 <- -log10(as.numeric(sumstats_chr4$`P-value`))
#Recalculate -log10 for small p-values (< 1e-300) using mpfr in parallel
threshold <- 1e-300
small_pvals <- which(as.numeric(sumstats_chr4$`P-value`) < threshold)
length(small_pvals)

num_cores <- detectCores() - 1
# Perform mpfr-based -log10 calculation for small p-values
mpfr_results <- unlist(mclapply(sumstats_chr4$`P-value`[small_pvals], function(x) {
  p_val_mpfr <- mpfr(x, precBits = 300)
  log_val <- -log10(p_val_mpfr)
  as.numeric(log_val)  # Convert mpfr result to numeric!!!!
}, mc.cores = num_cores))
sumstats_chr4$nlog10[small_pvals] <- mpfr_results
#sumstats_chr4$nlog10 <- round(sumstats_chr4$nlog10, 3) #Round results

sumstats_chr4[which.max(sumstats_chr4$nlog10), ]

sumstats_chr4$pval=as.numeric(sumstats_chr4$pval)
sumstats_chr4$nlog10_p= -log10(sumstats_chr4$pval)

sumstats_1 <- read_sumstats(sumstats_chr4,
                            Name="MarkerName",
                            rsID = "MarkerName",
                            CHR = "Chr",
                            POS = "Pos_b38",
                            A1 = "Allele1", #re-check
                            A2 = "Allele2", #re-check
                            BETA = "Effect",
                            SE = "StdErr",
                            nlog10p_value = "nlog10",
                            AF = "Freq1")
sumstats_1$N=sumstats_chr4$n_total_sum_sum #N is needed for magma input
colnames(sumstats_1)[2] <- "rsID"

cat("number of SNPs ", nrow(sumstats_1))

unique_rows <- sumstats_1 %>%
  group_by(Name) %>%
  slice(which.max(nlog10P)) %>%
  ungroup()
cat("number of SNPs after removing duplicated 'Name'", nrow(unique_rows))

unique_rows$Name <- paste0("chr", unique_rows$Name)

saveRDS(sumstats_1, "/data/programs/pipelines/GWASannotation/02_output/test_Sahar_2025/test_chr4_hg38.RDS")
saveRDS(unique_rows, "/data/programs/pipelines/GWASannotation/02_output/test_Sahar_2025/test_chr4_hg38_dedup.RDS")


#Check if overlaps with coloc datasets
c <- data.table::fread("/data/public_resources/CKDGen/preprocessing/Wuttke2019/CKD_overall_ALL_JW_20180223_nstud30.dbgap.txt_hg38.gz")
length(unique(unique_rows$Name))
#21283
#Number of sumstats SNPs in coloc
int <- unique_rows[which(unique_rows$Name %in% c$Name_hg38),]
length(unique(int$Name))

#Number of coloc SNPs in sumstats (sanity check)
nrow(c[which(c$Name_hg38 %in% unique_rows$Name),])


split_data <- do.call(rbind, strsplit(unique_rows$Name, ":"))
unique_rows$A1 <- split_data[, 3]
unique_rows$A2 <- split_data[, 4]
unique_rows$Name_inv <- paste0("chr",unique_rows$CHR,":", unique_rows$POS,":", unique_rows$A2,":", unique_rows$A1)
int_inv <- unique_rows[which(unique_rows$Name_inv %in% c$Name_hg38),]
length(unique(int_inv$Name))
#6469
nrow(c[which(c$Name_hg38 %in% unique_rows$Name_inv),])
#6469

unique_rows$Name <- unique_rows$Name_inv
unique_rows <- unique_rows[, !colnames(unique_rows) %in% "Name_inv"]

saveRDS(unique_rows, "/data/programs/pipelines/GWASannotation/02_output/test_Sahar_2025/test_chr4_hg38_dedup_invAlelle_to_coloc.RDS")
