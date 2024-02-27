library(readxl)
library(writexl)
library(dplyr)
options(scipen=999)
suppressMessages(library(optparse))

option_list = list(
  make_option("--tophit.file", action="store", default=NA, type='character', help="tophit.file full path [required]"),
  make_option("--progem.file", action="store", default=NA, type='character', help="progem.file summary full path [required]"),
  make_option("--out.file", action="store", default=NA, type='character', help="out.file full path [required]")
)

opt = parse_args(OptionParser(option_list=option_list))


tophit.file <- opt$tophit.file
progem.file <- opt$progem.file
out.file <- opt$out.file

#tophit <- as.data.frame(read_excel(tophit.file))
#tophit=read.table(opt$tophit.file, header=T)
load(opt$tophit.file)
tophit=peaks500$peak_regions_full
dim(tophit)
head(tophit)
tophit<- tophit %>%
    mutate(ID = ifelse(startsWith(ID, "rs"), ID, paste0(CHROM, ":", GENPOS)))


#tophit <- tophit[, c("ID", "CHROM","GENPOS", "ALLELE0", "ALLELE1", "beta", "SE", "pval", "AF_coded_all", "n_total", "FUNC", "GENES", "MarkerName")]

progem <- read.table(progem.file, header=TRUE, sep="\t")
dim(progem)
head(progem)

m <- merge(tophit, progem, by.x="ID", by.y="rsID", sort=FALSE)
dim(m)
names(m)

write_xlsx(x=m, path = out.file, format_headers = TRUE)
print(out.file)
