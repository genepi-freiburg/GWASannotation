library(readxl)

inc.datasets <- c("GTEXv8", "Kidney_eQTL")
in.dir <- "summary_individual_dataset"
out.txt <- "input_Progem_eQTL.txt"
inc.cols <- c("rsID", "trait", "sumstats_2_min_P", "PP.H4.abf")
inc.cols2 <- c("_Tissue", "_gene_id", "_gene_type", "_gene_name", "_cis_trans")

out.df <- data.frame()

for (dataset in inc.datasets) {
    print(dataset)
    dataset2 <- paste0(dataset, "_")

    infile <- file.path(in.dir, paste0("summary_coloc_", dataset, ".xlsx"))
    print(file.exists(infile))

    data <- as.data.frame(read_excel(infile))
    print(dim(data))
    data$rsID <- paste(data$chr, data$position, sep = ":")

    print(all(inc.cols %in% names(data)))
    print(all(paste0(dataset, inc.cols2) %in% names(data)))

    cols <- c(inc.cols, paste0(dataset, inc.cols2))
    data <- data[, cols]
    print(dim(data))

    names(data) <- sub(pattern=dataset2, replacement="", x=names(data))
    data$gene_id <- sub(pattern="\\.[0-9]*$", replacement="", x=data$gene_id)
    data <- data[data$cis_trans =="cis", ]
    data <- data[complete.cases(data$gene_id), ]
    # df2[complete.cases(df2[c(1,3)]),]

    print(dim(data))

    out.df <- rbind(out.df, data)
}

print(dim(out.df))
write.table(x=out.df, file = out.txt, quote=FALSE, sep = "\t", row.names=FALSE)
print(out.txt)
