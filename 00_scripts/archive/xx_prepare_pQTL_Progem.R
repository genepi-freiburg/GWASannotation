library(readxl)
library(biomaRt)

inc.datasets <- c("ARIC_pGWAS", "Icelanders_pGWAS", "UKB_PPP_EUR")
in.dir <- "summary_individual_dataset"
out.txt <- "input_Progem_pQTL.txt"
inc.cols <- c("rsID", "trait", "sumstats_2_min_P", "PP.H4.abf")

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl")

inc.cols2 <- c("_cis_trans")
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
    # print(all(paste0(dataset, inc.cols2) %in% names(data)))
    print(names(data))

    if (dataset == "ARIC_pGWAS") {
        gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
            filters = "hgnc_symbol",
            values = unique(data$ARIC_pGWAS_entrezgenesymbol), mart = ensembl)
        print(gene_info)
        data <- merge(data, gene_info, by.x="ARIC_pGWAS_entrezgenesymbol",
            by.y="hgnc_symbol")

        names(data)[names(data) == "ensembl_gene_id"] <- "gene_id"

    } else if (dataset == "Icelanders_pGWAS") {
        names(data)[names(data) == "Icelanders_pGWAS_Ensembl.Gene.ID"] <-
            "gene_id"
    } else if (dataset == "UKB_PPP_EUR") {
        names(data)[names(data) == "UKB_PPP_EUR_ensembl_id"] <- "gene_id"
    }

    cols <- c(inc.cols, "gene_id", paste0(dataset, inc.cols2))
    data <- data[, cols]
    print(dim(data))

    names(data) <- sub(pattern=dataset2, replacement="", x=names(data))
    # data$gene_id <- sub(pattern="\\.[0-9]*$", replacement="", x=data$gene_id)

    data <- data[data$cis_trans =="cis", ]
    data <- data[complete.cases(data$gene_id), ]

    print(dim(data))
    out.df <- rbind(out.df, data)
}

print(dim(out.df))
write.table(x=out.df, file = out.txt, quote=FALSE, sep = "\t", row.names=FALSE)
print(out.txt)
