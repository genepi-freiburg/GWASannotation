library(readxl)
library(biomaRt)

library(readxl)
library(writexl)

inc.datasets <- c("GTEXv8", "Kidney_eQTL","ARIC_pGWAS", "Icelanders_pGWAS", "UKB_PPP_EUR")

merge_tophit <- function(dataset, tophit) {
    summary.file <- paste0(dataset, ".RDS")
    print(file.exists(summary.file))

    summary <- readRDS(summary.file)

    summary <- summary[!is.na(summary$PP.H4.abf) & summary$PP.H4.abf> 0.5, ]
    print(dim(summary))

    tophit=read.table("/data/programs/pipelines/GWASannotation/02_output/test/test_sentinel.txt", header=T)
    summary.m2 <- merge(tophit, summary, by.x=c("CHR_var", "START"), by.y=c("CHR", "START"),  sort = FALSE)
    print(dim(summary.m2))

    outfile <- paste0("summary_coloc_", dataset, ".xlsx")
    write_xlsx(x = summary.m2, path = file.path(out.dir, outfile))
    print(outfile)
}


process_dataset_eQTL <- function(dataset) {
    inc.cols <- c("rsID", "trait", "sumstats_2_min_P", "PP.H4.abf")
    inc.cols2 <- c("_Tissue", "_gene_id", "_gene_type", "_gene_name", "_cis_trans")
    print(dataset)
    dataset2 <- paste0(dataset, ".RDS")
    print(file.exists(infile))
    data <- as.data.frame(read_excel(in.file))
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
    print(dim(data))
    return(data)
}

process_dataset_pQTL <- function(dataset) {
    inc.cols <- c("rsID", "trait", "sumstats_2_min_P", "PP.H4.abf")
    inc.cols2 <- c("_cis_trans")
    print(dataset)
    dataset2 <- paste0(dataset, ".RDS")
    infile <- file.path(in.file, paste0("summary_coloc_", dataset, ".xlsx"))
    print(file.exists(infile))
    data <- as.data.frame(read_excel(infile))
    print(dim(data))
    data$rsID <- paste(data$chr, data$position, sep = ":")
    print(all(inc.cols %in% names(data)))
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
        names(data)[names(data) == "Icelanders_pGWAS_Ensembl.Gene.ID"] <- "gene_id"
    } else if (dataset == "UKB_PPP_EUR") {
        names(data)[names(data) == "UKB_PPP_EUR_ensembl_id"] <- "gene_id"
    }

    cols <- c(inc.cols, "gene_id", paste0(dataset, inc.cols2))
    data <- data[, cols]
    print(dim(data))
    names(data) <- sub(pattern=dataset2, replacement="", x=names(data))
    data <- data[data$cis_trans =="cis", ]
    data <- data[complete.cases(data$gene_id), ]
    print(dim(data))
    return(data)
}

# eQTL
inc.datasets_eQTL <- c("GTEXv8", "Kidney_eQTL")
#in.file_eQTL <- "summary_individual_dataset"
out.txt_eQTL <- "input_Progem_eQTL.txt"
out.df_eQTL <- data.frame()

for (dataset in inc.datasets_eQTL) {
    dataset_data <- process_dataset_eQTL(dataset)
    out.df_eQTL <- rbind(out.df_eQTL, dataset_data)
}

print(dim(out.df_eQTL))
write.table(x=out.df_eQTL, file = out.txt_eQTL, quote=FALSE, sep = "\t", row.names=FALSE)
print(out.txt_eQTL)

# pQTL
inc.datasets_pQTL <- c("ARIC_pGWAS", "Icelanders_pGWAS", "UKB_PPP_EUR")
#in.file_pQTL <- "summary_individual_dataset"
out.txt_pQTL <- "input_Progem_pQTL.txt"
out.df_pQTL <- data.frame()

for (dataset in inc.datasets_pQTL) {
    dataset_data <- process_dataset_pQTL(dataset)
    out.df_pQTL <- rbind(out.df_pQTL, dataset_data)
}

print(dim(out.df_pQTL))
write.table(x=out.df_pQTL, file = out.txt_pQTL, quote=FALSE, sep = "\t", row.names=FALSE)
print(out.txt_pQTL)
