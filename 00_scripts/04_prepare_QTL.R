#ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",
#                  dataset="hsapiens_gene_ensembl")


input_path=output_path
tophit.file=paste0(input_path, "_lead.txt")
regions.file=paste0(input_path, "_coloc_regions.RDS")

#inc.datasets <- c("GTEXv8", "Kidney_eQTL","ARIC_pGWAS", "Icelanders_pGWAS", "UKB_PPP_EUR")
inc.datasets <- datasets_coloc
print(inc.datasets)
merge_coloc_lead <- function(coloc_path, dataset, regions.file, tophit.file) {
    summary.file <- paste0(coloc_path, dataset, "_annot_unfilt.RDS")
    # Check if files exist
    if (!file.exists(summary.file) | !file.exists(regions.file) | !file.exists(tophit.file)) {
        stop("One or more files do not exist.")
    }
    
    # Reading necessary files
    summary <- readRDS(summary.file)
    regions <- readRDS(regions.file)
    regions <- regions[["coloc_regions"]]
    regions <- regions[regions$comment == "PASS", c("CHR_var", "BP_START_var", "BP_STOP_var", "rsID")]
    tophit <- read.table(tophit.file, header = TRUE)
    tophit <- tophit[, c("rsID", "CHR", "START")]
    colnames(tophit) <- c("rsID", "CHR_var", "snpPOS")
    
    summary <- summary[!is.na(summary$PP.H4.abf) & summary$PP.H4.abf > eqtl_PP.H4.abf_thresh, ]
    cat("nr of rows with PP.H4 > ", eqtl_PP.H4.abf_thresh," : ", nrow(summary), "\n")
    
    # Merge coloc results with regions to get rsID
    summary.m1 <- merge(summary, regions, sort = FALSE, by = c("CHR_var", "BP_START_var", "BP_STOP_var"))
    # Print(dim(summary.m1))
    
    # Merge with lead file
    summary.m2 <- merge(tophit, summary.m1, sort = FALSE)
    cat("nr rows merged file: ", nrow(summary.m2), "\n")
    return(summary.m2)
}

process_dataset_eQTL <- function(coloc_path, dataset, regions.file, tophit.file) {
    inc.cols <- c("rsID","trait", "Name", "sumstats_2_max_nlog10P", "PP.H4.abf")
    inc.cols2 <- c("_Tissue", "_gene_id", "_gene_type", "_gene_name", "_cis_trans")
    print(dataset)
    
    data <- merge_coloc_lead(coloc_path, dataset, regions.file, tophit.file)
    if(nrow(data>0)){
        data$trait=dataset
        #data$rsID <- paste(data$chr, data$position, sep = ":")
        data$Name <- paste(data$CHR_var, data$snpPOS, sep = ":")
        print(all(inc.cols %in% names(data)))
        print(all(paste0(dataset, inc.cols2) %in% names(data)))
        cols <- c(inc.cols, paste0(dataset, inc.cols2))
        data <- data[, cols]
        print(dim(data))
        names(data) <- sub(pattern=paste0(dataset, "_"), replacement="", x=names(data))
        data$gene_id <- sub(pattern="\\.[0-9]*$", replacement="", x=data$gene_id)
        data <- data[data$cis_trans =="cis", ]
        data <- data[complete.cases(data$gene_id), ]
        print(dim(data))
        # Modify "_Tissue" column for eQTLGen dataset
        if (dataset == "eQTLGen_Tissue") {
            data$Tissue <- "eQTLGen_blood"
        }
        return(data)
    }
}

process_dataset_pQTL <- function(coloc_path, dataset, regions.file, tophit.file) {
    inc.cols <- c("rsID", "trait", "Name", "sumstats_2_max_nlog10P", "PP.H4.abf")
    inc.cols2 <- c("_cis_trans")
    data <- merge_coloc_lead(coloc_path, dataset, regions.file, tophit.file)
    if(nrow(data>0)){
        data$trait=dataset
        data$Name <- paste(data$CHR_var, data$snpPOS, sep = ":")
        print(all(inc.cols %in% names(data)))
        if (dataset == "ARIC_pGWAS") {
            #gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
            #                   filters = "hgnc_symbol",
            #                   values = unique(data$ARIC_pGWAS_entrezgenesymbol), mart = ensembl)
            gene_info <- ensembl_genes.df[unique(data$ARIC_pGWAS_entrezgenesymbol), c("gene_id", "hgnc_symbol")]
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
        names(data) <- sub(pattern=paste0(dataset, "_"), replacement="", x=names(data))
        data <- data[data$cis_trans =="cis", ]
        data <- data[complete.cases(data$gene_id), ]
        print(dim(data))
        return(data)
    }
}

# eQTL
inc.datasets_eQTL <- eQTL_datasets_coloc
out.txt_eQTL <- paste0(coloc_path, "input_Anno_eQTL.txt")
out.df_eQTL <- data.frame()

for (dataset in inc.datasets_eQTL) {
    dataset_data <- process_dataset_eQTL(coloc_path, dataset,  regions.file, tophit.file)
    out.df_eQTL <- rbind(out.df_eQTL, dataset_data)
}

print(dim(out.df_eQTL))
write.table(x=out.df_eQTL, file = out.txt_eQTL, quote=FALSE, sep = "\t", row.names=FALSE)
print(out.txt_eQTL)

# pQTL
inc.datasets_pQTL <- pQTL_datasets_coloc
out.txt_pQTL <- paste0(coloc_path, "input_Anno_pQTL.txt")
out.df_pQTL <- data.frame()

for (dataset in inc.datasets_pQTL) {
    dataset_data <- process_dataset_pQTL(coloc_path, dataset, regions.file, tophit.file)
    out.df_pQTL <- rbind(out.df_pQTL, dataset_data)
}

print(dim(out.df_pQTL))
write.table(x=out.df_pQTL, file = out.txt_pQTL, quote=FALSE, sep = "\t", row.names=FALSE)
print(out.txt_pQTL)
