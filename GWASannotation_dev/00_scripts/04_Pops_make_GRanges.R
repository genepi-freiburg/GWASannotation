prefix <- sub(".preds*", "", basename(PoPS_results))
print(prefix)

rdata.file <- paste0(output_path, "_PoPS.RData")

res <- read.table(PoPS_results, header=TRUE, sep="\t")
print("Pops results (head)")
print(head(res))

load(file.path(gene_model_filename))
# ensembl_genes is loaded, create a data frame for function gene_annotator
# To make it accessible to the function - use <<-
ensembl_genes.df <<- as.data.frame(ensembl_genes)
gene_info <- ensembl_genes.df[unique(res$ENSGID), c("gene_id", "hgnc_symbol",
   "gene_biotype", "seqnames", "start", "end")]

names(gene_info)[names(gene_info) == "gene_id"] <- "ensembl_gene_id"
names(gene_info)[names(gene_info) == "seqnames"] <- "chromosome_name"
gene_info <- data.table(gene_info)


# format the gene_info object slightly
if(length(grep("CHR", gene_info$chromosome_name))) {
    gene_info <- gene_info[-grep("CHR", gene_info$chromosome_name),]
}
gene_info[is.na(gene_info)] <- "-"
gene_info[gene_info == ""] <- "-"


# Remove rows with duplicate Ensembl gene id
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]


# Add PoPS_Score
print("Add PoPS_Score")
gene_info <- merge(gene_info, res[, c("ENSGID", "PoPS_Score")],
    by.x="ensembl_gene_id", by.y="ENSGID")
gene_info$start=as.numeric(gene_info$start)
gene_info$end=as.numeric(gene_info$end)

ranges_object <- GRanges(seqnames = gene_info$chromosome_name,
                         ranges = IRanges(start = gene_info$start,
                                          end = gene_info$end),
                 GRCh38_ensembl_gene_id = gene_info$ensembl_gene_id,
                 GRCh38_hgnc_symbol = gene_info$hgnc_symbol,
                 gene_biotype = gene_info$gene_biotype,
                 PoPS_Score = gene_info$PoPS_Score)

names(ranges_object) <- gene_info$ensembl_gene_id
annotated_PoPS_genes <- list(PoPS=ranges_object)

save(annotated_PoPS_genes, file = rdata.file)
print(rdata.file)
cat("\n## Pops make GRanges finished ##\n")
