library(GenomicRanges)
library(biomaRt)
library(data.table)

preds.file <- "UMOD_meta_UKB50k.preds"
rdata.file <- "UMOD_meta_UKB50k_pops.RData"

print(file.exists(preds.file))

res <- read.table(preds.file, header=TRUE, sep="\t")
print(dim(res))

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org",
    path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
gene_info <- data.table(getBM(attributes = c("ensembl_gene_id", "hgnc_symbol",
                                            "gene_biotype", "chromosome_name",
                                            "start_position", "end_position"),
        filters = "ensembl_gene_id", values = res$ENSGID, mart = ensembl))

# format the gene_info object slightly
if (length(grep("CHR", gene_info$chromosome_name))) {
    gene_info <- gene_info[-grep("CHR", gene_info$chromosome_name),]
}
gene_info[is.na(gene_info)] <- "-"
gene_info[gene_info == ""] <- "-"
dim(gene_info)

# Remove rows with duplicate Ensembl gene id
gene_info <- gene_info[!duplicated(gene_info$ensembl_gene_id), ]
dim(gene_info)

# Add PoPS_Score
gene_info <- merge(gene_info, res[, c("ENSGID", "PoPS_Score")],
    by.x="ensembl_gene_id", by.y="ENSGID")
dim(gene_info)

ranges_object <- GRanges(seqnames = gene_info$chromosome_name,
                         ranges = IRanges(start = gene_info$start_position,
                                          end = gene_info$end_position),
                 GRCh38_ensembl_gene_id = gene_info$ensembl_gene_id,
                 GRCh38_hgnc_symbol = gene_info$hgnc_symbol,
                 GRCh37_ensembl_gene_id = gene_info$ensembl_gene_id,
                 GRCh37_hgnc_symbol = gene_info$hgnc_symbol,
                 gene_biotype = gene_info$gene_biotype,
                 PoPS_Score = gene_info$PoPS_Score)

names(ranges_object) <- gene_info$ensembl_gene_id
annotated_metabolic_genes <- list(PoPS=ranges_object)

save(annotated_metabolic_genes, file = rdata.file)
print(rdata.file)
