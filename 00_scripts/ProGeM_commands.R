
########################################   CODE FOR EXECUTION   #######################################

### IMPORT RELEVANT PACKAGES
library(GenomicAlignments)
library(GenomicFeatures)
library(GenomicRanges)
library(biomaRt)
library(data.table)

cat("\nImported required packages.\n")


### READ IN SENTINEL AND PROXY DATA FILES AND FORMAT APPROPRIATELY

# read in sentinel and proxy data files
sentinel_data <- data.table(read.table(file = file.path(sentinel_filename), header = TRUE,
                                       quote = NULL, sep = "\t",
                                       colClasses = c("character", "character", "numeric", "numeric")))
sentinel_data=unique(sentinel_data)
proxy_data <- data.table(read.table(file = file.path(proxy_filename), header = TRUE,
                                    quote = NULL, sep = "\t", stringsAsFactors = FALSE,
                                    colClasses = c("character", "character", "numeric", "numeric",
                                                   "character", "numeric")))

# filter proxies based on LD (if required) and remove sentinel records from proxy_data
if(filtering_required == TRUE) {
  proxy_data <- proxy_data[proxy_data$r2 >= r2_thresh,]
}
proxy_data <- proxy_data[-which(proxy_data$PROXY_rsID == proxy_data$LEAD_rsID),]


# add a new column containing chromsome and position information in the format: "chr_position_" (this 
# will be a search term for GTEX cis-eQTL data further down the script)
GTEx_search_term <- paste(sentinel_data$CHR, "_", sentinel_data$START, "_", sep = "")
sentinel_data <- cbind(sentinel_data, GTEx_search_term)
sentinel_data$GTEx_search_term <- as.character(sentinel_data$GTEx_search_term)

GTEx_search_term <- paste(proxy_data$PROXY_CHR, "_", proxy_data$PROXY_START, "_", sep = "")
proxy_data <- cbind(proxy_data, GTEx_search_term)
proxy_data$GTEx_search_term <- as.character(proxy_data$GTEx_search_term)

cat("Sentinel and proxy information read and formatted.\n\n")


### PULL OUT ALL GENES THAT RESIDE WITHIN +/-"X"kb FROM EACH SENTINEL VARIANT

# create a genomic ranges object for the sentinel variants
sentinel_ranges <- GR_object_creator(sentinel_data$CHR, sentinel_data$START, sentinel_data$END,
                                     sentinel_data$rsID)

# load a TxDb object containing ensembl genes according to build GRCh38 - this is based on the following
load(file.path(gene_model_filename))

# identify all genes that reside +/- X kb from the sentinel variants
local_genes <- find_overlapping_genes(sentinel_ranges, ensembl_genes, interval_kb)
# determine distance between each local gene and the corresponding sentinel variant - rank accordingly
local_genes <- gene_variant_distance_finder(sentinel_ranges, local_genes)

# annotate the genes in the above local_genes object
local_genes_annotation <- gene_annotator(unique(names(local_genes)))


### IDENTIFY ALL BOTTOM-UP CANDIDATE GENES

cat(paste("-------------------------------------------------------------------\n",
          "Beginning bottom-up annotation...\n\n", sep = ""))

## pull out all genes from "GRCh37_genes" that overlap the LD region around each sentinel variant

# create a genomic ranges object containing LD regions
LD_region_ranges <- LD_region_range_finder(sentinel_data, proxy_data)

# identify all genes that overlap the LD regions (+/-Xkb)
LD_region_genes <- find_overlapping_genes(LD_region_ranges, ensembl_genes, 
                                          LD_region_overhang_kb)

# annotate the genes in the above LD_region_overlapping_genes object
LD_region_genes_annotated <- gene_annotator(unique(names(LD_region_genes)))

cat("\t- Genes overlapping an LD region have been identified.\n")


## pull out the nearest genes to each sentinel variant

nearest_genes <- nearest_genes_selector(local_genes, local_genes_annotation, biotype_of_interest, 
                                        number_of_nearest)
print("nearest_genes")
print(head(nearest_genes))

cat(paste("\t- The ", number_of_nearest, " genes nearest to each sentinel have been identified.\n",
          sep = ""))


## identify all cis-eQTL targets of the sentinel and proxy variants from GTEx data v6p

cat("\t- Searching for cis-eQTL targets...\n")

# search cis-eQTL files (NOTE: will take approx 10-15 minutes if all tissues have been chosen)
cis_eQTL_hits <- cis_eQTL_target_finder(eQTLdata_dir, tissues_of_interest, sentinel_data, proxy_data)
print(cis_eQTL_hits)

# format the cis_eQTL_hits list; i.e., add rsids and remove .* from ensembl_ids
cis_eQTL_hits <- cis_eQTL_formatter(sentinel_data, proxy_data, cis_eQTL_hits)

cis_eQTL_targets <- cis_eQTL_gene_extractor(sentinel_data, proxy_data, cis_eQTL_hits)

cis_eQTL_target_annotation <- gene_annotator(unique(cis_eQTL_targets$ensembl_id))
cat("\t- Cis-eQTL targets have been identified.\n")


### SUMARISE BOTTOM-UP INFORMATION
bottom_up_summary <- bottom_up_summariser(sentinel_data, nearest_genes, LD_region_genes, LD_region_genes_annotated, cis_eQTL_targets,
                                          cis_eQTL_target_annotation)
print("bottom_up_summariser done")
# filter by biotype_of_interest
bottom_up_summary <- bottom_up_summary[bottom_up_summary$gene_biotype == biotype_of_interest,]

print("IDENTIFY GENES WITH INCREASED LIKELIHOOD OF HAVING A FUNCTIONAL IMPACT USING VEP ANNOTATION")
### IDENTIFY GENES WITH INCREASED LIKELIHOOD OF HAVING A FUNCTIONAL IMPACT USING VEP ANNOTATION
# note: to save resources with particularly large datasets it may be worthwhile filtering for only coding variants when running VEP,
# however, bear in mind that this would limit your variant consequences of interest; i.e., the moderate impact consequence "regulatory
# region ablation" would be lost

# read in VEP output file
VEP_annotations <- data.table(read.table(file = file.path(VEP_filename), header = FALSE,
                                         quote = NULL, sep = "\t", stringsAsFactors = FALSE))
#print(head(VEP_annotations))
print(dim(VEP_annotations))
if(filtering_required == TRUE & length(r2_column) > 0) {
  VEP_annotations <- VEP_annotations[VEP_annotations[[r2_column]] >= r2_thresh,]
}

# filter out all records that do not have either a "HIGH" or "MODERATE" IMPACT
VEP_IMPACTS_of_interest <- c("HIGH", "MODERATE")
dim(VEP_annotations)
VEP_annotations <- VEP_annotations[VEP_annotations[[IMPACT_column]] %in% VEP_IMPACTS_of_interest,]
dim(VEP_annotations)
# identify leads and proxies most likely to have a functional impact as per VEP annotation
SNP_IMPACT_annotations <- VEP_integrator(sentinel_data, bottom_up_summary, VEP_annotations, sentinel_rsID_column, proxy_rsID_column,
                                         ensembl_gene_id_column, IMPACT_column)

# combine bottom_up_summary with the SNP_IMPACT_annotations
bottom_up_summary <- cbind(bottom_up_summary, SNP_IMPACT_annotations)

cat("\t- VEP annotations have been integrated.\n\n")

# COLOC data
# Coloc with eQTL

  # read in coloc eQTL file
coloc_eqtl_data <- data.table(read.table(file =
                                         file.path(COLOC_EQTL_filename),
                                       header = TRUE, quote = NULL, sep = "\t", stringsAsFactors = FALSE))

if (length(eqtl_tissue_col) > 0) {
coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[eqtl_tissue_col]] %in% eqtl_tissues_of_interest, ]
}

if (length(eqtl_sumstats_2_max_nlog10P_col) > 0) {
coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[eqtl_sumstats_2_max_nlog10P_col]] > eqtl_sumstats_2_max_nlog10P_thresh, ]
}

if (length(eqtl_PP.H4.abf_col) > 0) {
coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[eqtl_PP.H4.abf_col]] > eqtl_PP.H4.abf_thresh, ]
}

if (length(coloc_eqtl_gene_type_col) > 0) {
coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[coloc_eqtl_gene_type_col]] %in% biotype_of_interest, ]
}

# filter out all records that do not have a "cis"
if (length(eqtl_cis_trans_col)) {
coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[eqtl_cis_trans_col]] %in% eqtl_cis_trans_sel, ]
}

coloc_eqtl_data <- as.data.frame(coloc_eqtl_data)
coloc_eqtl_data <- coloc_eqtl_data[!duplicated(coloc_eqtl_data[,
                                                             c(coloc_eqtl_sentinel_rsID_col, coloc_eqtl_ensembl_gene_id_col)]), ]
coloc_eqtl_data <- data.table(coloc_eqtl_data)

# formatting so that it can be added to the bottom-up summary
COLOC_EQTL_annotations <- COLOC_integrator(sentinel_data, bottom_up_summary,
                                         coloc_eqtl_data, coloc_eqtl_sentinel_rsID_col,
                                         coloc_eqtl_ensembl_gene_id_col)

# combine bottom_up_summary with the COLOC EQTL annotation
bottom_up_summary <- cbind(bottom_up_summary, COLOC_EQTL_annotations)


cat("\t- Coloc eQTL information have been integrated.\n")

  # read in coloc pQTL file
coloc_pqtl_data <- data.table(read.table(file =
                                         file.path(COLOC_PQTL_filename),
                                       header = TRUE, quote = NULL, sep = "\t", stringsAsFactors = FALSE))

if (length(pqtl_sumstats_2_max_nlog10P_col) > 0) {
coloc_pqtl_data <- coloc_pqtl_data[coloc_pqtl_data[[pqtl_sumstats_2_max_nlog10P_col]] > eqtl_sumstats_2_max_nlog10P_thresh, ]
}

if (length(pqtl_PP.H4.abf_col) > 0) {
coloc_pqtl_data <- coloc_pqtl_data[coloc_pqtl_data[[pqtl_PP.H4.abf_col]] > pqtl_PP.H4.abf_thresh, ]
}

# filter out all records that do not have a "cis"
if (length(pqtl_cis_trans_col)) {
coloc_pqtl_data <- coloc_pqtl_data[coloc_pqtl_data[[pqtl_cis_trans_col]] %in% pqtl_cis_trans_sel, ]
}

coloc_pqtl_data <- as.data.frame(coloc_pqtl_data)
coloc_pqtl_data <- coloc_pqtl_data[!duplicated(coloc_pqtl_data[,
                                                             c(coloc_pqtl_sentinel_rsID_col, coloc_pqtl_ensembl_gene_id_col)]), ]
coloc_pqtl_data <- data.table(coloc_pqtl_data)

# formatting so that it can be added to the bottom-up summary
COLOC_PQTL_annotations <- COLOC_integrator(sentinel_data, bottom_up_summary,
                                         coloc_pqtl_data, coloc_pqtl_sentinel_rsID_col,
                                         coloc_pqtl_ensembl_gene_id_col, QTL_type = "pQTL")

# combine bottom_up_summary with the COLOC EQTL annotation
bottom_up_summary <- cbind(bottom_up_summary, COLOC_PQTL_annotations)


cat("\t- Coloc pQTL information have been integrated.\n\n")

# classify bottom_up evidence according to source
# 1 = gene contains a lead or proxy with either a HIGH or MODERATE IMPACT (proximity and cis-eQTL evidence may also be present)
# 2 = gene is proximal to lead; i.e., nearest or LD overlapping (cis-eQTL evidence may also be present)
# 3 = gene is a cis-eQTL target of either the lead or a proxy variant
bottom_up_summary <- bottom_up_classifier(bottom_up_summary)

cat("Bottom-up annotation complete!\n\n")


### IDENTIFY ALL TOP-DOWN CANDIDATES

cat(paste("-------------------------------------------------------------------\n",
          "Beginning top-down annotation...\n\n", sep = ""))

# load in PoPS result list of genomicRanges objects
load(file.path(filename_PoPS))

# identify top-down candidates and score each candidate accordingly
PoPS_overlapping_genes <- top_down_candidate_identifier(sentinel_ranges, annotated_PoPS_genes, 
                                                             interval_kb)

write.table(PoPS_overlapping_genes, file = file.path(output_dir, "PoPS_overlapping_genes.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
#print("PoPS_overlapping_genes")
#print(head(PoPS_overlapping_genes))

top_down_scores <- top_down_scorer(sentinel_data, PoPS_overlapping_genes)
#print("top_down_scores")
#print(head(top_down_scores))

cat("\t- Top-down candidates have been identified and scored.\n\n")


# filter by biotype of interest
top_down_scores <- top_down_scores[top_down_scores$gene_biotype == biotype_of_interest,]

cat("Top-down annotation complete.\n\n")


### IDENTIFY CO-OCCURRING GENES
cooccurring_candidates <- converging_candidate_finder(sentinel_data, bottom_up_summary, top_down_scores)

cat(paste("-------------------------------------------------------------------\n",
          "Cooccurring candidates identified.\n", sep = ""))


### WRITE OUT PRIMARY INFORMATION TO FILE

# bottom-up summary
write.table(bottom_up_summary, file = file.path(output_dir, "OUTPUT_bottom_up_summary.txt"), 
            quote = FALSE, row.names = FALSE, sep = "\t")

# top-down scores
write.table(top_down_scores, file = file.path(output_dir, "OUTPUT_top_down_scores.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# cooccurring candidates
write.table(cooccurring_candidates, file = file.path(output_dir, "OUTPUT_cooccurring_candidates.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# HIGH and MODERATE VEP dataframe
write.table(VEP_annotations, file = file.path(output_dir, "OUTPUT_VEP_IMPACT_annotations.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# lead eQTLs
write.table(cis_eQTL_hits[[1]], file = file.path(output_dir, "OUTPUT_lead_cis_eQTLs.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# proxy eQTLs
write.table(cis_eQTL_hits[[2]], file = file.path(output_dir, "OUTPUT_proxy_cis_eQTLs.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# nearest genes
write.table(nearest_genes, file = file.path(output_dir, "OUTPUT_nearest_genes.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

cat(paste0("Output files written to", output_dir, "\n"))
cat("Annotation complete!!\n\n")

#######################################################################################################
