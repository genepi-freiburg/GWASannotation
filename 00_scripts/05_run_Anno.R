
########################################   CODE FOR EXECUTION   #######################################

### READ IN LEAD AND PROXY DATA FILES AND FORMAT APPROPRIATELY

# read in lead and proxy data files
lead_data <- data.table(read.table(file = file.path(lead_filename), header = TRUE,
                                       quote = NULL, sep = "\t",
                                       colClasses = c("character", "character", "numeric", "numeric")))
lead_data=unique(lead_data)
proxy_data <- data.table(read.table(file = file.path(proxy_filename), header = TRUE,
                                    quote = NULL, sep = "\t", stringsAsFactors = FALSE,
                                    colClasses = c("character", "character", "numeric", "numeric",
                                                   "character", "numeric")))


#remove the lead from the proxy_data dataset
proxy_data <- proxy_data[-which(proxy_data$PROXY_rsID == proxy_data$LEAD_rsID),]


# add a new column containing chromsome and position information in the format: "chr_position_" (this 
# will be a search term for GTEX cis-eQTL data further down the script)
GTEx_search_term <- paste(lead_data$CHR, "_", lead_data$START, "_", sep = "")
lead_data <- cbind(lead_data, GTEx_search_term)
lead_data$GTEx_search_term <- as.character(lead_data$GTEx_search_term)

GTEx_search_term <- paste(proxy_data$PROXY_CHR, "_", proxy_data$PROXY_START, "_", sep = "")
proxy_data <- cbind(proxy_data, GTEx_search_term)
proxy_data$GTEx_search_term <- as.character(proxy_data$GTEx_search_term)

cat("lead and proxy information read and formatted.\n\n")


### PULL OUT ALL GENES THAT RESIDE WITHIN +/-"X"kb FROM EACH lead VARIANT

# create a genomic ranges object for the lead variants
lead_ranges <- GR_object_creator(lead_data$CHR, lead_data$START, lead_data$END,
                                     lead_data$rsID)

# load a TxDb object containing ensembl genes according to build GRCh38 - this is based on the following
load(file.path(gene_model_filename))
# ensembl_genes is loaded, create a data frame for function gene_annotator
# To make it accessible to the function - use <<-
ensembl_genes.df <<- as.data.frame(ensembl_genes)

# identify all genes that reside +/- X kb from the lead variants
local_genes <- find_overlapping_genes(lead_ranges, ensembl_genes, interval_kb)
# determine distance between each local gene and the corresponding lead variant - rank accordingly
local_genes <- gene_variant_distance_finder(lead_ranges, local_genes)

# annotate the genes in the above local_genes object
local_genes_annotation <- gene_annotator(unique(names(local_genes)))


### IDENTIFY ALL BOTTOM-UP CANDIDATE GENES

cat(paste("-------------------------------------------------------------------\n",
          "Beginning bottom-up annotation...\n\n", sep = ""))

## pull out all genes from "GRCh37_genes" that overlap the LD region around each lead variant

# create a genomic ranges object containing LD regions
LD_region_ranges <- LD_region_range_finder(lead_data, proxy_data)
save(LD_region_ranges,file = file.path(output_dir, "LD_region_ranges.RData"))
# identify all genes that overlap the LD regions (+/-Xkb)
LD_region_genes <- find_overlapping_genes(LD_region_ranges, ensembl_genes, 
                                          LD_region_overhang_kb)
save(LD_region_genes,file = file.path(output_dir, "LD_region_genes.RData"))
# annotate the genes in the above LD_region_overlapping_genes object
LD_region_genes_annotated <- gene_annotator(unique(names(LD_region_genes)))
save(LD_region_genes_annotated,file = file.path(output_dir, "LD_region_genes_annotated.RData"))
cat("\t- Genes overlapping an LD region have been identified.\n")


## pull out the nearest genes to each lead variant
#print("local_genes")
#print(head(local_genes))
#print("local_genes_annotation")
#print(head(local_genes_annotation))
#print("biotype_of_interest")
#print(biotype_of_interest)
#print("number_of_nearest")
#print(number_of_nearest)

nearest_genes <- nearest_genes_selector(local_genes, local_genes_annotation, biotype_of_interest, 
                                        number_of_nearest)
#print("nearest_genes")
#print(head(nearest_genes))

cat(paste("\t- The ", number_of_nearest, " genes nearest to each lead have been identified.\n",
          sep = ""))


## identify all cis-eQTL targets of the lead and proxy variants from GTEx data v6p

cat("\t- Searching for cis-eQTL targets...\n")

# search cis-eQTL files (NOTE: will take approx 10-15 minutes if all tissues have been chosen)
cis_eQTL_hits <- cis_eQTL_target_finder(eQTLdata_dir, tissues_eQTL_associations, lead_data, proxy_data)
print(cis_eQTL_hits)

# format the cis_eQTL_hits list; i.e., add rsids and remove .* from ensembl_ids
cis_eQTL_hits <- cis_eQTL_formatter(lead_data, proxy_data, cis_eQTL_hits)

cis_eQTL_targets <- cis_eQTL_gene_extractor(lead_data, proxy_data, cis_eQTL_hits)

cis_eQTL_target_annotation <- gene_annotator(unique(cis_eQTL_targets$ensembl_id))
cat("\t- Cis-eQTL targets have been identified.\n")


### SUMARISE BOTTOM-UP INFORMATION
anno_summary <- anno_summariser(lead_data, nearest_genes, LD_region_genes, LD_region_genes_annotated, cis_eQTL_targets,
                                          cis_eQTL_target_annotation)
print("anno_summariser done")
# filter by biotype_of_interest
anno_summary <- anno_summary[anno_summary$gene_biotype == biotype_of_interest,]

print("IDENTIFY GENES WITH INCREASED LIKELIHOOD OF HAVING A FUNCTIONAL IMPACT USING VEP ANNOTATION")
### IDENTIFY GENES WITH INCREASED LIKELIHOOD OF HAVING A FUNCTIONAL IMPACT USING VEP ANNOTATION
# note: to save resources with particularly large datasets it may be worthwhile filtering for only coding variants when running VEP,
# however, bear in mind that this would limit your variant consequences of interest; i.e., the moderate impact consequence "regulatory
# region ablation" would be lost

# read in VEP output file
file_conn <- file(VEP_filename, "r")



# Read the first non-commented line
first_line <- NULL
while (TRUE) {
  line <- readLines(file_conn, n = 1)
  if (length(line) == 0) break
  if (!startsWith(line, "#")) {
    first_line <- line
    break
  }
}
# Ensure the file connection is closed after reading
on.exit(close(file_conn))

if (is.null(first_line)) {
  print("No non-commented lines found")
}
 
 # Check if the first non-commented line is empty or contains only whitespace
 if (nchar(trimws(first_line)) == 0) {
   cat("No VEP annotation found.\n")
   # If the file does not exist, create an empty data table with NA values
   # Define the column names
   column_names <- c(
     "Uploaded_variation", "Location", "Allele", "Gene", "Feature", "Feature_type",
     "Consequence", "cDNA_position", "CDS_position", "Protein_position",
     "Amino_acids", "Codons", "Existing_variation", "IMPACT", "DISTANCE",
     "STRAND", "FLAGS", "SYMBOL", "SYMBOL_SOURCE", "HGNC_ID", "BIOTYPE",
     "CANONICAL", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF",
     "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF",
     "gnomAD_SAS_AF", "CLIN_SIG", "SOMATIC", "PHENO", "MOTIF_NAME",
     "MOTIF_POS", "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "TRANSCRIPTION_FACTORS", "LEAD_rsID"
   )
   
   # Create the data table with NA values and set column names
   VEP_annotations <- data.table(matrix(NA, nrow = 0, ncol = length(column_names)))
   colnames(VEP_annotations) <- column_names
 } else {
   VEP_annotations <- data.table(read.table(file = file.path(VEP_filename), header = FALSE,
                                            quote = NULL, sep = "\t", stringsAsFactors = FALSE))
   #add collumn with lead snp
   p <- paste0(output_path,"_proxies.txt")
   proxies <- data.table(read.table(file = file.path(p), header = TRUE, quote = NULL, sep = "\t", stringsAsFactors = FALSE))
   proxies =as.data.frame(proxies)
   VEP_annotations <- merge(VEP_annotations, proxies[, c("PROXY_rsID", "LEAD_rsID")], by.x ="V1", by.y = "PROXY_rsID", all.x = TRUE) }
 
 
#print(head(VEP_annotations))
print(dim(VEP_annotations))

# filter out all records that do not have either a "HIGH" or "MODERATE" IMPACT
VEP_IMPACTS_of_interest <- c("HIGH", "MODERATE")
dim(VEP_annotations)
VEP_annotations <- VEP_annotations[VEP_annotations[[IMPACT_column]] %in% VEP_IMPACTS_of_interest,]
cat("VEP results that have a HIGH or MODERATE IMPACT :", nrow(VEP_annotations), "\n")
# identify leads and proxies most likely to have a functional impact as per VEP annotation
SNP_IMPACT_annotations <- VEP_integrator(lead_data, anno_summary, VEP_annotations, lead_rsID_column, proxy_rsID_column,
                                         ensembl_gene_id_column, IMPACT_column)
cat("SNP_IMPACT_annotations :", nrow(VEP_annotations), "\n")

head(SNP_IMPACT_annotations)
write.table(SNP_IMPACT_annotations, file = file.path(output_dir, "SNP_IMPACT_annotations.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
# combine anno_summary with the SNP_IMPACT_annotations
anno_summary <- cbind(anno_summary, SNP_IMPACT_annotations)

cat("\t- VEP annotations have been integrated.\n\n")

# COLOC data
# Coloc with eQTL

# read in coloc eQTL file if not empty
file_conn <- file(COLOC_EQTL_filename, "r")
 
 # Read the first non-commented line
 first_line <- readLines(file_conn, n = 1, warn = FALSE)
 
 # Close the file connection
 close(file_conn)
 
 # Check if the first non-commented line is empty or contains only whitespace
 if (nchar(trimws(first_line)) == 0) {
     cat("No coloc_eqtl  found.\n")
     column_names <- c("rsID", "trait", "Name", "sumstats_2_max_nlog10P", "PP.H4.abf", "Tissue", "gene_id", "gene_type", "gene_name","cis_trans")
     coloc_eqtl_data <- data.table(matrix(NA, nrow = 0, ncol = length(column_names)))
     colnames(coloc_eqtl_data) <- column_names
 } else {
    coloc_eqtl_data <- data.table(read.table(file =
                                             file.path(COLOC_EQTL_filename),
                                           header = TRUE, quote = NULL, sep = "\t", stringsAsFactors = FALSE))

    #if (length(eqtl_tissue_col) > 0) {
    #coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[eqtl_tissue_col]] %in% eqtl_tissues_of_interest, ]
    #}

    coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[eqtl_sumstats_2_max_nlog10P_col]] > eqtl_sumstats_2_max_nlog10P_thresh, ]
    coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[eqtl_PP.H4.abf_col]] > eqtl_PP.H4.abf_thresh, ]
    coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[coloc_eqtl_gene_type_col]] %in% biotype_of_interest, ]
    # filter out all records that do not have a "cis"
    coloc_eqtl_data <- coloc_eqtl_data[coloc_eqtl_data[[eqtl_cis_trans_col]] %in% eqtl_cis_trans_sel, ]

    coloc_eqtl_data <- as.data.frame(coloc_eqtl_data)
    coloc_eqtl_data <- coloc_eqtl_data[!duplicated(coloc_eqtl_data[,
                                                                 c(coloc_eqtl_lead_rsID_col, coloc_eqtl_ensembl_gene_id_col)]), ]
    coloc_eqtl_data <- data.table(coloc_eqtl_data)
}
# formatting so that it can be added to the bottom-up summary
COLOC_EQTL_annotations <- COLOC_integrator(lead_data, anno_summary,
                                         coloc_eqtl_data, coloc_eqtl_lead_rsID_col,
                                         coloc_eqtl_ensembl_gene_id_col)

# combine anno_summary with the COLOC EQTL annotation
anno_summary <- cbind(anno_summary, COLOC_EQTL_annotations)

cat("\t- Coloc eQTL information have been integrated.\n")

if (any(!is.na(tissues_interest))){
    print("eQTL tissues of interest for coloc")
    print(tissues_interest)
    #create a data.frame with tissues of interest
    coloc_eqtl_tissue_interest  <- coloc_eqtl_data[which(coloc_eqtl_data$Tissue %in% tissues_interest),]
    # formatting so that it can be added to the bottom-up summary
    COLOC_EQTL_tissues_int_annotations <- COLOC_integrator(lead_data, anno_summary, coloc_eqtl_tissue_interest, coloc_eqtl_lead_rsID_col,
                                             coloc_eqtl_ensembl_gene_id_col, QTL_type = "eQTL_tissues_interest")
    
    write.table(coloc_eqtl_tissue_interest, file = file.path(output_dir, "coloc_eqtl_tissue_interest.txt"),
                                                                                                  quote = FALSE, row.names = FALSE, sep = "\t")
    
    # combine anno_summary with the COLOC EQTL annotation
    anno_summary <- cbind(anno_summary, COLOC_EQTL_tissues_int_annotations)

    cat("\t- Coloc eQTL tissues of interest information have been integrated.\n")
}
  # read in coloc pQTL file
# read in coloc pQTL file if not empty
file_conn <- file(COLOC_PQTL_filename, "r")
 
 # Read the first non-commented line
 first_line <- readLines(file_conn, n = 1, warn = FALSE)
 
 # Close the file connection
 close(file_conn)
 
 # Check if the first non-commented line is empty or contains only whitespace
 if (nchar(trimws(first_line)) == 0) {
     cat("No coloc_pqtl  found.\n")
     column_names <- c("rsID", "trait", "Name", "sumstats_2_max_nlog10P", "PP.H4.abf", "gene_id","cis_trans")
     # Create the data table with NA values and set column names
     coloc_pqtl_data <- data.table(matrix(NA, nrow = 0, ncol = length(column_names)))
     colnames(coloc_pqtl_data) <- column_names
 } else {
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
                                                                 c(coloc_pqtl_lead_rsID_col, coloc_pqtl_ensembl_gene_id_col)]), ]
    coloc_pqtl_data <- data.table(coloc_pqtl_data)
}
# formatting so that it can be added to the bottom-up summary
COLOC_PQTL_annotations <- COLOC_integrator(lead_data, anno_summary,
                                         coloc_pqtl_data, coloc_pqtl_lead_rsID_col,
                                         coloc_pqtl_ensembl_gene_id_col, QTL_type = "pQTL")

# combine anno_summary with the COLOC EQTL annotation
anno_summary <- cbind(anno_summary, COLOC_PQTL_annotations)


cat("\t- Coloc pQTL information have been integrated.\n\n")

# classify anno evidence according to source
# 1 = gene contains a lead or proxy with either a HIGH or MODERATE IMPACT (proximity and cis-eQTL evidence may also be present)
# 2 = gene is proximal to lead; i.e., nearest or LD overlapping (cis-eQTL evidence may also be present)
# 3 = gene is a cis-eQTL target of either the lead or a proxy variant



cat("Annotation complete!\n\n")



### WRITE INFORMATION TO FILE

# bottom-up summary
write.table(anno_summary, file = file.path(output_dir, "OUTPUT_anno_summary.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")


# HIGH and MODERATE VEP dataframe
write.table(VEP_annotations, file = file.path(output_dir, "OUTPUT_VEP_IMPACT_annotations.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# coloc_eQTL
write.table(coloc_eqtl_data, file = file.path(output_dir, "OUTPUT_COLOC_EQTL.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
# coloc_pQTL
write.table(coloc_pqtl_data, file = file.path(output_dir, "OUTPUT_COLOC_PQTL.txt"),
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
