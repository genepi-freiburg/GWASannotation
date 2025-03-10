library(optparse)
library(writexl)
library(dplyr)
library(readxl)
library(data.table)


option_list <- list(
  make_option("--output_path", action="store", default=NA, type='character', help="output_path for GWASAnno full path [required]"),
  make_option("--output_file_name", action="store", default="GWASAnno_summary.txt", type='character', help="output file name [default='GWASAnno_summary.txt]'"),
  make_option("--eQTL_tissues_interest_coloc", action="store", default=NA, type='character', help="Comma-separated tissues of interest previously selected from the eQTL datasets [default=NA]"),
  make_option("--nearest", action="store", default="1", type='numeric', help="nearest score [default=1]"),
  make_option("--second_nearest", action="store", default="0.99", type='numeric', help="second_nearest score [default=0.99]"),
  make_option("--third_nearest", action="store", default="0.98", type='numeric', help="third_nearest score [default=0.98]"),
  make_option("--LD_overlapping", action="store", default="1", type='numeric', help="LD_overlapping score [default=1"),
  make_option("--lead_IMPACT", action="store", default="1.02", type='numeric', help="lead_IMPACT score [default=1.02]"),
  make_option("--coloc_eQTL_tissues_interest", action="store", default="1.04", type='numeric', help="coloc_eQTL_tissues_interest score (if --eQTL_tissues_interest_coloc are not defined, --coloc_eQTL_tissues_interest is not used) [default=1.04]"),
  make_option("--coloc_eQTL", action="store", default="1.03", type='numeric', help="coloc_eQTL score (if --eQTL_tissues_interest_coloc is defined, this score is only applied if NO coloc_eQTL_tissues_interest is identified) [default=1.03]"),
  make_option("--lead_eQTL", action="store", default="0.75", type='numeric', help="lead_eQTL or proxie_eQTL score (this score is only applied if NO coloc_eQTL is identified) [default=0.75]"),
  make_option("--coloc_pQTL", action="store", default="1.04", type='numeric', help="coloc_pQTL score [default=1.04]")
)


opt_used <- FALSE
tryCatch({
  opt <- parse_args(OptionParser(option_list=option_list))
  output_path <- opt$output_path
  folder_path <- dirname(opt$output_path)
  output_dir <- paste0(folder_path, "/GWASAnno/")
  print(output_path)
  tissues_interest <- strsplit(opt$eQTL_tissues_interest_coloc, ",")[[1]]

  # Reading and processing for anno with optparse
  anno.file <- paste0(output_dir, "OUTPUT_anno_summary.txt")
  tophit.file <- paste0(opt$output_path,"_lead.txt")
  summary <- paste0(output_dir, opt$output_file_name)
  #scoring
  nearest_score <-  opt$nearest
  second_nearest_score <-  opt$second_nearest
  third_nearest_score <-  opt$third_nearest
  LD_overlapping_score <-  opt$LD_overlapping
  lead_IMPACT_score <- opt$lead_IMPACT
  coloc_eQTL_tissues_interest_score <- opt$coloc_eQTL_tissues_interest
  coloc_eQTL <- opt$coloc_eQTL
  lead_eQTL_score <- opt$lead_eQTL
  coloc_pQTL_score <- opt$coloc_pQTL
  opt_used <- TRUE
}, error = function(e) {
  cat("")
}, finally = {
    if (!opt_used) {
        # Reading and processing for anno
        anno.file <- paste0(output_dir, "OUTPUT_anno_summary.txt")
        tophit.file <- lead_filename
        summary <- paste0(output_dir, "GWASAnno_summary.txt")
        #scoring
        nearest_score <-  1
        second_nearest_score <-  0.99
        third_nearest_score <-  0.98
        LD_overlapping_score <- 1
        lead_IMPACT_score <- 1.02
        coloc_eQTL_tissues_interest_score <- 1.04
        coloc_eQTL <- 1.03
        lead_eQTL_score <- 0.75
        coloc_pQTL_score <- 1.04
        anno <- read.table(anno.file, header=TRUE, sep="\t")
        anno[anno$hgnc_symbol != "-", ]
        empty_symbol <- which(anno$hgnc_symbol == "-")
        anno$hgnc_symbol[empty_symbol] <- anno$ensembl_id[empty_symbol]

    }
})


anno <- read.table(anno.file, header=TRUE, sep="\t")
anno[anno$hgnc_symbol != "-", ]
empty_symbol <- which(anno$hgnc_symbol == "-")
anno$hgnc_symbol[empty_symbol] <- anno$ensembl_id[empty_symbol]

#If tissues_interest are selected"
if (any(!is.na(tissues_interest))){
    print("#check if it already exist a file with tissues_interest info: coloc_eqtl_tissue_interest")
    if (file.exists(paste0(output_dir, "coloc_eqtl_tissue_interest.txt"))){
        ##check if tissues_interest are the same
        tissues_done <- read.table(paste0(output_dir, "coloc_eqtl_tissue_interest.txt"), header=T)
        tissues_done <- unique(tissues_done$Tissue)
        if (!all(tissues_interest %in% tissues_done)) {
          cat("tissues of interest not found in coloc_eqtl_tissue_interest.txt - re-check. \n")
          args <- commandArgs(trailingOnly = FALSE)
          script_path <- sub("--file=", "", args[grep("--file=", args)])
          script_path <- normalizePath(script_path)
          script_path <- dirname(script_path)
          source(file = file.path(paste0(dirname(script_path),"/04_utils/Anno_functions.R")))
          ####create new tissues_interest
          coloc_eqtl_data <- read.table(paste0(output_dir, "OUTPUT_COLOC_EQTL.txt"), header=T)
          coloc_eqtl_tissue_interest  <- coloc_eqtl_data[which(coloc_eqtl_data$Tissue %in% tissues_interest),]
          if (nrow(coloc_eqtl_tissue_interest)>0){
              cat("Re-writing as coloc_eqtl_tissue_interest_v2.txt \n")
              write.table(coloc_eqtl_tissue_interest, file = file.path(output_dir, "coloc_eqtl_tissue_interest_v2.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
          }
          # formatting so that it can be added to the bottom-up summary
          lead_data <- read.table(tophit.file, header=T)
          COLOC_EQTL_tissues_int_annotations <- COLOC_integrator(lead_data, anno[, c("LEAD_rsID", "ensembl_id", "hgnc_symbol", "gene_biotype")], coloc_eqtl_tissue_interest, 1, 7, QTL_type = "eQTL_tissues_interest")
          anno$coloc_eQTL_tissues_interest <- COLOC_EQTL_tissues_int_annotations$coloc_eQTL_tissues_interest
        } else {
            cat("tissues of interest selected \n")
        }
    } else {
        cat("file with tissues_interest info doesnt exist. Creating it from OUTPUT_COLOC_EQTL (only works after Feb. 2025 update)! \n")
        args <- commandArgs(trailingOnly = FALSE)
        script_path <- sub("--file=", "", args[grep("--file=", args)])
        script_path <- normalizePath(script_path)
        script_path <- dirname(script_path)
        source(file = file.path(paste0(dirname(script_path),"/04_utils/Anno_functions.R")))
        ####create new tissues_interest
        coloc_eqtl_data <- read.table(paste0(output_dir, "OUTPUT_COLOC_EQTL.txt"), header=T)
        coloc_eqtl_tissue_interest  <- coloc_eqtl_data[which(coloc_eqtl_data$Tissue %in% tissues_interest),]
        if (nrow(coloc_eqtl_tissue_interest)>0){
            cat("Re-writing as coloc_eqtl_tissue_interest_v2.txt \n")
            write.table(coloc_eqtl_tissue_interest, file = file.path(output_dir, "coloc_eqtl_tissue_interest_v2.txt"), quote = FALSE, row.names = FALSE, sep = "\t")
        }
        # formatting so that it can be added to the bottom-up summary
        lead_data <- read.table(tophit.file, header=T)
        COLOC_EQTL_tissues_int_annotations <- COLOC_integrator(lead_data, anno[, c("LEAD_rsID", "ensembl_id", "hgnc_symbol", "gene_biotype")], coloc_eqtl_tissue_interest, 1, 7, QTL_type = "eQTL_tissues_interest")
        anno$coloc_eQTL_tissues_interest <- COLOC_EQTL_tissues_int_annotations$coloc_eQTL_tissues_interest
    }
}

if(nrow(anno)>0){
    anno$nearest[is.na(anno$nearest)] <- 0
    anno$second_nearest[is.na(anno$second_nearest)] <- 0
    anno$third_nearest[is.na(anno$third_nearest)] <- 0
    anno$LD_overlapping[is.na(anno$LD_overlapping)] <- 0
    anno$lead_eQTL[is.na(anno$lead_eQTL)] <- 0
    anno$proxy_eQTL[is.na(anno$proxy_eQTL)] <- 0
    anno$lead_IMPACT[is.na(anno$lead_IMPACT)] <- 0
    anno$proxy_IMPACT[is.na(anno$proxy_IMPACT)] <- 0


    anno$coloc_eQTL[is.na(anno$coloc_eQTL)] <- 0
    if (any(!is.na(tissues_interest))){
        #IF IT WASNT PREVIOUSLY SELECTED TISSUES OF INTEREST NEED TO BE CREATED NOW
        anno$coloc_eQTL_tissues_interest[is.na(anno$coloc_eQTL_tissues_interest)] <- 0
    }
    anno$coloc_pQTL[is.na(anno$coloc_pQTL)] <- 0
    print(head(anno))

    out.sum <- data.frame()

    for (rsid in unique(anno$LEAD_rsID)) {
        cat("rsid: ", rsid, "\n")
        genes.df <- data.frame()
        locus <- anno[anno$LEAD_rsID== rsid, ]
        head(locus)
        #print(dim(locus))
        for (i in 1:nrow(locus)) {
            gene <- locus[i, "hgnc_symbol"]
            # print(gene)

            evidence <- c()
            cnt <- 0
            if (locus[i, "nearest"] == 1) {
                evidence <- c(evidence, "nearest")
                cnt <- cnt + nearest_score
            }
            if (locus[i, "second_nearest"] == 1) {
                evidence <- c(evidence, "2nd_nearest")
                cnt <- cnt + second_nearest_score
            }
            if (locus[i, "third_nearest"] == 1) {
                evidence <- c(evidence, "3rd_nearest")
                cnt <- cnt + third_nearest_score
            }
            if (locus[i, "LD_overlapping"] == 1) {
                evidence <- c(evidence, "LD_overlap")
                cnt <- cnt + LD_overlapping_score
            }
            
            if (locus[i, "lead_IMPACT"]==1 || locus[i, "proxy_IMPACT"]==1) {
                evidence <- c(evidence, "IMPACT_moderate_high")
                cnt <- cnt + lead_IMPACT_score
            }
            #if including coloc_eQTL_tissue_interest
            if (any(!is.na(tissues_interest))){
                if (locus[i, "coloc_eQTL_tissues_interest"] == 1) {
                    evidence <- c(evidence, "coloc_eQTL_tissues_interest")
                    cnt <- cnt + coloc_eQTL_tissues_interest_score
                } else {
                    if (locus[i, "coloc_eQTL"] == 1) {
                        evidence <- c(evidence, "coloc_eQTL")
                        cnt <- cnt + coloc_eQTL
                    } else {
                      if (locus[i, "lead_eQTL"]==1 | locus[i, "proxy_eQTL"]==1) {
                        evidence <- c(evidence, "cis-eQTL")
                        cnt <- cnt + lead_eQTL_score
                      }
                    }
                }
            } else {
                if (locus[i, "coloc_eQTL"] == 1) {
                    evidence <- c(evidence, "coloc_eQTL")
                    cnt <- cnt + coloc_eQTL
                }else {
                  if (locus[i, "lead_eQTL"]==1 | locus[i, "proxy_eQTL"]==1) {
                    evidence <- c(evidence, "cis-eQTL")
                    cnt <- cnt + lead_eQTL_score
                  }
                }
            }
            if (locus[i, "coloc_pQTL"] == 1) {
                evidence <- c(evidence, "coloc_pQTL")
                cnt <- cnt + coloc_pQTL_score
            }
         gene.df <- data.frame(symbol=gene,
                evidences =paste(evidence, collapse=", "), count =cnt)
            genes.df <- rbind(genes.df, gene.df)
        }

        
        #custom_order <- c("nearest", "coloc_eQTL_tissues_interest" ,"coloc_eQTL", "coloc_pQTL", "LD_overlapping", "IMPACT_moderate_high")
     
        # Group by count and reorder within each group
       # genes.df <- genes.df %>%
       #   group_by(count) %>%
       #   arrange(factor(evidences, levels = custom_order))
        
        genes.df <- genes.df[order(genes.df$count, decreasing = TRUE), ]
        # Print the resulting dataframe
        #print(genes.df)
        
        # Order based on evidence type
        #subset_order <- order(factor(evidence_type, levels = custom_order))

        #genes.df[genes.df$count == 1, ] <- genes.df[genes.df$count == 1, ][subset_order, ]
        #print(genes.df[genes.df$count ==1,])
        
        
        cat("genes and evidences \n")
        print(genes.df)

        genes.df$new.text <- paste0(genes.df$symbol, " (", genes.df$evidences, ")", " *score: ", genes.df$count)
        new.df <- data.frame(rsID = rsid,
                  gene_evidences_all = paste(genes.df$new.text, collapse="; "),
                  gene_evidences_top1 = genes.df$new.text[1], gene_evidences_top2 = genes.df$new.text[2], gene_evidences_top3 = genes.df$new.text[3])
        out.sum <- rbind(out.sum, new.df)
        #print(out.sum)
        # print(paste(genes.info, collapse="; "))
    }
    regions_file <- paste0(output_path, "_coloc_regions.RDS")
    regions_list <- readRDS(regions_file)
    regions <- regions_list$coloc_regions
    regions <- regions[which(regions$comment=="PASS"),]
    regions$P <- 10^(-regions$nlog10P)
    regions <- regions[,c("CHR_var", "BP_START_var", "BP_STOP_var","rsID","P","BETA","SE","AF")]
    regions <-  as.data.table(regions)
    out.sum <-  as.data.table(out.sum)
    #print("test")
    #print(head(regions))
    #print(head(out.sum))
    out.sum <- merge(regions, out.sum, by="rsID", all.x = TRUE)

} else {
    cat("no gene annotation was found for this dataset\n")
    regions_file <- paste0(output_path, "_coloc_regions.RDS")
    regions_list <- readRDS(regions_file)
    regions <- regions_list$coloc_regions
    regions <- regions[which(regions$comment=="PASS"),]
    regions$P <- 10^(-regions$nlog10P)
    regions <- regions[,c("CHR_var", "BP_START_var", "BP_STOP_var","rsID","P","BETA","SE","AF")]
    out.sum <-  as.data.table(regions)
    out.sum$gene_evidences_all <- NA
    out.sum$gene_evidences_top1 <- NA
    out.sum$gene_evidences_top2 <- NA
    out.sum$gene_evidences_top3 <- NA
}


write.table(x= out.sum, file = summary, sep="\t", row.names = FALSE)
outxlsx <- sub("\\.txt$", ".xlsx", summary)

write_xlsx(x=out.sum, path = outxlsx, format_headers = TRUE)
cat("summary file done \n")


#write_xlsx(x=m, path = paste0(output_dir, "ProGeM_lead.xlsx"), format_headers = TRUE)
cat("## postprocessing finished ## \n")
