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
  make_option("--coloc_pQTL", action="store", default="1.04", type='numeric', help="coloc_pQTL score [default=1.04]"),
  make_option("--PoPS_top1", action="store", default="0", type='numeric', help="PoPS_top1 score [default=0]"),
  make_option("--PoPS_top2", action="store", default="0", type='numeric', help="PoPS_top2 score [default=0]"),
  make_option("--PoPS_top3", action="store", default="0", type='numeric', help="PoPS_top3 score [default=0]")
)


opt_used <- FALSE
tryCatch({
  opt <- parse_args(OptionParser(option_list=option_list))
  output_path <- opt$output_path
  folder_path <- dirname(opt$output_path)
  output_dir <- paste0(folder_path, "/GWASAnno/")
  print(output_path)
  tissues_interest <- opt$eQTL_tissues_interest_coloc
  
  # Reading and processing for bottom_up, top_down, and pops with optparse
  bottom_up.file <- paste0(output_dir, "OUTPUT_bottom_up_summary.txt")
  top_down.file <- paste0(output_dir,  "OUTPUT_top_down_scores.txt")
  pops.file <- paste0(opt$output_path,"_PoPS.RData")
  tophit.file <- paste0(opt$output_path,"_sentinel.txt")
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
  PoPS_top1_score <- opt$PoPS_top1
  PoPS_top2_score <- opt$PoPS_top2
  PoPS_top3_score <- opt$PoPS_top3
  opt_used <- TRUE
}, error = function(e) {
  cat("")
}, finally = {
    if (!opt_used) {
        # Reading and processing for bottom_up, top_down, and pops
        bottom_up.file <- paste0(output_dir, "OUTPUT_bottom_up_summary.txt")
        top_down.file <- paste0(output_dir, "OUTPUT_top_down_scores.txt")
        pops.file <- filename_PoPS
        tophit.file <- sentinel_filename
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
        PoPS_top1_score <- 0
        PoPS_top2_score <- 0
        PoPS_top3_score <- 0
    }
})

bottom_up <- read.table(bottom_up.file, header=TRUE, sep="\t")
bottom_up[bottom_up$hgnc_symbol != "-", ]
empty_symbol <- which(bottom_up$hgnc_symbol == "-")
bottom_up$hgnc_symbol[empty_symbol] <- bottom_up$ensembl_id[empty_symbol]
#dim(bottom_up)
print("head bottom_up")
print(head(bottom_up))
if(file.exists(top_down.file)){
    top_down <- read.table(top_down.file, header=TRUE, sep="\t")
    #dim(top_down)
    #head(top_down)

    top_down[top_down$hgnc_symbol == "-", ]
    empty_symbol <- which(top_down$hgnc_symbol == "-")
    top_down$hgnc_symbol[empty_symbol] <- top_down$ensembl_id[empty_symbol]

    top_down$PoPS_top3 <- top_down$PoPS_top2 <- top_down$PoPS_top1 <- 0
    #im(top_down)
    head(top_down)

    for (rsid in unique(top_down$LEAD_rsID)) {
        locus <- top_down[top_down$LEAD_rsID== rsid, ]
        #print(dim(locus))

        locus <- locus[order(locus$PoPS_score, decreasing = TRUE), ]
        locus <- head(locus, n=3)
        #print(dim(locus))

        top_down[top_down$LEAD_rsID== rsid &
            top_down$ensembl_id == locus$ensembl_id[1], "PoPS_top1"] <- 1

        if (nrow(locus) >=2) {
            top_down[top_down$LEAD_rsID== rsid &
                top_down$ensembl_id == locus$ensembl_id[2], "PoPS_top2"] <- 1
        }

        if (nrow(locus) >=3) {
            top_down[top_down$LEAD_rsID== rsid &
                top_down$ensembl_id == locus$ensembl_id[3], "PoPS_top3"] <- 1
        }
    }

    top_down <- top_down[top_down$PoPS_top1 != 0 | top_down$PoPS_top2 != 0 |
        top_down$PoPS_top3 != 0, ]
    #dim(top_down)

    intersect(names(bottom_up), names(top_down))
    m <- merge(bottom_up, top_down, all=TRUE)
    m$PoPS_top1[is.na(m$PoPS_top1)] <- 0
    m$PoPS_top2[is.na(m$PoPS_top2)] <- 0
    m$PoPS_top3[is.na(m$PoPS_top3)] <- 0
} else {
    cat("TOPdown annotation not found, including only Bottom-up \n")
    m <- bottom_up
    print("head m")
    print(head(m))
    m$PoPS_top1 <- 0
    m$PoPS_top2 <- 0
    m$PoPS_top3 <- 0
    #print(head(m))
}
#dim(m)
#names(m)

# m[is.na(m)] <- 0
# m$cis_eQTL <- ifelse(m$lead_eQTL==1 | m$proxy_eQTL==1, 1, 0)
# m$IMPACT_moderate_high <- ifelse(m$lead_IMPACT== 1 | m$proxy_IMPACT==1, 1, 0)
# m$category <- m$category - m$score

m$nearest[is.na(m$nearest)] <- 0
m$second_nearest[is.na(m$second_nearest)] <- 0
m$third_nearest[is.na(m$third_nearest)] <- 0
m$LD_overlapping[is.na(m$LD_overlapping)] <- 0
m$lead_eQTL[is.na(m$lead_eQTL)] <- 0
m$proxy_eQTL[is.na(m$proxy_eQTL)] <- 0
m$lead_IMPACT[is.na(m$lead_IMPACT)] <- 0
m$proxy_IMPACT[is.na(m$proxy_IMPACT)] <- 0


m$coloc_eQTL[is.na(m$coloc_eQTL)] <- 0
if (any(!is.na(tissues_interest))){
    m$coloc_eQTL_tissues_interest[is.na(m$coloc_eQTL_tissues_interest)] <- 0
}
m$coloc_pQTL[is.na(m$coloc_pQTL)] <- 0



# if pops is not NA and category is not NA: category - pops score (1)
# pops.ind <- which(m$category>= min.cat & m$PoPS==1)
# m[pops.ind,]
# (m$category[pops.ind] <- m$category[pops.ind] - m$PoPS[pops.ind])

# set category to the maximum category value when it only has PoPS
# m$category[is.na(m$category) & m$PoPS==1] <- max.cat

# df.order <- with(m, order(LEAD_rsID, category, -nearest, -second_nearest,
#     -third_nearest, - Score))
# m <- m[df.order, ]

out.sum <- data.frame()

for (rsid in unique(m$LEAD_rsID)) {
    cat("rsid: ", rsid, "\n")
    genes.df <- data.frame()
    locus <- m[m$LEAD_rsID== rsid, ]
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

        if (locus[i, "PoPS_top1"] == 1) {
            evidence <- c(evidence,
                paste0("PoPS_top1_", round(locus[i, "PoPS_score"], 2)))
            cnt <- cnt + PoPS_top1_score
        }
        if (locus[i, "PoPS_top2"] == 1) {
            evidence <- c(evidence,
                paste0("PoPS_top2_", round(locus[i, "PoPS_score"], 2)))
            cnt <- cnt + PoPS_top2_score
        }
        if (locus[i, "PoPS_top3"] == 1) {
            evidence <- c(evidence,
                paste0("PoPS_top3_", round(locus[i, "PoPS_score"], 2)))
            cnt <- cnt + PoPS_top3_score
        }

        gene.df <- data.frame(symbol=gene,
            evidences =paste(evidence, collapse=", "), count =cnt)
        genes.df <- rbind(genes.df, gene.df)
    }

    
    print("here")
    custom_order <- c("nearest", "coloc_eQTL_tissues_interest" ,"coloc_eQTL", "coloc_pQTL", "LD_overlapping", "IMPACT_moderate_high", "PoPS_top1", "PoPS_top2", "PoPS_top3") 
    #evidence_type <- sub("^PoPS_top1_([-0-9.]+)$", "PoPS_top1", genes.df[genes.df$count == 1, ]$evidences)
    #evidence_type <- sub("^PoPS_top2_([-0-9.]+)$", "PoPS_top2", evidence_type)
    #evidence_type <- sub("^PoPS_top3_([-0-9.]+)$", "PoPS_top3", evidence_type)
    
    # Group by count and reorder within each group
    genes.df <- genes.df %>%
      group_by(count) %>%
      arrange(factor(evidences, levels = custom_order))
    
    genes.df <- genes.df[order(genes.df$count, decreasing = TRUE), ]
    # Print the resulting dataframe
    print(genes.df)
    
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


write.table(x= out.sum, file = summary, sep="\t", row.names = FALSE)
outxlsx <- sub("\\.txt$", ".xlsx", summary)

write_xlsx(x=out.sum, path = outxlsx, format_headers = TRUE)
cat("summary file done \n")

#tophit=read.table(tophit.file, header=T)
#cat("sentinel file \n")
#dim(tophit)
#head(tophit)
#tophit<- tophit %>%
 #   mutate(rsID = ifelse(startsWith(rsID, "rs"), rsID, paste0(CHR, ":", START)))
#cat("merging ProGeM summary file with sentinel file \n")
#progem <- out.sum


#m <- merge(tophit, progem,  by="rsID", sort=FALSE)
#dim(m)
#names(m)

#write_xlsx(x=m, path = paste0(output_dir, "ProGeM_sentinel.xlsx"), format_headers = TRUE)
#cat("ProGeM_sentinel file done \n")
cat("## postprocessing finished ## \n")
