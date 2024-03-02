# Reading and processing for bottom_up, top_down, and pops
bottom_up.file <- paste0(output_dir, "OUTPUT_bottom_up_summary.txt")
top_down.file <- paste0(output_dir, "OUTPUT_top_down_scores.txt")
pops.file <- filename_PoPS
tophit.file <- sentinel_filename
summary <- paste0(output_dir, "GWASAnno_summary.txt")



bottom_up <- read.table(bottom_up.file, header=TRUE, sep="\t")
bottom_up[bottom_up$hgnc_symbol != "-", ]
empty_symbol <- which(bottom_up$hgnc_symbol == "-")
bottom_up$hgnc_symbol[empty_symbol] <- bottom_up$ensembl_id[empty_symbol]
#dim(bottom_up)
head(bottom_up)

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

m$PoPS_top1[is.na(m$PoPS_top1)] <- 0
m$PoPS_top2[is.na(m$PoPS_top2)] <- 0
m$PoPS_top3[is.na(m$PoPS_top3)] <- 0

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
            cnt <- cnt + 1
        }
        if (locus[i, "second_nearest"] == 1) {
            evidence <- c(evidence, "2nd_nearest")
            cnt <- cnt + 0.5
        }
        if (locus[i, "third_nearest"] == 1) {
            evidence <- c(evidence, "3rd_nearest")
            cnt <- cnt + 0.25
        }
        if (locus[i, "LD_overlapping"] == 1) {
            evidence <- c(evidence, "LD_overlap")
            cnt <- cnt + 1
        }
        
        if (locus[i, "lead_IMPACT"]==1 || locus[i, "proxy_IMPACT"]==1) {
            evidence <- c(evidence, "IMPACT_moderate_high")
            cnt <- cnt + 1
        }
        #if including coloc_eQTL_tissue_interest
        #if(!is.na(tissues_interest)){
        if (any(!is.na(tissues_interest))){
            if (locus[i, "coloc_eQTL_tissues_interest"] == 1) {
                evidence <- c(evidence, "coloc_eQTL_tissues_interest")
                cnt <- cnt + 1
            } else {
                if (locus[i, "coloc_eQTL"] == 1) {
                    evidence <- c(evidence, "coloc_eQTL")
                    cnt <- cnt + 0.5
                }
            }
        } else {
            if (locus[i, "coloc_eQTL"] == 1) {
                evidence <- c(evidence, "coloc_eQTL")
                cnt <- cnt + 1
            }
            if (locus[i, "coloc_eQTL"] == 0) {
                if (locus[i, "lead_eQTL"]==1 | locus[i, "proxy_eQTL"]==1) {
                    evidence <- c(evidence, "cis-eQTL")
                    cnt <- cnt + 0.25
                }
            }
        }
        if (locus[i, "coloc_pQTL"] == 1) {
            evidence <- c(evidence, "coloc_pQTL")
            cnt <- cnt + 1
        }

        if (locus[i, "PoPS_top1"] == 1) {
            evidence <- c(evidence,
                paste0("PoPS_top1_", round(locus[i, "PoPS_score"], 2)))
            cnt <- cnt + 1
        }
        if (locus[i, "PoPS_top2"] == 1) {
            evidence <- c(evidence,
                paste0("PoPS_top2_", round(locus[i, "PoPS_score"], 2)))
            cnt <- cnt + 1
        }
        if (locus[i, "PoPS_top3"] == 1) {
            evidence <- c(evidence,
                paste0("PoPS_top3_", round(locus[i, "PoPS_score"], 2)))
            cnt <- cnt + 1
        }

        gene.df <- data.frame(symbol=gene,
            evidences =paste(evidence, collapse=", "), count =cnt)
        genes.df <- rbind(genes.df, gene.df)
    }

    genes.df <- genes.df[order(genes.df$count, decreasing = TRUE), ]

    custom_order <- c("nearest", "coloc_eQTL_tissues_interest" ,"coloc_eQTL", "coloc_pQTL", "LD_overlapping", "IMPACT_moderate_high", "PoPS_top1", "PoPS_top2", "PoPS_top3") 
    evidence_type <- sub("^PoPS_top1_([-0-9.]+)$", "PoPS_top1", genes.df[genes.df$count == 1, ]$evidences)
    evidence_type <- sub("^PoPS_top2_([-0-9.]+)$", "PoPS_top2", evidence_type)
    evidence_type <- sub("^PoPS_top3_([-0-9.]+)$", "PoPS_top3", evidence_type)

    # Order based on evidence type
    subset_order <- order(factor(evidence_type, levels = custom_order))

    genes.df[genes.df$count == 1, ] <- genes.df[genes.df$count == 1, ][subset_order, ]
    print(genes.df[genes.df$count ==1,])
    
    
    cat("genes and evidences \n")
    print(genes.df)

    genes.df$new.text <- paste0(genes.df$symbol, " (", genes.df$evidences, ")", " *score: ", genes.df$count)
    new.df <- data.frame(rsID = rsid,
              gene_evidences_all = paste(genes.df$new.text, collapse="; "),
              gene_evidences_top1 = genes.df$new.text[1], gene_evidences_top2 = genes.df$new.text[2], gene_evidences_top2 = genes.df$new.text[3])
    out.sum <- rbind(out.sum, new.df)
    #print(out.sum)
    # print(paste(genes.info, collapse="; "))
}


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
