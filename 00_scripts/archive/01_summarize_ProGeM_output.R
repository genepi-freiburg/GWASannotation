suppressMessages(library(optparse))
library(writexl)
option_list = list(
  make_option("--bottom_up.file", action="store", default=NA, type='character', help="bottom_up.file full path [required]"),
  make_option("--top_down.file", action="store", default=NA, type='character', help="top_down.file full path [required]"),
  make_option("--pops.file", action="store", default=NA, type='character', help="pops.file full path [required]"),
  make_option("--out.file", action="store", default=NA, type='character', help="out.file full path [required]")
)

opt = parse_args(OptionParser(option_list=option_list))



bottom_up.file <- opt$bottom_up.file
top_down.file <- opt$top_down.file
pops.file <- opt$pops.file
out.file <- opt$out.file

#bottom_up.file <- "/data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/02_output/model1_qnorm_tkv.bsa/OUTPUT_bottom_up_summary.txt"
#top_down.file <- "/data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/02_output/model1_qnorm_tkv.bsa/OUTPUT_top_down_scores.txt"
#pops.file <- "/data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/02_Pops/02_output/model1_qnorm_tkv.bsa/model1_qnorm_tkv.bsa_all_chr.results"
#out.file <- "/data/studies/06_UKBB/02_Projects/14_MRI-kidney/07_ProGeM/02_output/model1_qnorm_tkv.bsa/model1_qnorm_tkv.bsa_ProGeM_summary.txt"


bottom_up <- read.table(bottom_up.file, header=TRUE, sep="\t")
dim(bottom_up)
head(bottom_up)

(min.cat <- min(bottom_up$category))
(max.cat <- max(bottom_up$category))

top_down <- read.table(top_down.file, header=TRUE, sep="\t")
dim(top_down)
head(top_down)

pops <- read.table(pops.file, header=TRUE, sep="\t")
dim(pops)
head(pops)

top_down <- merge(top_down, pops, by.x="ensembl_id", by.y="ENSGID")
dim(top_down)
head(top_down)
rm(pops)

intersect(names(bottom_up), names(top_down))
m <- merge(bottom_up, top_down, all=TRUE)
dim(m)
names(m)

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

# if pops is not NA and category is not NA: category - pops score (1)
pops.ind <- which(m$category>= min.cat & m$PoPS==1)
m[pops.ind,]
(m$category[pops.ind] <- m$category[pops.ind] -  m$PoPS[pops.ind])

# set category to the maximum category value when it only has PoPS
m$category[is.na(m$category) & m$PoPS==1] <- max.cat

df.order <- with(m, order(LEAD_rsID, category, -nearest, -second_nearest,
    -third_nearest, - Score))
m <- m[df.order, ]

out.df <- data.frame()

for (rsid in unique(m$LEAD_rsID)) {
    # genes.info <- c()
    genes.df <- data.frame()
    locus <- m[m$LEAD_rsID== rsid, ]
    print(dim(locus))

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
        if (locus[i, "lead_eQTL"]==1 | locus[i, "proxy_eQTL"]==1) {
            evidence <- c(evidence, "cis-eQTL")
            cnt <- cnt + 1
        }
        if (locus[i, "lead_IMPACT"]==1 || locus[i, "proxy_IMPACT"]==1) {
            evidence <- c(evidence, "IMPACT_moderate_high")
            cnt <- cnt + 1
        }
        if (!is.na(locus[i, "PoPS"]) & locus[i, "PoPS"] == 1) {
            evidence <- c(evidence,
                paste0("PoPS_", round(locus[i, "Score"], 2)))
            cnt <- cnt + 1
        }
        # print(paste0(gene, " (", paste(evidence, collapse=", "), ")"))

        gene.df <- data.frame(symbol=gene,
            evidences =paste(evidence, collapse=", "),
            count =cnt)
        genes.df <- rbind(genes.df, gene.df)
        # genes.info <- c(genes.info,
        #     paste0(gene, " (", paste(evidence, collapse=", "), ")"))
    }

    genes.df <- genes.df[order(genes.df$count, decreasing = TRUE), ]
    print(genes.df)

    genes.df$new.text <- paste0(genes.df$symbol, " (", genes.df$evidences, ")")
    new.df <- data.frame(rsID = rsid,
              gene_evidences_all = paste(genes.df$new.text, collapse="; "),
              gene_evidences_top1 = genes.df$new.text[1])
    out.df <- rbind(out.df, new.df)
    print(out.df)

    # print(paste(genes.info, collapse="; "))
}

write.table(x= out.df, file = out.file, sep="\t", row.names = FALSE)
outxlsx <- sub("\\.txt$", ".xlsx", out.file)

write_xlsx(x=out.df, path = outxlsx, format_headers = TRUE)
print(out.file)
