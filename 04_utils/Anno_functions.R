GR_object_creator <- function(CHR, START, END, NAMES) {
  ranges_object <- GRanges(seqnames = CHR,
                           ranges = IRanges(start = START, end = END))
  names(ranges_object) <- NAMES
  return(ranges_object)
}


LD_region_range_finder <- function(lead_table, proxy_table, overhang) {
  unique_lead <- unique(lead_table$rsID)
  LD_ranges <- matrix(nrow = length(unique_lead), ncol = 3,
                      dimnames = list(unique_lead, c("CHR", "START", "END")))
  for(i in 1:length(unique_lead)) {
    #print(i)
    LD_ranges[i, 1] <- lead_table$CHR[which(lead_table$rsID == unique_lead[i])]
    LD_ranges[i, 2] <- min(proxy_table$PROXY_START[which(proxy_table$LEAD_rsID == unique_lead[i])],
                           lead_table$START[which(lead_table$rsID == unique_lead[i])])
    LD_ranges[i, 3] <- max(proxy_table$PROXY_END[which(proxy_table$LEAD_rsID == unique_lead[i])],
                           lead_table$END[which(lead_table$rsID == unique_lead[i])])
  }
  GR_object_creator(LD_ranges[,1], as.numeric(LD_ranges[,2]), as.numeric(LD_ranges[,3]), 
                    rownames(LD_ranges))
}


find_overlapping_genes <- function(ranges, gene_model, overhang) {
  hits <- findOverlaps(ranges, gene_model, maxgap = overhang * 1000) # convert kb to bp
  overlapping_genes <- gene_model[subjectHits(hits)]
  mcols(overlapping_genes) <- names(ranges[queryHits(hits)])
  names(mcols(overlapping_genes)) <- "LEAD_rsID"
  return(overlapping_genes)
}


gene_annotator <- function(ids) {
    # Rewritten - use ensembl_genes.df for lookup instead of biomaRt
  gene_info <- ensembl_genes.df[unique(ids), c("gene_id", "hgnc_symbol",
     "gene_biotype", "seqnames")]
  names(gene_info)[names(gene_info) == "gene_id"] <- "ensembl_gene_id"
  names(gene_info)[names(gene_info) == "seqnames"] <- "chromosome_name"
  gene_info <- data.table(gene_info)
  
  # format the gene_info object slightly
  if(length(grep("CHR", gene_info$chromosome_name))) {
   gene_info <- gene_info[-grep("CHR", gene_info$chromosome_name),]
  }
  gene_info[is.na(gene_info)] <- "-"
  gene_info[gene_info == ""] <- "-"
  gene_info <- gene_info[,c(1:3)]
  return(gene_info)
 }


gene_variant_distance_finder <- function(lead_variant_ranges, local_gene_ranges) {
  distance <- vector(mode = "numeric", length = length(local_gene_ranges))
  rank <- vector(mode = "numeric", length = length(local_gene_ranges))
  for(i in 1:length(lead_variant_ranges)) {
    indices <- which(mcols(local_gene_ranges)[[1]] == names(lead_variant_ranges)[i])
    if(length(indices) > 0) {
      distance[indices] <- distance(lead_variant_ranges[i],
                                    local_genes[which(mcols(local_gene_ranges)[[1]] == names(lead_variant_ranges)[i])])
      temp <- vector(mode = "numeric", length = length(indices))
      temp[order(distance[indices])] <- 1:length(indices)
      rank[indices] <- temp
    }
  }
  mcols(local_gene_ranges)[2] <- distance
  names(mcols(local_gene_ranges))[2] <- "distance_to_lead"
  mcols(local_gene_ranges)[3] <- rank
  names(mcols(local_gene_ranges))[3] <- "distance_ranking"
  return(local_gene_ranges)
}


nearest_genes_selector <- function(local_gene_ranges, annotation, biotype, number) {
  #nearest_genes <- NULL
  nearest_genes <- data.frame(
      LEAD_rsID = character(),
      distance_to_lead = integer(),
      distance_ranking = integer(),
      ensembl_id = character(),
      hgnc_symbol = character(),
      gene_biotype = character(),
      stringsAsFactors = FALSE
    )
  
  annotation <- annotation[which(annotation$gene_biotype %in% biotype),]
  unique_lead <- unique(mcols(local_gene_ranges)[[1]])
  for(i in 1:length(unique_lead)) {
    temp <- local_gene_ranges[which(mcols(local_gene_ranges)[[1]] == unique_lead[i])]
    temp <- temp[which(names(temp) %in% annotation$ensembl_gene_id)]
    temp_length <- length(temp)
    if(temp_length > 0) {
      for(x in 1:min(number, temp_length)) {
        index <- which(mcols(temp)[[3]] == min(mcols(temp)[[3]]))
        nearest_genes <- rbind(nearest_genes, c(as.character(mcols(temp[index])[[1]]), 
                                                as.character(mcols(temp[index])[[2]]), x, 
                                                as.character(annotation[which(annotation$ensembl_gene_id == names(temp)[index]),])))
        temp <- temp[-index]
      }
    }
  }

  colnames(nearest_genes) <- c("LEAD_rsID", "distance_to_lead", "distance_ranking", "ensembl_id",
                               "hgnc_symbol", "gene_biotype")
  return(data.table(nearest_genes))
}



cis_eQTL_target_finder <- function(dir, tissues, lead, proxies) {
  
  lead_eQTL_hits <- NULL
  proxy_eQTL_hits <- NULL
  for(i in tissues) {
    # current_tissue <- unlist(strsplit(i[1], ".signifpairs.txt"))
    current_tissue <- unlist(strsplit(i[1], ".signif_variant_gene_pairs.txt.gz"))
    temp_tissue_data <- data.table(read.table(file = file.path(dir, i), header = TRUE, quote = NULL,
                                              sep = "\t", stringsAsFactors = FALSE))
    #temp_search_terms <- gsub("^([^_]*_[^_]*_).*$", "\\1", temp_tissue_data$variant_id)
    temp_search_terms <- gsub("^chr([^_]*_)([^_]*_).*$", "\\1\\2", temp_tissue_data$variant_id)
    lead_indices <- which(temp_search_terms %in% lead$GTEx_search_term)
    tissue <- rep(current_tissue, length(lead_indices))
    search_term <- temp_search_terms[lead_indices]
    lead_eQTL_hits <- rbind(lead_eQTL_hits,
                                cbind(search_term, 
                                      temp_tissue_data[lead_indices, c(1:4, 9)],
                                      tissue))
    
    proxy_indices <- which(temp_search_terms %in% proxies$GTEx_search_term)
    tissue <- rep(current_tissue, length(proxy_indices))
    search_term <- temp_search_terms[proxy_indices]
    
    proxy_eQTL_hits <- rbind(proxy_eQTL_hits,
                             cbind(search_term,
                                   temp_tissue_data[proxy_indices, c(1:4, 9)], 
                                   tissue))
    
    names(lead_eQTL_hits) <- c("search_term", names(temp_tissue_data)[c(1:4, 9)], "tissue")
    names(proxy_eQTL_hits) <- c("search_term", names(temp_tissue_data)[c(1:4, 9)], "tissue")
    
    cat(paste("\t\tSearch of ", current_tissue, " complete!\n", sep = ""))
  }
  return(list(lead_eQTL_hits, proxy_eQTL_hits))
}

   

cis_eQTL_formatter <- function(lead, proxies, cis_eQTLs) {
  cis_eQTLs[[1]]$gene_id <- gsub("[.].*$", "", cis_eQTLs[[1]]$gene_id)
  cis_eQTLs[[2]]$gene_id <- gsub("[.].*$", "", cis_eQTLs[[2]]$gene_id)
  
  LEAD_rsID <- vector(mode = "character", length = length(cis_eQTLs[[1]]$search_term))
  for(i in 1:length(lead$rsID)) {
    indices <- which(cis_eQTLs[[1]]$search_term == lead$GTEx_search_term[i])
    LEAD_rsID[indices] <- lead$rsID[i]
  }
  
  cis_eQTLs[[1]] <- cbind(LEAD_rsID, cis_eQTLs[[1]])
  
  PROXY_rsID <- vector(mode = "character", length = length(cis_eQTLs[[2]]$search_term))
  for(x in 1:length(proxies$PROXY_rsID)) {
    indices <- which(cis_eQTLs[[2]]$search_term == proxies$GTEx_search_term[x])
    PROXY_rsID[indices] <- proxies$PROXY_rsID[x]
  }
  
  cis_eQTLs[[2]] <- cbind(PROXY_rsID, cis_eQTLs[[2]])

  return(cis_eQTLs)
}



cis_eQTL_gene_extractor <- function(lead, proxies, cis_eQTLs) {
  cis_eQTL_targets <- NULL
  for(i in 1:length(lead$rsID)) {
    
    lead_genes <- NULL
    if(lead$rsID[i] %in% cis_eQTLs[[1]]$LEAD_rsID) {
      lead_indices <- which(cis_eQTLs[[1]]$LEAD_rsID == lead$rsID[i])
      lead_genes <- unique(cis_eQTLs[[1]]$gene_id[lead_indices])
    }
    proxy_variants <- proxies$PROXY_rsID[which(proxies$LEAD_rsID == lead$rsID[i])]
    
    proxy_genes <- unique(cis_eQTLs[[2]]$gene_id[cis_eQTLs[[2]]$PROXY_rsID %in% proxy_variants])
    if(length(lead_genes) > 0 | length(proxy_genes) > 0) {
      all_genes <- unique(unlist(c(lead_genes, proxy_genes)))
      lead_proxy <- vector(mode = "character", length = length(all_genes))
      for(x in 1:length(all_genes)) {
        if(all_genes[x] %in% lead_genes & all_genes[x] %in% proxy_genes) {
          lead_proxy[x] <- "lead_and_proxy"
        }
        else if(all_genes[x] %in% proxy_genes) {
          lead_proxy[x] <- "proxy"
        }
        else { lead_proxy[x] <- "lead" }
      }
      cis_eQTL_targets <- rbind(cis_eQTL_targets, 
                                cbind(rep(lead$rsID[i], length(all_genes)),
                                      all_genes, lead_proxy))
      colnames(cis_eQTL_targets) <- c("rsID", "ensembl_id", "lead_or_proxy")
    }
  }
  return(data.table(cis_eQTL_targets))
}


anno_summariser <- function(lead, nearest, LD_overlapping, LD_annotation, cis_eQTLs,
                                 cis_eQTL_annotation) {
  anno_summary <- NULL
  for(i in 1:length(lead$rsID)) {
    #print(lead$rsID[i])
    temp_nearest <- nearest$ensembl_id[which(nearest$LEAD_rsID == lead$rsID[i])]
    if (length(temp_nearest) == 0) {
      temp_nearest <- NA
    }
    
    temp_LD_overlapping <- names(LD_overlapping[which(mcols(LD_overlapping)[[1]] == lead$rsID[i])])
    temp_cis_eQTLs <- cis_eQTLs$ensembl_id[which(cis_eQTLs$rsID == lead$rsID[i])]
    
 
    lead_or_proxy <- cis_eQTLs$lead_or_proxy[which(cis_eQTLs$rsID == lead$rsID[i])]
    
    all_genes <- unique(c(temp_nearest, temp_LD_overlapping, temp_cis_eQTLs))
    hgnc <- vector(mode = "character", length = length(all_genes))
    biotype <- vector(mode = "character", length = length(all_genes))
    if(length(all_genes)>0){
    for(z in 1:length(all_genes)) {

      if(!is.na(all_genes[z]) && all_genes[z] %in% nearest$ensembl_id) {
        hgnc[z] <- nearest$hgnc_symbol[match(all_genes[z], nearest$ensembl_id)]
        biotype[z] <- nearest$gene_biotype[match(all_genes[z], nearest$ensembl_id)]
      }
      else if(all_genes[z] %in% LD_annotation$ensembl_gene_id) {
        hgnc[z] <- LD_annotation$hgnc_symbol[match(all_genes[z], LD_annotation$ensembl_gene_id)]
        biotype[z] <- LD_annotation$gene_biotype[match(all_genes[z], LD_annotation$ensembl_gene_id)]
      }
      else if(all_genes[z] %in% cis_eQTL_annotation$ensembl_gene_id) {
        hgnc[z] <- cis_eQTL_annotation$hgnc_symbol[match(all_genes[z], cis_eQTL_annotation$ensembl_gene_id)]
        biotype[z] <- cis_eQTL_annotation$gene_biotype[match(all_genes[z], cis_eQTL_annotation$ensembl_gene_id)]
      } else {
        hgnc[z] <- NA
        biotype[z] <- NA
      }
    }
    temp_matrix <- matrix(data = 0,
                          nrow = length(all_genes), ncol = 10,
                          dimnames = list(NULL, c("LEAD_rsID", "ensembl_id",
                                                  "hgnc_symbol", "gene_biotype", "nearest",
                                                  "second_nearest", "third_nearest",
                                                  "LD_overlapping", "lead_eQTL", "proxy_eQTL")))
    for(x in 1:length(all_genes)) {
      if(!is.na(all_genes[x]) && all_genes[x] %in% temp_nearest) {
        index <- which(temp_nearest == all_genes[x])
        if(index == 1) {
          temp_matrix[x, 5] <- 1
        }
        else if(index == 2) {
          temp_matrix[x, 6] <- 1
        }
        else if(index == 3) {
          temp_matrix[x, 7] <- 1
        }
     # } else if (is.na(temp_nearest)) {
      } else if (all(is.na(temp_nearest))) {
        temp_matrix[x, 5] <- NA
      }
      if(all_genes[x] %in% temp_LD_overlapping) {
        temp_matrix[x, 8] <- 1
      }
      if(all_genes[x] %in% temp_cis_eQTLs) {
        index <- which(temp_cis_eQTLs == all_genes[x])
        if(lead_or_proxy[index] == "lead") {
          temp_matrix[x, 9] <- 1
        }
        else if(lead_or_proxy[index] == "proxy") {
          temp_matrix[x, 10] <- 1
        }
        else { temp_matrix[x, c(9, 10)] <- 1 }
      }
    }
    temp_matrix[,1] <- rep(lead$rsID[i], length(all_genes))
    temp_matrix[,2] <- all_genes
    temp_matrix[,3] <- hgnc
    temp_matrix[,4] <- biotype
    anno_summary <- rbind(anno_summary, temp_matrix)
  }
  }
  return(data.table(anno_summary))
}


VEP_integrator <- function(lead, anno_summary, VEP, lead_col, proxy_col, ensembl_col, IMPACT_col) {
  VEP_matrix <- matrix(data = 0, nrow = length(anno_summary$LEAD_rsID), ncol = 2,
                       dimnames = list(NULL, c("lead_IMPACT", "proxy_IMPACT")))
  for(i in 1:length(lead$rsID)) {
    l <- lead$rsID[i]

    if(l %in% VEP[[lead_col]]) {
      temp <- VEP[which(VEP[[lead_col]] == l),]
      for(x in unique(temp[[ensembl_col]])) {
        temp1 <- temp[which(temp[[ensembl_col]] == x),]
        if("MODERATE" %in% temp1[[IMPACT_col]]) {
          temp2 <- temp1[which(temp1[[IMPACT_col]] == "MODERATE"),]
          boolean <- temp2[[lead_col]] == temp2[[proxy_col]]
          if(TRUE %in% boolean) {
            VEP_matrix[which(anno_summary$LEAD_rsID == l & anno_summary$ensembl_id == x), 1] <- 1
          }
          if(FALSE %in% boolean) {
            VEP_matrix[which(anno_summary$LEAD_rsID == l & anno_summary$ensembl_id == x), 2] <- 1
          }
        }
        if("HIGH" %in% temp1[[IMPACT_col]]) {
          temp3 <- temp1[which(temp1[[IMPACT_col]] == "HIGH"),]
          boolean <- temp3[[lead_col]] == temp3[[proxy_col]]
          if(TRUE %in% boolean) { 
              VEP_matrix[which(anno_summary$LEAD_rsID == l & anno_summary$ensembl_id == x), 1] <- 2
            }
          if(FALSE %in% boolean) {
            VEP_matrix[which(anno_summary$LEAD_rsID == l & anno_summary$ensembl_id == x), 2] <- 2
            }
          }
        }
      }
    }
  return(data.table(VEP_matrix))
}

COLOC_integrator <- function(lead, anno_summary, COLOC, lead_col, ensembl_col, QTL_type = "eQTL") {
  col.n <- paste0("coloc_", QTL_type)
  COLOC_matrix <- matrix(data = 0, nrow = length(anno_summary$LEAD_rsID), ncol = 1,
                         dimnames = list(NULL, c(col.n)))
  for (i in 1:length(lead$rsID)) {
    l <- lead$rsID[i]
    if (l %in% COLOC[[lead_col]]) {
      temp <- COLOC[which(COLOC[[lead_col]] == l), ]
      for(x in unique(temp[[ensembl_col]])) {
        COLOC_matrix[which(anno_summary$LEAD_rsID == l & anno_summary$ensembl_id == x), 1] <- 1
      }
    }
  }
  return(data.table(COLOC_matrix))
}

