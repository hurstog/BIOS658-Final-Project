
#### Perform Functional Analysis ####

#Create name of the output file and set settings
data_dir <- getwd()
p_val_cutoff   <- 0.05 
FDR_cutoff     <- 0.3
msigdb_all     <- FALSE 
degs_sheet <- "proj_deg"
fileNameOut1 <- file.path(data_dir, paste0("GSEA2_", make.names(degs_sheet), ".xlsx"))
run_gsea <- TRUE 
min_kegg_genes <- 20
max_kegg_genes <- 2000 
p_adj_cutoff   <- 1
nperm          <- 1000
num_gseaplots  <- 3 
ntop_symbols   <- 24
human_analysis <- TRUE
mouse_analysis <- FALSE
get_latest_EnsDb <- function(species, genome) {
  ah <- AnnotationHub::AnnotationHub()
  q <- AnnotationHub::query(ah, c("EnsDb", species, genome))
  if (length(q) == 0L) stop("No EnsDb found for: ", species, " ", genome)
  q[[which.max(as.Date(mcols(q)$rdatadateadded))]]
}
build_gene_data <- function(ensdb, orgdb, canonical_chrs, species_code, kegg_set, msigdbr_org) {
  g <- genes(
    ensdb,
    columns = c("gene_id", "gene_name", "gene_biotype", "seq_name",
                "gene_seq_start", "gene_seq_end")
  )
  df <- as_tibble(g) %>%
    transmute(
      ensgene = gene_id,
      symbol = gene_name,
      biotype = gene_biotype,
      chr = as.character(seqnames),
      start = as.integer(start),
      end = as.integer(end)
    ) %>%
    dplyr::filter(chr %in% canonical_chrs)
  
  # Map description (GENENAME) and ENTREZID from org.*.eg.db using ENSEMBL keys
  desc_map <- AnnotationDbi::select(
    orgdb,
    keys = unique(df$ensgene),
    keytype = "ENSEMBL",
    columns = c("GENENAME", "ENTREZID")
  ) %>%
    dplyr::rename(ensgene = ENSEMBL, description = GENENAME, entrezid = ENTREZID)
  
  df <- df %>%
    left_join(desc_map, by = "ensgene")
  
  # Final gene_annotations: protein_coding, non-NA symbol/description
  gene_annotations <- df %>%
    # dplyr::filter(
    # biotype == "protein_coding",
    # !is.na(symbol),
    # !is.na(description),
    # description != ""
    # ) %>%
    distinct(ensgene, symbol, biotype, description, entrezid, .keep_all = FALSE)
  
  # Gene length as gene span (end - start)
  gene_length <- df %>%
    transmute(Geneid = ensgene, Length = end - start) %>%
    distinct()
  
  # Background symbols
  all.symbol <- unique(gene_annotations$symbol)
  
  list(
    gene_annotations = gene_annotations,
    gene_length = gene_length,
    all.symbol = all.symbol,
    OrgDb = orgdb,
    species = species_code,
    KEGG = kegg_set,
    msigdbr_org = msigdbr_org
  )
}

if (exists("human_analysis") && isTRUE(human_analysis)) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
  
  edb_hs <- get_latest_EnsDb("Homo sapiens", "GRCh38")
  
  human <- build_gene_data(
    ensdb = edb_hs,
    orgdb = org.Hs.eg.db,
    canonical_chrs = c(as.character(1:22), "X", "Y", "MT"),
    species_code = "hsa",
    kegg_set = "KEGG_2019_Human",
    msigdbr_org = "Homo sapiens"
  )
  
  gene_annotations <- human$gene_annotations
  gene_length <- human$gene_length
  all.symbol <- human$all.symbol
  OrgDb <- "org.Hs.eg.db"
  species <- "hsa"
  KEGG <- human$KEGG
  msigdbr_org <- human$msigdbr_org
}

#Load data and prepare for analysis
mtx <- read.csv("DGE_results.csv")
#Remove genes that do not have a Gene Symbol
mtx <- mtx[-c(which(is.na(mtx$Gene_Symbol))),]
#Order the genes from most to least significant
mtx <- mtx[order(mtx$FDR), ] 
#Remove repeated genes (keep first instance of the gene)
mtx <- mtx %>%
  distinct(Gene_Symbol, .keep_all = TRUE)
res <- data.frame(symbol = mtx$Gene_Symbol, logFC = mtx$logFC, p.val = mtx$PValue, p.adj = mtx$FDR)
res <- res[ order(res$p.val, decreasing = FALSE), ]

##### Perform Hypergeometric analysis #####
#For Kegg Pathways:
kegg_enrich <- function(compartment_genes, p_adj_cutoff = p_adj_cutoff) {
  # stop if min_kegg is not met 
  stopifnot(length(compartment_genes) >= min_kegg_genes)
  res.kegg <- enrichr(unique(compartment_genes), databases = KEGG) # KEGG results only
  # If significant results are present, save them
  if (nrow(res.kegg[[KEGG]]) > 0 & sum(res.kegg[[KEGG]]$Adjusted.P.value < p_adj_cutoff) > 0) {
    res.kegg <- as.data.frame(res.kegg[[KEGG]])
    res.kegg <- res.kegg[res.kegg$Adjusted.P.value < p_adj_cutoff, , drop = FALSE]
    compartment_genes <- res.kegg
    # reorder the genes alphabetically 
    compartment_genes <- compartment_genes %>% 
      # separate the rows by splitting by the delimiter and expanding the rows 
      separate_rows(Genes, convert = TRUE, sep = ";") %>% 
      # group the genes by their Term & other columns to keep them 
      group_by(Term, Overlap, P.value, Adjusted.P.value, Old.P.value, Old.Adjusted.P.value,
               Odds.Ratio, Combined.Score) %>% 
      # sort the genes for each term 
      arrange(Genes) %>% 
      summarise(Genes = paste(Genes, collapse="/")) %>% 
      arrange(P.value, Adjusted.P.value)
  } else {
    compartment_genes <- as.data.frame(matrix(data = "Nothing significant", nrow = 1, ncol = 9))
    colnames(compartment_genes) <- c("Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")
    compartment_genes$`P.value` = 0
    compartment_genes$`Adjusted.P.value` = 0
    compartment_genes$`Old.P.value` = 0
    compartment_genes$`Old.Adjusted.P.value` = 0
    compartment_genes$`Odds.Ratio` = 0
    compartment_genes$`Combined.Score` = 0
  }
  return(compartment_genes)
}

# run unranked KEGG analysis 
websiteLive <- TRUE # Check if EnrichR is up
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if(websiteLive) {
  # Subset the number of DEGs for KEGG analysis to the maximum
  if (nrow(res[res$p.val < p_val_cutoff, ]) > max_kegg_genes) {
    degs_subset <- res[1:max_kegg_genes, ]
  } else {
    degs_subset <- res[res$p.adj < FDR_cutoff, ]
  }
  # Get list of up- and downregulated genes
  up.genes <- sort(unique(degs_subset$symbol[ degs_subset$logFC > 0 ]))
  dn.genes <- sort(unique(degs_subset$symbol[ degs_subset$logFC < 0 ]))
  res.kegg  <- NULL # Initially, empty value
  print(paste0("KEGG pathway run on ", length(unique(c(up.genes, dn.genes))),
               " genes without distinguishing them by directionality."))
  # run hypergeomtric KEGG analysis on combination of genes 
  res.kegg <- kegg_enrich(unique(c(up.genes, dn.genes)), 
                          p_adj_cutoff = p_adj_cutoff)
}

#Now for the three gene ontologies
entrez_to_symbols <- function(entrez_list, gene_annotations, ntop = 24) {
  # Create a fast lookup table (named vector) from entrezid to symbol
  # Filter out NA entrezids, keep first occurrence per entrezid for O(1) hash-based lookups
  entrez_map_df <- gene_annotations %>%
    dplyr::filter(!is.na(entrezid)) %>%
    dplyr::distinct(entrezid, symbol)
  entrez_lookup <- setNames(entrez_map_df$symbol, entrez_map_df$entrezid)
  
  # Function to convert a single slash-separated string
  convert_single <- function(entrez_string) {
    entrez_ids <- unlist(strsplit(entrez_string, "/"))
    # Use named vector indexing for O(1) lookup with proper NA handling
    matched_symbols <- entrez_lookup[entrez_ids]
    symbols <- ifelse(is.na(matched_symbols), entrez_ids, matched_symbols)
    # Sort alphabetically for consistent output across runs
    paste(sort(symbols), collapse = "/")
  }
  
  # Initialize result as original input
  result <- entrez_list
  
  # Process top ntop entries (if ntop is Inf, process all entries)
  n_to_process <- if (is.infinite(ntop)) length(entrez_list) else min(ntop, length(entrez_list))
  result[1:n_to_process] <- vapply(entrez_list[1:n_to_process], convert_single, character(1))
  
  return(result)
}
if (msigdb_all) {
  m_df <- msigdbr(species = msigdbr_org)
  m_df_gs_cat <- unique(m_df$gs_collection) %>% sort()
} else {
  # Use most informative ones
  m_df_gs_cat <- c("C2", "C5", "H")
}

# Function to perform enrichment analysis using MSigDb signatures
msigdb_enrich <- function(dataset, p_adj_cutoff = p_adj_cutoff) {
  # Top DEGs for enrichr
  # res.all <- dataset[dataset$p.val < p_val_cutoff, ]
  # Subset the number of DEGs for KEGG analysis to the maximum
  if (nrow(dataset[dataset$p.val < p_val_cutoff, ]) > max_kegg_genes) {
    res.all <- dataset[1:max_kegg_genes, ]
  } else {
    res.all <- dataset[dataset$p.adj < p_val_cutoff, ]
  }
  # Convert symbols to entrezids
  eid <- bitr(res.all$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb=OrgDb)
  # Attach converted entrezids
  res.all <- left_join(res.all, eid, by = c("symbol" = "SYMBOL"))
  res.all <- res.all[ !is.na(res.all$ENTREZID), ]
  # List of t-statistics
  geneList_significant <- res.all$logFC
  # Make it named
  names(geneList_significant) <- res.all$ENTREZID
  # And decreasing sorted
  geneList_significant <- sort(geneList_significant, decreasing = TRUE)
  res.msigdf.all <- list()
  # gs_cat="C2" # For testing
  for (gs_cat in m_df_gs_cat) {
    # Term to gene
    m_t2g <- msigdbr(species = msigdbr_org, collection = gs_cat) %>% 
      dplyr::distinct(gs_name, ncbi_gene)
    # Term to description 
    m_t2d <- msigdbr(species = msigdbr_org, collection = gs_cat) %>% 
      dplyr::distinct(gs_name, gs_description)
    # Enrichment analysis
      em <- enricher(names(geneList_significant), TERM2GENE=m_t2g, pvalueCutoff = p_adj_cutoff)
      if (!is.null(em) | !any(is.na(em@result$qvalue))) {
        em@result <- em@result %>% dplyr::mutate(Direction = "ALL") # add new Direction column 
      }
      #  # Check if the results are non-empty
      if (!(is.null(em))) {
        res.msigdf.em <- rbind(em@result)
        res.msigdf.em <- res.msigdf.em[res.msigdf.em$p.adjust < p_adj_cutoff, , drop = FALSE]
        res.msigdf.em$core_enrichment <- entrez_to_symbols(res.msigdf.em$geneID, gene_annotations, ntop = ntop_symbols)
        # sort the genes alphabetically
        res.msigdf.em <- res.msigdf.em %>% 
          # separate the rows by splitting by the delimiter and expanding the rows 
          separate_rows(geneID, convert = TRUE, sep = "/") %>% 
          # group the genes by their other columns to keep them 
          group_by(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count, Direction, core_enrichment) %>%
          # sort the genes for each term 
          arrange(geneID) %>% 
          summarise(geneID = paste(geneID, collapse="/")) %>% 
          arrange(pvalue, p.adjust)
        # Append description
        res.msigdf.em <- left_join(res.msigdf.em, m_t2d, by = c("ID" = "gs_name"))
        res.msigdf.em$Description <- res.msigdf.em$gs_description
        res.msigdf.em$gs_description <- NULL
      } else {
        res.msigdf.em <- as.data.frame(matrix(data = "Nothing significant", nrow = 1,
                                              ncol = 10))
        colnames(res.msigdf.em) <- c("ID", "Description", "GeneRatio", "BgRatio",
                                     "pvalue", "p.adjust", "qvalue", "geneID", "Count", "Direction")
      }
    # Combine the results and add names
    res.msigdf.all <- c(res.msigdf.all, list(res.msigdf.em))
    names(res.msigdf.all)[length(res.msigdf.all)] <- paste0("Enrich.", gs_cat)
  }
  res.msigdf.all
}
# run unranked MSigDb analysis 
res.msigdb.enrich <- msigdb_enrich(dataset = res, p_adj_cutoff = p_adj_cutoff)

##### Now run GSEA #####
#First for Kegg pathways
gsea_kegg_enrich <- function(dataset = res, p_adj_cutoff = 0.05) {
  res.all <- dataset %>%
    dplyr::group_by(symbol) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  eid <- bitr(unique(res.all$symbol), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  res.all <- dplyr::left_join(res.all, eid, by = c("symbol" = "SYMBOL"))
  res.all <- res.all[!is.na(res.all$ENTREZID), ]
  geneList <- res.all$logFC
  names(geneList) <- res.all$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  set.seed(1)
  ego3 <- gseKEGG(
    geneList     = geneList,
    organism     = species,
    minGSSize    = 10,
    pvalueCutoff = 1,
    verbose      = FALSE
  )
  ego3 <- setReadable(ego3, OrgDb = OrgDb, keyType = "ENTREZID")
  return(ego3)
}

res.kegg.ego3 <- gsea_kegg_enrich(res, p_adj_cutoff = p_adj_cutoff)

# Function to reformat gseaResult object into processed data frame
gsea_kegg_results_to_df <- function(ego3, p_adj_cutoff = 0.05) {
  res.kegg.gsea <- as.data.frame(ego3)
  if (nrow(res.kegg.gsea) > 0) {
    res.kegg.gsea <- res.kegg.gsea[, c("ID", "Description", "NES", "pvalue", "p.adjust", "core_enrichment")]
    res.kegg.gsea <- res.kegg.gsea %>%
      tidyr::separate_rows(core_enrichment, convert = TRUE, sep = "/") %>%
      dplyr::group_by(ID, Description, NES, pvalue, p.adjust) %>%
      dplyr::arrange(core_enrichment) %>%
      dplyr::summarise(core_enrichment = paste(core_enrichment, collapse = "/"), .groups = "drop")
    res.kegg.gsea <- res.kegg.gsea[order(abs(res.kegg.gsea$NES), decreasing = TRUE), ]
    res.kegg.gsea <- res.kegg.gsea[res.kegg.gsea$p.adjust < p_adj_cutoff, ]
    res.kegg.gsea$NES       <- round(res.kegg.gsea$NES, digits = 2)
    res.kegg.gsea$pvalue    <- formatC(res.kegg.gsea$pvalue, format = "e", digits = 2)
    res.kegg.gsea$p.adjust  <- formatC(res.kegg.gsea$p.adjust, format = "e", digits = 2)
    rownames(res.kegg.gsea) <- NULL
  } else {
    res.kegg.gsea <- as.data.frame(matrix(data = "Nothing significant", nrow = 1, ncol = 6))
    colnames(res.kegg.gsea) <- c("ID", "Description", "NES", "pvalue", "p.adjust", "core_enrichment")
  }
  return(res.kegg.gsea)
}

res.kegg.gsea <- gsea_kegg_results_to_df(res.kegg.ego3, p_adj_cutoff = p_adj_cutoff)

# Plot the top enriched pathways
for (i in 1:num_gseaplots) {
  pathway_id <- res.kegg.gsea$ID[i]
  pathway_desc <- res.kegg.gsea$Description[i]
  print(gseaplot2(res.kegg.ego3, geneSetID = pathway_id, title = pathway_desc))
}

#Now for gene ontologies 
gsea_msigdb_enrich <- function(dataset, p_adj_cutoff = 0.05) {
  # Preprocess dataset: get the most significant gene per symbol
  res.all <- dataset %>%
    dplyr::group_by(symbol) %>%
    dplyr::slice_max(order_by = abs(logFC), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
  # Convert gene symbols to entrez IDs
  eid <- bitr(res.all$symbol, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  res.all <- dplyr::left_join(res.all, eid, by = c("symbol" = "SYMBOL"))
  res.all <- res.all[!is.na(res.all$ENTREZID), ]
  # Named gene list
  geneList <- res.all$logFC
  names(geneList) <- res.all$ENTREZID
  geneList <- sort(geneList, decreasing = TRUE)
  # GSEA for each MSigDB category
  res.msigdb.all <- list()
  for (gs_cat in m_df_gs_cat) {
    m_t2g <- msigdbr(species = msigdbr_org, collection = gs_cat) %>%
      dplyr::distinct(gs_name, ncbi_gene)
    em2 <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = p_adj_cutoff)
    res.msigdb.all[[paste0("GSEA.", gs_cat)]] <- em2
  }
  return(res.msigdb.all)
}

res.msigdb.em2 <- gsea_msigdb_enrich(res, p_adj_cutoff = p_adj_cutoff)


# Function to reformat a list of GSEA results into a list of data frames
gsea_msigdb_results_to_df <- function(res.msigdb.all, gene_annotations, msigdbr_org) {
  res.msigdf.all <- list()
  for (nm in names(res.msigdb.all)) {
    em2 <- res.msigdb.all[[nm]]
    # Term to description 
    res.msigdf.em2 <- NULL
    if (nrow(em2@result) > 0) {
      gs_cat <- sub("GSEA.", "", nm)
      m_t2d <- msigdbr(species = msigdbr_org, collection = gs_cat) %>%
        dplyr::distinct(gs_name, gs_description)
      res.msigdf.em2 <- em2@result
      res.msigdf.em2$core_enrichment <- entrez_to_symbols(res.msigdf.em2$core_enrichment, gene_annotations, ntop = ntop_symbols)
      res.msigdf.em2 <- res.msigdf.em2 %>%
        tidyr::separate_rows(core_enrichment, convert = TRUE, sep = "/") %>%
        dplyr::group_by(ID, Description, setSize, enrichmentScore, NES, pvalue, p.adjust, qvalue,
                        rank, leading_edge) %>%
        dplyr::arrange(core_enrichment) %>%
        dplyr::summarise(core_enrichment = paste(core_enrichment, collapse = "/"), .groups = "drop") %>%
        dplyr::arrange(pvalue, p.adjust)
      # Append description
      res.msigdf.em2 <- dplyr::left_join(res.msigdf.em2, m_t2d, by = c("ID" = "gs_name"))
      res.msigdf.em2$Description <- res.msigdf.em2$gs_description
      res.msigdf.em2$gs_description <- NULL
    } else {
      res.msigdf.em2 <- as.data.frame(matrix(data = "Nothing significant", nrow = 1, ncol = 11))
      colnames(res.msigdf.em2) <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues", "rank", "leading_edge", "core_enrichment")
    }
    res.msigdf.all[[nm]] <- res.msigdf.em2
  }
  return(res.msigdf.all)
}

res.msigdb <- gsea_msigdb_results_to_df(res.msigdb.all = res.msigdb.em2, gene_annotations = gene_annotations, msigdbr_org = msigdbr_org)

# Plot the top enriched pathways
for (catname in names(res.msigdb.em2)) {
  gsea_res <- res.msigdb.em2[[catname]]
  gsea_df <- as.data.frame(gsea_res)
  for (i in seq_len(min(num_gseaplots, nrow(gsea_df)))) {
    pathway_id <- gsea_df$ID[i]
    pathway_desc <- gsea_df$Description[i]
    print(
      gseaplot2(
        gsea_res,
        geneSetID = pathway_id,
        title = paste0(catname, ": ", pathway_desc)
      )
    )
  }
}
interleave_lists <- function(list1, list2) {
  n <- length(list1)
  if (length(list2) != n) stop("Lists must have the same length")
  names1 <- names(list1)
  names2 <- names(list2)
  out <- vector("list", 2 * n)
  out_names <- character(2 * n)
  out[seq(1, 2*n, by=2)] <- list1
  out[seq(2, 2*n, by=2)] <- list2
  out_names[seq(1, 2*n, by=2)] <- names1
  out_names[seq(2, 2*n, by=2)] <- names2
  names(out) <- out_names
  out
}

# Example usage:
res.msigdf.all <- interleave_lists(res.msigdb.enrich, res.msigdb)
x <- c(list(Enrich.KEGG = res.kegg), list(GSEA.KEGG = res.kegg.gsea), res.msigdf.all) # 
# names(x)[1:2] <- c("Enrich.KEGG", "GSEA.KEGG")
write_xlsx(x, path = fileNameOut1)
