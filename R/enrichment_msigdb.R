#' Run MSigDB Gene Set Enrichment Analysis
#'
#' @param ranked_symbols Named numeric vector (gene symbols as names)
#' @param collection MSigDB collection category
#' @param subcollection MSigDB subcollection (optional)
#' @param strip_prefix Remove common prefixes from gene set names?
#' @param species Species name for msigdbr
#' @param ... Additional GSEA parameters
#'
#' @return gseaResult object
#' @export
run_gsea_msigdb <- function(
  ranked_symbols,
  collection = "H",
  subcollection = NULL,
  strip_prefix = TRUE,
  species = "Homo sapiens",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 1e-10,
  nPermSimple = 10000,
  verbose = FALSE
) {
  # Input validation
  if (!is.numeric(ranked_symbols)) {
    stop("ranked_symbols must be a numeric vector")
  }
  if (is.null(names(ranked_symbols))) {
    stop("ranked_symbols must be a named vector (gene SYMBOLs as names)")
  }
  if (any(is.na(ranked_symbols)) || any(!is.finite(ranked_symbols))) {
    warning("Removing non-finite values from ranked_symbols")
    ranked_symbols <- ranked_symbols[is.finite(ranked_symbols)]
  }

  # Check for duplicates
  if (any(duplicated(names(ranked_symbols)))) {
    warning("Duplicate gene symbols detected. Keeping first occurrence only.")
    ranked_symbols <- ranked_symbols[!duplicated(names(ranked_symbols))]
  }

  # Fetch MSigDB gene sets
  if (verbose) {
    msg <- paste0("Fetching MSigDB gene sets (", collection)
    if (!is.null(subcollection)) msg <- paste0(msg, ":", subcollection)
    msg <- paste0(msg, ")...")
    message(msg)
  }

  msig_sets <- msigdbr::msigdbr(
    species = species,
    collection = collection,
    subcollection = subcollection
  )

  # Optionally strip common prefixes for cleaner term names
  if (strip_prefix) {
    # Common prefixes by collection
    prefix_map <- c(
      "H" = "HALLMARK_",
      "C2" = "^(KEGG_MEDICUS_|REACTOME_|BIOCARTA_|PID_|WIKIPATHWAYS_)",
      "C5" = "^(GO_|GOBP_|GOCC_|GOMF_)",
      "C6" = "^",
      "C7" = "^",
      "C8" = "^"
    )

    if (collection %in% names(prefix_map)) {
      msig_sets <- msig_sets |>
        dplyr::mutate(gs_name = stringr::str_replace_all(gs_name, prefix_map[[collection]], ""))
    }
  }

  # Build TERM2GENE data frame
  term2gene <- msig_sets |>
    dplyr::select(term = gs_name, gene = gene_symbol) |>
    dplyr::distinct() |>
    as.data.frame()

  if (verbose) {
    message("Running GSEA on ", length(unique(term2gene$term)), " gene sets...")
  }

  # Run GSEA
  gsea_result <- clusterProfiler::GSEA(
    geneList     = ranked_symbols,
    TERM2GENE    = term2gene,
    minGSSize    = minGSSize,
    maxGSSize    = maxGSSize,
    eps          = eps,
    nPermSimple  = nPermSimple,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "BH",
    verbose      = verbose
  )

  if (verbose) {
    n_sig <- sum(gsea_result@result$p.adjust < 0.05)
    message("Completed: ", n_sig, " significant gene sets (padj < 0.05)")
  }

  return(gsea_result)
}


#' Run Multiple MSigDB Collections
#'
#' @param ranked_symbols Named numeric vector
#' @param collections Vector of collection names
#' @param kegg_variant Which KEGG variant to use
#' @param species Species name
#' @param parallel Use parallel processing?
#' @param n_workers Number of parallel workers
#' @param ... Additional parameters
#'
#' @return Named list of gseaResult objects
#' @export
run_gsea_msigdb_suite <- function(
  ranked_symbols,
  collections = c("hallmarks", "reactome"),
  kegg_variant = "KEGG_MEDICUS",
  species = "Homo sapiens",
  parallel = FALSE,
  n_workers = NULL,
  verbose = FALSE,
  ...
) {
  
  # Define collection specifications
  collection_map <- list(
    hallmarks     = list(cat = "H",  subcat = NULL),
    kegg          = list(cat = "C2", subcat = kegg_variant),
    reactome      = list(cat = "C2", subcat = "REACTOME"),
    biocarta      = list(cat = "C2", subcat = "BIOCARTA"),
    wikipathways  = list(cat = "C2", subcat = "WIKIPATHWAYS"),
    pid           = list(cat = "C2", subcat = "PID"),
    go_bp         = list(cat = "C5", subcat = "GO:BP"),
    go_mf         = list(cat = "C5", subcat = "GO:MF"),
    go_cc         = list(cat = "C5", subcat = "GO:CC"),
    oncogenic     = list(cat = "C6", subcat = NULL),
    immunologic   = list(cat = "C7", subcat = "IMMUNESIGDB"),
    cell_type     = list(cat = "C8", subcat = NULL)
  )
  
  # Validate requested collections
  available <- names(collection_map)
  collections <- match.arg(collections, available, several.ok = TRUE)
  
  # Get specs for requested collections
  specs_to_run <- collection_map[collections]
  
  # Setup parallel processing if requested
  if (parallel) {
    if (!requireNamespace("furrr", quietly = TRUE)) {
      stop("Package 'furrr' is required for parallel processing. Install with: install.packages('furrr')")
    }
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required for parallel processing. Install with: install.packages('future')")
    }
    
    # Set up workers
    if (is.null(n_workers)) {
      n_workers <- future::availableCores() - 1
    }
    
    if (verbose) message("Setting up parallel processing with ", n_workers, " workers...")
    future::plan(future::multisession, workers = n_workers)
    
    map_fn <- furrr::future_map
  } else {
    map_fn <- purrr::map
  }
  
  # Run GSEA for each collection
  results <- map_fn(
    names(specs_to_run),
    function(coll_name) {
      spec <- specs_to_run[[coll_name]]
      
      if (verbose) {
        msg <- paste0("Running ", coll_name, " (", spec$cat)
        if (!is.null(spec$subcat)) msg <- paste0(msg, ":", spec$subcat)
        msg <- paste0(msg, ")...")
        message(msg)
      }
      
      run_gsea_msigdb(
        ranked_symbols = ranked_symbols,
        collection = spec$cat,
        subcollection = spec$subcat,
        species = species,
        verbose = verbose,
        ...
      )
    },
    .progress = if (parallel && verbose) TRUE else FALSE
  )
  
  names(results) <- names(specs_to_run)
  
  # Clean up parallel workers
  if (parallel) {
    future::plan(future::sequential)
  }
  
  return(results)
}