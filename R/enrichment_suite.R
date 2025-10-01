#' Run Complete GSEA Analysis Suite
#'
#' @param ranked_entrez Named numeric vector (Entrez IDs)
#' @param ranked_symbols Named numeric vector (gene symbols)
#' @param analyses Which analyses to run: "GO", "KEGG", "MSigDB"
#' @param ... Analysis-specific parameters
#'
#' @return gsea_suite_results object (list)
#' @export
run_gsea_suite <- function(
  ranked_entrez = NULL,
  ranked_symbols = NULL,
  analyses = c("GO", "KEGG", "MSigDB"),
  # GO-specific
  go_ont = "BP",
  go_simplify = FALSE,
  go_simplify_cutoff = 0.7,
  go_simplify_measure = "Wang",
  # KEGG-specific
  kegg_organism = "hsa",
  kegg_use_internal = FALSE,
  # MSigDB-specific
  msigdb_collections = "hallmarks",
  msigdb_kegg_variant = "KEGG_MEDICUS",
  msigdb_strip_prefix = TRUE,
  msigdb_species = "Homo sapiens",
  msigdb_parallel = FALSE,
  msigdb_n_workers = NULL,
  # Shared parameters
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 1e-10,
  nPermSimple = 10000,
  verbose = FALSE
) {
  
  # Validate analyses argument
  valid_analyses <- c("GO", "KEGG", "MSigDB")
  analyses <- match.arg(analyses, valid_analyses, several.ok = TRUE)
  
  # Check that required inputs are provided for requested analyses
  needs_entrez <- any(c("GO", "KEGG") %in% analyses)
  needs_symbols <- "MSigDB" %in% analyses
  
  if (needs_entrez && is.null(ranked_entrez)) {
    stop("ranked_entrez is required for GO and/or KEGG analyses")
  }
  if (needs_symbols && is.null(ranked_symbols)) {
    stop("ranked_symbols is required for MSigDB analyses")
  }
  
  # Initialize results list
  results <- list()
  
  # Run GO GSEA
  if ("GO" %in% analyses) {
    if (verbose) message("\n=== Running GO GSEA ===")
    results$GO <- run_gsea_go(
      ranked_entrez = ranked_entrez,
      ont = go_ont,
      simplify = go_simplify,
      simplify_cutoff = go_simplify_cutoff,
      simplify_measure = go_simplify_measure,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      eps = eps,
      nPermSimple = nPermSimple,
      verbose = verbose
    )
  }
  
  # Run KEGG GSEA
  if ("KEGG" %in% analyses) {
    if (verbose) message("\n=== Running KEGG GSEA ===")
    results$KEGG <- run_gsea_kegg(
      ranked_entrez = ranked_entrez,
      organism = kegg_organism,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      eps = eps,
      nPermSimple = nPermSimple,
      use_internal_data = kegg_use_internal,
      verbose = verbose
    )
  }
  
  # Run MSigDB GSEA (now supports multiple collections)
  if ("MSigDB" %in% analyses) {
    if (verbose) message("\n=== Running MSigDB GSEA ===")
    results$MSigDB <- run_gsea_msigdb_suite(
      ranked_symbols = ranked_symbols,
      collections = msigdb_collections,
      kegg_variant = msigdb_kegg_variant,
      species = msigdb_species,
      parallel = msigdb_parallel,
      n_workers = msigdb_n_workers,
      strip_prefix = msigdb_strip_prefix,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pvalueCutoff = pvalueCutoff,
      eps = eps,
      nPermSimple = nPermSimple,
      verbose = verbose
    )
  }
  
  # Add metadata
  attr(results, "timestamp") <- Sys.time()
  attr(results, "analyses") <- analyses
  attr(results, "parameters") <- list(
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff,
    eps = eps,
    nPermSimple = nPermSimple
  )
  
  class(results) <- c("gsea_suite_results", "list")
  
  if (verbose) message("\n=== GSEA Suite Complete ===")
  
  return(results)
}

#' @export
print.gsea_suite_results <- function(x, ...) {
  cat("GSEA Suite Results\n")
  cat("==================\n")
  cat("Analyses run:", paste(attr(x, "analyses"), collapse = ", "), "\n")
  cat("Timestamp:", format(attr(x, "timestamp")), "\n\n")
  
  for (analysis in names(x)) {
    if (analysis == "MSigDB" && is.list(x[[analysis]]) && !inherits(x[[analysis]], "gseaResult")) {
      # MSigDB suite results (multiple collections)
      cat(sprintf("%s:\n", analysis))
      for (coll_name in names(x[[analysis]])) {
        n_total <- nrow(x[[analysis]][[coll_name]]@result)
        n_sig <- sum(x[[analysis]][[coll_name]]@result$p.adjust < 0.05)
        cat(sprintf("  %s: %d significant / %d total terms\n", coll_name, n_sig, n_total))
      }
    } else {
      # Single analysis results (GO, KEGG, or single MSigDB)
      n_total <- nrow(x[[analysis]]@result)
      n_sig <- sum(x[[analysis]]@result$p.adjust < 0.05)
      cat(sprintf("%s: %d significant / %d total terms\n", analysis, n_sig, n_total))
    }
  }
  
  invisible(x)
}
