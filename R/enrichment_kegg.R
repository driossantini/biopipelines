#' Run KEGG Gene Set Enrichment Analysis
#'
#' @param ranked_entrez Named numeric vector (Entrez IDs as names)
#' @param organism KEGG organism code (default "hsa")
#' @param use_internal_data Use internal KEGG data instead of API?
#' @param ... Additional parameters passed to gseKEGG
#'
#' @return gseaResult object
#' @export
run_gsea_kegg <- function(
  ranked_entrez,
  organism = "hsa",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 1e-10,
  nPermSimple = 10000,
  use_internal_data = FALSE,
  verbose = FALSE
) {

  # Input validation
  if (!is.numeric(ranked_entrez)) {
    stop("ranked_entrez must be a numeric vector")
  }
  if (is.null(names(ranked_entrez))) {
    stop("ranked_entrez must be a named vector (ENTREZ IDs as names)")
  }
  if (any(is.na(ranked_entrez)) || any(!is.finite(ranked_entrez))) {
    warning("Removing non-finite values from ranked_entrez")
    ranked_entrez <- ranked_entrez[is.finite(ranked_entrez)]
  }

  # Check for duplicates
  if (any(duplicated(names(ranked_entrez)))) {
    warning("Duplicate ENTREZ IDs detected. Keeping first occurrence only.")
    ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]
  }

  # Run KEGG GSEA
  if (verbose) message("Running KEGG GSEA (organism: ", organism, ")...")

  gsea_result <- clusterProfiler::gseKEGG(
    geneList     = ranked_entrez,
    organism     = organism,
    minGSSize    = minGSSize,
    maxGSSize    = maxGSSize,
    eps          = eps,
    nPermSimple  = nPermSimple,
    pvalueCutoff = pvalueCutoff,
    use_internal_data = use_internal_data,
    verbose      = verbose
  )

  if (verbose) {
    n_sig <- sum(gsea_result@result$p.adjust < 0.05)
    message("Completed: ", n_sig, " significant pathways (padj < 0.05)")
  }

  return(gsea_result)
}
