#' Run GO Gene Set Enrichment Analysis
#'
#' @param ranked_entrez Named numeric vector (Entrez IDs as names)
#' @param ont GO ontology: "BP", "MF", or "CC"
#' @param simplify Use semantic similarity to reduce redundancy?
#' @param simplify_cutoff Similarity cutoff for simplification
#' @param simplify_measure Similarity measure: "Wang", "Rel", etc.
#' @param minGSSize Minimum gene set size
#' @param maxGSSize Maximum gene set size
#' @param pvalueCutoff P-value cutoff
#' @param eps Adaptive permutation threshold
#' @param nPermSimple Number of permutations
#' @param verbose Print progress?
#'
#' @return gseaResult object from clusterProfiler
#' @export
run_gsea_go <- function(
  ranked_entrez,
  ont = "BP",
  simplify = FALSE,
  simplify_cutoff = 0.7,
  simplify_measure = "Wang",
  minGSSize = 15,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  eps = 1e-10,
  nPermSimple = 10000,
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

  # Check for duplicates (shouldn't happen if prep was done properly)
  if (any(duplicated(names(ranked_entrez)))) {
    warning("Duplicate ENTREZ IDs detected. Keeping first occurrence only.")
    ranked_entrez <- ranked_entrez[!duplicated(names(ranked_entrez))]
  }

  # Run GO GSEA
  if (verbose) message("Running GO GSEA (", ont, ")...")

  gsea_result <- clusterProfiler::gseGO(
    geneList     = ranked_entrez,
    OrgDb        = org.Hs.eg.db::org.Hs.eg.db,
    ont          = ont,
    minGSSize    = minGSSize,
    maxGSSize    = maxGSSize,
    eps          = eps,
    nPermSimple  = nPermSimple,
    pvalueCutoff = pvalueCutoff,
    verbose      = verbose
  )

  # Optional simplification
  if (simplify && nrow(gsea_result) > 0) {
    if (verbose) message("Simplifying GO results with ", simplify_measure, " similarity...")

    # Load semantic similarity data
    hsGO <- GOSemSim::godata(annoDb = org.Hs.eg.db, ont = ont)

    gsea_result <- clusterProfiler::simplify(
      gsea_result,
      cutoff     = simplify_cutoff,
      by         = "p.adjust",
      select_fun = min,
      measure    = simplify_measure,
      semData    = hsGO
    )
  }

  if (verbose) {
    n_sig <- sum(gsea_result@result$p.adjust < 0.05)
    message("Completed: ", n_sig, " significant terms (padj < 0.05)")
  }

  return(gsea_result)
}