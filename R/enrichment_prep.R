#' Prepare Gene Rankings for GSEA
#'
#' Maps gene symbols to Entrez IDs and creates a ranked dataframe.
#'
#' @param deseq_dds DESeq2 results data frame
#' @param gene_col Column name containing gene symbols
#' @param stat_col Column name containing ranking statistic
#' @param add_noise Add tiny noise to break ties?
#' @param noise_sd Standard deviation of noise
#' @param keep_duplicates How to handle duplicate mappings
#' @param verbose Print progress messages?
#'
#' @return Tibble with SYMBOL, stat, ENTREZID columns
#' @export
#' @importFrom dplyr select filter inner_join group_by ungroup arrange mutate
#' @importFrom tibble as_tibble
prep_ranks <- function(
  deseq_dds,
  gene_col = "Gene",
  stat_col = "stat",
  add_noise = TRUE,
  noise_sd = 1e-12,
  keep_duplicates = "max_abs",
  verbose = TRUE
) {
  
  # Input validation
  if (!gene_col %in% colnames(deseq_dds)) {
    stop("Column '", gene_col, "' not found in input data")
  }
  if (!stat_col %in% colnames(deseq_dds)) {
    stop("Column '", stat_col, "' not found in input data")
  }
  
  # Map SYMBOL to ENTREZ in one go
  if (verbose) message("Mapping gene symbols to Entrez IDs...")
  
  gene_mappings <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = deseq_dds[[gene_col]],
    keytype = "SYMBOL",
    columns = c("ENTREZID", "SYMBOL")
  ) |>
    as_tibble() |>
    # Remove NA mappings immediately
    filter(!is.na(ENTREZID), !is.na(SYMBOL))
  
  # Join back to original data
  rank_df <- deseq_dds |>
    select(SYMBOL = all_of(gene_col), stat = all_of(stat_col)) |>
    filter(is.finite(stat)) |>
    inner_join(gene_mappings, by = "SYMBOL", relationship = "many-to-many") |>
    filter(!is.na(ENTREZID))
  
  if (verbose) {
    n_original <- nrow(deseq_dds)
    n_mapped <- n_distinct(rank_df$SYMBOL)
    n_lost <- n_original - n_mapped
    message(sprintf("  %d/%d genes mapped successfully (%.1f%% loss)",
                    n_mapped, n_original, 100 * n_lost / n_original))
  }
  
  # Handle duplicate mappings (one symbol -> multiple Entrez or vice versa)
  if (keep_duplicates == "max_abs") {
    rank_df <- rank_df |>
      group_by(ENTREZID) |>
      slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) |>
      ungroup() |>
      # Also deduplicate by SYMBOL in case multiple Entrez map to same symbol
      group_by(SYMBOL) |>
      slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) |>
      ungroup()
  } else if (keep_duplicates == "first") {
    rank_df <- rank_df |>
      distinct(ENTREZID, .keep_all = TRUE) |>
      distinct(SYMBOL, .keep_all = TRUE)
  } else {
    stop("keep_duplicates must be 'max_abs' or 'first'")
  }
  
  # Add tiny noise to break ties if requested
  if (add_noise) {
    rank_df <- rank_df |>
      mutate(stat = stat + rnorm(n(), mean = 0, sd = noise_sd))
  }
  
  # Sort by stat (descending)
  rank_df <- rank_df |>
    arrange(desc(stat))
  
  if (verbose) {
    message("Final dataset: ", nrow(rank_df), 
            " unique genes with both SYMBOL and ENTREZID")
  }
  
  # Return tibble with both ID types
  return(rank_df)
}


#' Extract Named Ranking Vector
#'
#' @param rank_df Data frame from prep_ranks()
#' @param id_type Either "ENTREZID" or "SYMBOL"
#' @param value Column to use as values (default "stat")
#'
#' @return Named numeric vector
#' @export
#' @importFrom tibble deframe
deframe_ranks <- function(rank_df, id_type = c("ENTREZID", "SYMBOL"), 
                          value = "stat") {
  
  id_type <- match.arg(id_type)
  
  rank_df |>
    select({{id_type}}, {{value}}) |>
    deframe()
}
