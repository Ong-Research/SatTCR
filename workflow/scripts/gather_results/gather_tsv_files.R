
library(magrittr)
library(tidyverse)

stopifnot(snakemake@params[["method"]] %in% c("trust4", "mixcr"))

results <- tibble::tibble(file = unlist(snakemake@input, vroom::vroom))
results %<>%
  dplyr::mutate(report = purrr::map(file, vroom::vroom)) %>%
  tidyr::unnest(cols = c(report))

if (stringr::str_to_lower(snakemake@params[["method"]]) == "trust4") {

  results %<>%
    dplyr::rename_with(snakecase::to_snake_case) %>%
    dplyr::rename(cdr3_nt = cdr_3_nt, cdr3_aa = cdr_3_aa, v_gene = v,
      d_gene = d, j_gene = j, c_gene = c) %>%
    dplyr::select(file, count, starts_with("cdr3"), ends_with("gene")) %>%
    tidyr::nest(report = c(count, cdr3_nt, cdr3_aa, v_gene, j_gene,
      d_gene, c_gene))

} else {

  remove_score <- function(gene_col) {

    split_genes <- gene_col %>%
      stringr::str_split("\\,")
    split_genes <- purrr::map(split_genes, stringr::str_split, "\\(")
    split_gene_names <- purrr::map(split_genes,
      ~ purrr::map(., 1))
    purrr::map_chr(split_gene_names,
      ~ stringr::str_c(., collapse = ";"))

  }

    # dplyr::mutate(
    #   across(ends_with("gene"),
    #     list(~ remove_score(.)), .names = "{.col}")) %>%


  results %<>%
    dplyr::rename_with(snakecase::to_snake_case) %>%
    dplyr::rename(
      count = read_count,
      cdr3_nt = n_seq_cdr_3,
      cdr3_aa = aa_seq_cdr_3,
      v_gene = all_v_hits_with_score,
      d_gene = all_d_hits_with_score,
      j_gene = all_j_hits_with_score,
      c_gene = all_c_hits_with_score) %>%
    dplyr::select(file, count, starts_with("cdr3"), ends_with("gene")) %>%
    tidyr::nest(report = c(count, cdr3_nt, cdr3_aa, v_gene, j_gene,
      d_gene, c_gene))



}

results %>%
  qs::qsave(snakemake@output[[1]])
