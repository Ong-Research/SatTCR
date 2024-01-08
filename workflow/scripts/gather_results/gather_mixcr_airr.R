

library(magrittr)
library(tidyverse)
library(vroom)

library(future)
library(furrr)

plan("multicore", workers = snakemake@threads)

airr_file <- snakemake@input[["airr"]]
summary_file <- snakemake@input[["summary"]]

all_files <- list.files(snakemake@params[["bootdir"]], full.names = TRUE,
  pattern = "airr.tsv", recursive = TRUE)
sample <- snakemake@params[["sample"]]

sample_files <- str_subset(all_files, sample)

airr <- vroom::vroom(airr_file)
boot_summary <- summary_file %>%
  qs::qread()

boot_summary %<>%
  mutate(
    prefix = basename(r1),
    prefix = str_split(prefix, "_R1"),
    prefix = str_c("/", map_chr(prefix, 1), "/"),
    prefix = map_chr(prefix, ~ str_subset(sample_files, .)),
    clono = map(prefix, vroom::vroom), .progress = TRUE) %>%
  select(prefix, everything()) %>%
  select(-r1, -r2, -seq_blocks, -.progress)

# general statistics

clean_clono <- function(clono) {

  clono %<>%
    select(-ends_with("score"), -contains("alignment"),
      -contains("germline"), -ends_with("cigar"),
      -starts_with("fwr"), -ends_with("end"), -ends_with("start"),
      -starts_with("cdr1"), -starts_with("cdr2"), -starts_with("np1"),
      -starts_with("np2"), -rev_comp, -complete_vdj)

  clono %>%
    rename(complete = productive, count = duplicate_count)

}



## nclonotypes
count_clonotypes <- function(clono, filter_complete = TRUE) {

  if (filter_complete) {
    clono %<>%
      filter(complete)
  }

  clono %>%
    distinct(cdr3, v_call, j_call) %>%
    nrow()

}

## D50
calculate_d50 <- function(clono, filter_complete = TRUE) {

  if (filter_complete) {
    clono %<>%
      filter(complete)
  }

  clono %>%
    group_by(cdr3, v_call, j_call) %>%
    summarize(
      count = sum(count), .groups = "drop") %>%
    arrange(desc(count)) %>%
    mutate(
      id = row_number(),
      acc_count = cumsum(count) / sum(count)) %>%
    filter(acc_count >= 0.5) %>%
    pull(id) %>%
    min()

}

## Shannon diversity
calculate_shannon <- function(clono, filter_complete = TRUE) {

  if (filter_complete) {
    clono %<>%
      filter(complete)
  }

  clono %<>%
    group_by(cdr3, v_call, j_call) %>%
    summarize(
      count = sum(count), .groups = "drop") %>%
    arrange(desc(count)) %>%
    mutate(
      p = count / sum(count))

  p <- clono[["p"]]
  - sum(p * log(p))


}

boot_summary %<>%
  mutate(
    clono = furrr::future_map(clono, clean_clono)) %>%
  crossing(filter_complete = c(TRUE, FALSE), min_count = c(1, 5, 10)) %>%
  mutate(
    clono = furrr::future_map2(clono, min_count, ~ filter(.x, count >= .y),
      .progress = TRUE),
    nclonotypes = furrr::future_map2_dbl(clono, filter_complete,
      count_clonotypes, .progress = TRUE),
    d50 = furrr::future_map2_dbl(clono, filter_complete, calculate_d50),
    shannon = furrr::future_map2_dbl(clono, filter_complete,
      calculate_shannon, .progress = TRUE),
    exp_shannon = exp(shannon))

# comparison with full clonotypes

full_clono <- airr_file %>%
  vroom::vroom() %>%
  clean_clono()

filter_clono <- function(clono, filter_complete, min_count) {

  if (filter_complete) {
    clono %<>%
      filter(complete)
  }

  clono %>%
    filter(count >= min_count)

}

## jaccard comparison
boot_summary %<>%
  mutate(
    full_clono = map2(filter_complete, min_count,
      ~ filter_clono(full_clono, .x, .y)))

## jaccard
compute_jaccard <- function(clono, full_clono) {

  clono %<>%
    select(v_call, j_call, cdr3)
  full_clono %<>%
    select(v_call, j_call, cdr3)


  all_clono <- clono %>%
    full_join(full_clono, by = c("v_call", "j_call", "cdr3"),
      relationship = "many-to-many") %>%
    distinct(v_call, j_call, cdr3) %>%
    nrow()

  both_clono <- clono %>%
    inner_join(full_clono, by = c("v_call", "j_call", "cdr3"),
      relationship = "many-to-many") %>%
    distinct(v_call, j_call, cdr3) %>%
    nrow()

  both_clono / all_clono

}

boot_summary %<>%
  mutate(
    jaccard = furrr::future_map2_dbl(clono, full_clono, compute_jaccard,
      .progress = TRUE))

compute_spearman_corr <- function(clono, full_clono) {

  clono %<>%
    select(v_call, j_call, cdr3, count)

  full_clono %<>%
    select(v_call, j_call, cdr3, count)

  both_clono <- clono %>%
    inner_join(full_clono, by = c("v_call", "j_call", "cdr3"),
      relationship = "many-to-many")

  cc <- cor(both_clono["count.x"], both_clono[["count.y"]],
    method = "spearman")

  rlang::set_names(cc[1,1], NULL)

}

boot_summary %<>%
  mutate(
    spearman_cor = furrr::future_map2_dbl(clono, full_clono,
      compute_spearman_corr, .progress = TRUE)
  )


boot_summary %>%
  select(-full_clono) %>%
  qs::qsave(snakemake@output[["out"]])
