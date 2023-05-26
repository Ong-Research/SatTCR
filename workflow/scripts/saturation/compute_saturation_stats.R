
library(magrittr)
library(tidyverse)
library(Biostrings)
library(repertoire)
library(future)
library(furrr)

plan(multicore, workers = snakemake@threads)

conflicted::conflict_prefer_all("dplyr")
conflicted::conflicts_prefer(Biostrings::pattern)
conflicted::conflicts_prefer(purrr::set_names)


nreads <- snakemake@input[["nseqs"]] %>%
  vroom::vroom()

# clean simple report
report <- snakemake@input[["simple_report"]] %>%
  vroom::vroom() %>%
  rename_with(snakecase::to_snake_case) %>%
  rename(cdr3aa = cdr_3_aa, cdr3nt = cdr_3_nt) %>%
  is_complete(cdr3nt) %>%
  filter(str_detect(v, regex("^TRBV"))) %>%
  filter(str_detect(j, regex("^TRBJ"))) %>%
  mutate(
    libsize = sum(count),
    seq_depth = nreads[["original"]])

unique_seqs <- report %>%
  filter(is_complete == "yes") %>%
  filter(v != "." & j != ".") %>%
  distinct(cdr3nt, v, j, c) %>%
  nest(seqs = c(cdr3nt))

imgt_file <- glue::glue("output/trust4/imgt_annots/{species}/imgt.fa",
  species = snakemake@config[["species"]]) %>%
  as.character()

imgt <- readDNAStringSet(imgt_file)
imgt <- imgt[str_detect(names(imgt), regex("^TRB"))]

extract_table <- function(seqs, algn) {

  seqs %>%
    mutate(
      v_seq = as.character(
        as(pattern(algn[["v_align"]]), "DNAStringSet")),
      v_annot = as.character(
        as(subject(algn[["v_align"]]), "DNAStringSet")),
      j_seq = as.character(
        as(pattern(algn[["j_align"]]), "DNAStringSet")),
      j_annot = as.character(
        as(subject(algn[["j_align"]]), "DNAStringSet")),
      v_width = nchar(v_seq),
      j_width = nchar(j_seq),
      v_pid = pid(algn[["v_align"]]),
      j_pid = pid(algn[["j_align"]]),
      v_miss = nmismatch(algn[["v_align"]]),
      j_miss = nmismatch(algn[["j_align"]]),
      v_score = score(algn[["v_align"]]),
      j_score = score(algn[["j_align"]]))

}

unique_seqs %<>%
  mutate(
    algn = future_pmap(list(v, j, seqs), align_seqs, imgt),
    seqs = future_map2(seqs, algn, extract_table))

report %<>%
  left_join(
    unique_seqs %>%
      select(-algn) %>%
      unnest(cols = c(seqs)), by = c("v", "j", "c", "cdr3nt")) %>%
  mutate(width_nt = nchar(cdr3nt), width_aa = nchar(cdr3aa))


compare_seqs <- function(seq, annot, bp) {

  DNAStringSet(str_sub(seq, 1, round(bp, 0))) ==
    DNAStringSet(str_sub(annot, 1, round(bp, 0)))

}

ov_perc <- snakemake@config[["saturation"]][["perc_overlap_width"]]
med_perc <- snakemake@config[["saturation"]][["perc_med_width"]]

report %<>%
  filter(is_complete == "yes") %>%
  select(-cid_full_length) %>%
  mutate(
    v_match = compare_seqs(v_seq, v_annot, ov_perc * median(v_width)),
    j_match = compare_seqs(j_seq, j_annot, ov_perc * median(j_width))) %>%
  filter(v_width > med_perc * quantile(v_width, .5) &
    j_width > med_perc * quantile(j_width, .5)) %>%
  filter(v_match) %>%
  filter(j_match)

compute_libsize <- function(clono_str, min_count, report) {

  report %>%
    group_by(!! sym(clono_str)) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    filter(count >= min_count) %>%
    pull(count) %>%
    sum(na.rm = TRUE)

}


count_clonotype <- function(clono_str, min_count, report) {

  report %>%
    group_by(!! sym(clono_str)) %>%
    summarise(count = sum(count), .groups = "drop") %>%
    filter(count >= min_count) %>%
    nrow()

}

compute_shannon <- function(clono_str, min_count, report) {

  report %>%
    group_by(!! sym(clono_str)) %>%
    summarize(count = sum(count), .groups = "drop") %>%
    mutate(freq = count / sum(count)) %>%
    summarize(
      shannon = - (1 / log(n())) * sum(freq * log(freq)),
        .groups = "drop") %>%
    pull(shannon)

}

crossing(min_count = c(1, 5, 10), cdr3 = c("cdr3nt", "cdr3aa")) %>%
  mutate(seq_depth = nreads[["original"]]) %>%
  mutate(
    libsize = map2_dbl(cdr3, min_count, compute_libsize, report),
    nclonotype = map2_dbl(cdr3, min_count, count_clonotype, report),
    cpk = nclonotype / (seq_depth / 1000),
    shannon = map2_dbl(cdr3, min_count, compute_shannon, report)) %>%
  qs::qsave(snakemake@output[["stats"]])
