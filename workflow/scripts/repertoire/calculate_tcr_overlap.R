#!/usr/bin/Rscript

options(bspm.sudo = TRUE)

"Calculates summary statistics for a single TCR sample from an AIRR file

Usage:
calculate_tcr_overlap.R [<outfile>] [--airr1=<airr_file1> --airr2=<airr_file2>] [--min_count=<mc>]
calculate_tcr_overlap.R (-h|--help)
calculate_tcr_overlap.R --version

Options:
-h --help    show this screen
--airr1=<airr_file1>    A set of MIXCR results in a AIRR file format.
--airr2=<airr_file2>    A different set of MIXCR results in a AIRR file format.
--min_count=<mc>    The min. # of repeats per clonotypes [default: 1]" -> doc

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(
  doc, args = my_args, version = "calculate TCR stats V1")

stopifnot(file.exists(arguments$airr1), file.exists(arguments$airr2))

library(magrittr)
library(tidyverse)

source(here::here("workflow/src/tcr_funs.R"))

arguments[["min_count"]] %<>% as.numeric()

tcr1 <- vroom::vroom(arguments[["airr1"]])
tcr2 <- vroom::vroom(arguments[["airr2"]])

# remove unused columns
tcr1 %<>% clean_clono()
tcr2 %<>% clean_clono()


compute_overlap <- function(tcr1, tcr2, mc = 1) {

  # define a unique clonotype by having shared
  # V gene call, J gene call and CDR3 sequence
  keys <- c("v_call", "j_call", "cdr3")
  tcr1 %<>%
    filter(count >= mc) %>%
    distinct(!!! rlang::syms(keys))
  tcr2 %<>%
    filter(count >= mc) %>%
    distinct(!!! rlang::syms(keys))

tibble::tribble(
  ~ total1, ~ total2, ~ only1, ~ only2, ~ both,
  ~ universe,
  nrow(tcr1), nrow(tcr2),
  nrow(setdiff(tcr1, tcr2)),
  nrow(setdiff(tcr2, tcr1)),
  nrow(intersect(tcr1, tcr2)),
  nrow(union(tcr1, tcr2))) %>%
  mutate(min_count = mc)

}

tcr1_complete <- tcr1 %>%
  filter(complete)
tcr2_complete <- tcr2 %>%
  filter(complete)


# Use set-operations to check overlaps
bind_rows(
  compute_overlap(tcr1, tcr2, 1) %>%
    mutate(clones = "all"),
  compute_overlap(
    tcr1_complete, tcr2_complete, 1) %>%
    mutate(clones = "complete"),
  compute_overlap(tcr1, tcr2, 5) %>%
    mutate(clones = "all"),
  compute_overlap(
    tcr1_complete, tcr2_complete, 5) %>%
    mutate(clones = "complete")) %>%
  mutate(
    tcr1 = arguments$airr1,
    tcr2 = arguments$airr2) %>%
  select(starts_with("tcr"), everything()) %>%
  readr::write_tsv(arguments[["outfile"]])
