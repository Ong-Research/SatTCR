#!/usr/bin/Rscript

options(bspm.sudo = TRUE)

"Calculates summary statistics for a single TCR sample from an AIRR file

Usage:
calculate_tcr_stats.R [<outfile>] [--airr=<airr_file>] [--min_count=<mc>]
calculate_tcr_stats.R (-h|--help)
calculate_tcr_stats.R --version

Options:
-h --help    show this screen
--airr=<airr_file>    MIXCR results in a AIRR file format.
--min_count=<mc>    The min. # of repeats per clonotypes [default: 1]" -> doc

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(
  doc, args = my_args, version = "calculate TCR stats V1")

stopifnot(file.exists(arguments$airr))

library(magrittr)
library(tidyverse)

source(here::here("workflow/src/tcr_funs.R"))

arguments[["min_count"]] %<>% as.numeric()


results <- vroom::vroom(arguments[["airr"]])



results %<>% clean_clono()

results_complete <- results %>%
  filter(complete)

results_mc <- results %>%
  filter(count >= arguments[["min_count"]])

results_complete_mc <- results_complete %>%
  filter(count >= arguments[["min_count"]])


list(
  results %>%
    summarize_clonotypes(),
  results_complete %>%
    summarize_clonotypes() %>%
    rename_with(~ str_c(., "complete", sep = "_")),
  results_mc %>%
    summarize_clonotypes() %>%
    rename_with(~ str_c(., "mc", sep = "_")),
  results_complete_mc %>%
    summarize_clonotypes() %>%
    rename_with(~ str_c(., "complete_mc", sep = "_"))) %>%
  bind_cols() %>%
  select(
    starts_with("nclono"),
    starts_with("shannon"),
    starts_with("clonality"),
    starts_with("d50")) %>%
  readr::write_tsv(arguments[["outfile"]])
