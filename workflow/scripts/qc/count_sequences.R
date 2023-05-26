
library(magrittr)
library(tidyverse)
library(ShortRead)

orig <- ShortRead::readFastq(snakemake@input[["original"]])

trim_flag <- is.null(snakemake@input[["trimmed"]])
if (! is.null(snakemake@input[["trimmed"]])) {
  trim <- ShortRead::readFastq(snakemake@input[["trimmed"]])
}

tibble::tibble(
  sample = snakemake@params[["sample"]],
  original = length(orig),
  trimmed = ifelse(!trim_flag, length(trim), 0)) %>%
  readr::write_tsv(snakemake@output[["summary"]])
