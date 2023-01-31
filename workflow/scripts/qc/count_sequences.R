
library(magrittr)
library(tidyverse)
library(ShortRead)

orig <- ShortRead::readFastq(snakemake@input[["original"]])
trim <- ShortRead::readFastq(snakemake@input[["trimmed"]])


tibble::tibble(
  sample = snakemake@params[["sample"]],
  original = length(orig),
  trimmed = length(trim)) %>%
  readr::write_tsv(snakemake@output[["summary"]])
