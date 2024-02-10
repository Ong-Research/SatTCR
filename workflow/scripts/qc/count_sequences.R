#!/usr/bin/Rscript

 options(bspm.sudo = TRUE)

"Count # of sequences in raw and trimmed sequence files

Usage:
count_sequences.R [<outfile>] [--input=<input> --trimmed=<trimmed>] [--sample=<sample>]
count_sequences.R (-h|--help)
count_sequences.R --version

Options:
-h --help    show this screen
--input=<input>    Name of the raw input file
--trimmed=<trimmed>    Name of the trimmed file
--sample=<sample>    Name of the sample" -> doc

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(
  doc, args = my_args, version = "count sequences V1")

stopifnot(
  file.exists(arguments$input),
  is.null(arguments$trimmed) | file.exists(arguments$trimmed))

library(magrittr)
library(tidyverse)
library(ShortRead)

orig <- ShortRead::readFastq(arguments$input)

trim_flag <- is.null(arguments$trimmed)
if (! trim_flag) {
  trim <- ShortRead::readFastq(arguments$trimmed)
}

tibble::tibble(
  sample = arguments$sample,
  original = length(orig),
  trimmed = ifelse(!trim_flag, length(trim), 0)) %>%
  readr::write_tsv(arguments$outfile)
