#!/usr/bin/Rscript

options(bspm.sudo = TRUE)

"Sample sequences from a pair of fastq files under a sequential bootstrap approach

Usage:
sample_sequences.R [<outdir>] [--summary=<summary>] [--r1=<r1_file> --r2=<r2_file>] [--block_size=<block_size> --nblocks=<nblocks>] [--nreplicates=<nreps> --seed=<seed> --sample=<sample>]
sample_sequences.R (-h|--help)
sample_sequences.R --version

Options:
-h --help    show this screen
--summary=<summary>    location of the summary file with the # of sequence per pair of files
--r1=<r1_file>    name of the R1 end file
--r2=<r2_file>    name of the R2 end file
--block_size=<block_size>   # of sequences sampled per block or
--nblocks=<nblocks>    # of blocks to divide the R1 and R2 files
--nreplicates=<nreps>    # of bootstrap replicates
--seed=<seed>    seed to initialize random numbers
--sample=<sample>   sample prefix used to name the files" -> doc

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(
  doc, args = my_args, version = "sample sequences V1")

if (arguments$nblocks == "None") arguments$nblocks <- NULL
if (arguments$block_size == "None") arguments$block_size <- NULL

stopifnot(
  file.exists(arguments$r1),
  file.exists(arguments$r2),
  ! is.null(arguments$seed),
  ! is.null(arguments$nreplicates),
  ! is.null(arguments$sample),
  ! is.null(arguments$block_size) | ! is.null(arguments$nblocks))

if (! is.null(arguments$block_size) & ! is.null(arguments$nblocks)) {
  stop("both the block size and the # of blocks are non-null, plese pick one")
}


library(magrittr)
library(tidyverse)
library(ShortRead)
library(yaml)

r1 <- arguments[["r1"]]
r2 <- arguments[["r2"]]
seed <- as.numeric(arguments$seed)
set.seed(as.numeric(seed))

message("reading files")
r1_seqs <- ShortRead::readFastq(r1)
r2_seqs <- ShortRead::readFastq(r2)

n <- length(r1_seqs)
message("sequencing depth: ", scales::comma(n))

boots_replicates <- as.numeric(arguments$nreplicates)

stopifnot(identical(n, length(r2_seqs)))

if (! is.null(arguments$nblocks)) {

  nblocks <- as.numeric(arguments$nblocks)
} else {

  block_size <- as.numeric(arguments$block_size)
  nblocks <- floor(n / block_size) + 1
}

partition <- sample(nblocks, size = n, replace = TRUE)

# get blocks
idx_split <- split(seq_len(n), partition)
r1_splits <- map(idx_split, ~ r1_seqs[.])
r2_splits <- map(idx_split, ~ r2_seqs[.])

seq_blocks <- replicate(boots_replicates, {
  random_blocks <- sample(nblocks)
  purrr::map(seq_along(random_blocks), ~ random_blocks[seq_len(.)])
}, simplify = FALSE)

seq_blocks %<>% unlist(recursive = FALSE)

merge_shortreadq_list <- function(shortq_list) {

  sr_list <- purrr::map(shortq_list, ShortRead::sread)
  qu_list <- purrr::map(shortq_list, Biostrings::quality)
  id_list <- purrr::map(shortq_list, ShortRead::id)

  ShortRead::ShortReadQ(
    sread = Reduce(c, sr_list),
    quality = Reduce(append, qu_list),
    id = Reduce(c, id_list))
}

if (is.null(arguments[["outdir"]])) {
  outdir <- "."
} else {
  outdir <- arguments[["outdir"]]
}

fs::dir_create(outdir)

sample <- arguments[["sample"]]

files <- glue::glue(
  "{outdir}/{sample}_{id}_",
  outdir = outdir, sample = sample,
  id = seq_along(seq_blocks)) %>%
  as.character()

files_r1 <- str_c(files, "R1.fastq.gz")
files_r2 <- str_c(files, "R2.fastq.gz")

partition_counts <- table(partition)

summary_tibble <- tibble(
  r1 = files_r1, r2 = files_r2, seq_blocks) %>%
  mutate(
    boot_id = rep(seq_len(boots_replicates), each = nblocks),
    block_id = rep(seq_len(nblocks), boots_replicates),
    expected_seq_depth = map_dbl(seq_blocks,
      ~ sum(partition_counts[as.character(.)])))

message("Generating ", nrow(summary_tibble), " pairs of files")

write_and_save <- function(block, file_r1, file_r2,
  r1_splits, r2_splits) {

  r1_out <- r1_splits[as.character(block)] %>%
    merge_shortreadq_list()
  r2_out <- r2_splits[as.character(block)] %>%
    merge_shortreadq_list()

  ShortRead::writeFastq(r1_out, file_r1)
  ShortRead::writeFastq(r2_out, file_r2)

  # double check
  read_out <- ShortRead::readFastq(file_r1)
  length(read_out)

}

summary_tibble %<>%
  mutate(
    file_depth = pmap_dbl(
      list(seq_blocks, r1, r2),
      write_and_save, r1_splits, r2_splits, .progress = TRUE))

summary_tibble %>%
  qs::qsave(arguments[["summary"]])
