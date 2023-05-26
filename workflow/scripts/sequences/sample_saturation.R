#-------------------------------------------------------------------------------
#
#   Script to sample without replacement a proportion of the totall of the
#   # of sequences in a pair of FASTQ files
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

library(magrittr)
library(tidyverse)
library(ShortRead)
library(yaml)

r1 <- snakemake@input[["r1"]]
r2 <- snakemake@input[["r2"]]

stopifnot(file.exists(r1), file.exists(r2))

seed <- snakemake@params[["seed"]]

message("reading files")
r1_seqs <- ShortRead::readFastq(r1)
r2_seqs <- ShortRead::readFastq(r2)

n <- length(r1_seqs)
boots_replicates <- snakemake@config[["saturation"]][["bootstrap_replicates"]]

stopifnot(identical(n, length(r2_seqs)))
block_size <- snakemake@config[["saturation"]][["seq_depth_blocks"]]

nblocks <- floor(n / block_size)

set.seed(as.numeric(seed))
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



fs::dir_create(outdir <<- snakemake@output[["outdir"]])
files <- glue::glue(
  "{outdir}/{sample}_{id}_",
  outdir = outdir, sample = snakemake@params[["sample"]],
  id = seq_along(seq_blocks)) %>%
  as.character()

files_r1 <- str_c(files, "R1.fastq.gz")
files_r2 <- str_c(files, "R2.fastq.gz")

write_and_save <- function(block, file_r1, file_r2,
  r1_splits, r2_splits) {

  r1_out <- r1_splits[block] %>%
    merge_shortreadq_list()
  r2_out <- r2_splits[block] %>%
    merge_shortreadq_list()

  ShortRead::writeFastq(r1_out, file_r1)
  ShortRead::writeFastq(r2_out, file_r2)

}

pmap(list(seq_blocks, files_r1, files_r2),
  write_and_save, r1_splits, r2_splits)
