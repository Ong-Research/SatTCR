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

r1 <- snakemake@input[["r1"]]
r2 <- snakemake@input[["r2"]]

stopifnot(file.exists(r1), file.exists(r2))

sample_perc <- as.numeric(snakemake@params[["sample_perc"]])
seed <- snakemake@params[["seed"]]

message("reading files")
r1_seqs <- ShortRead::readFastq(r1)
r2_seqs <- ShortRead::readFastq(r2)

n <- length(r1_seqs)
stopifnot(identical(n, length(r2_seqs)))

set.seed(seed)
m <- floor(sample_perc * n)
message("sampling ", m, " pairs of reads")
idx <- sample(n, m, replace = FALSE)

r1_out <- r1_seqs[idx]
r2_out <- r2_seqs[idx]

message("saving files")
ShortRead::writeFastq(r1_out, snakemake@output[["r1"]])
ShortRead::writeFastq(r2_out, snakemake@output[["r2"]])
