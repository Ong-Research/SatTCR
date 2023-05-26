
library(magrittr)
library(tidyverse)
library(GenomicAlignments)
library(ShortRead)


bam_file <- snakemake@input[["sort_bam"]]
r1 <- snakemake@input[["r1"]]
r2 <- snakemake@input[["r2"]]
imgt <- snakemake@input[["imgt"]]
window <- as.numeric(snakemake@params[["window"]])

# bam_file <- "output/align/K9013-02_S112_L002.sorted.bam"
# r1 <- "data/K9013-02_S112_L002_R1_001.fastq"
# r2 <- "data/K9013-02_S112_L002_R2_001.fastq"
# imgt <- list.files("output/align_imgt", pattern = "sorted",
#   full.names = TRUE) %>%
#   stringr::str_subset("bai", negate = TRUE) %>%
#   stringr::str_subset("TRB")
# window <- 100

stopifnot(
  file.exists(bam_file),
  file.exists(r1), file.exists(r2),
  all(file.exists(imgt)))

message("reading alignments")
alignments <- GenomicAlignments::readGAlignments(bam_file, use.names = TRUE)
alignments <- subset(alignments, str_detect(seqnames(alignments), "CM"))

imgt_algn <- imgt %>%
  purrr::map(GenomicAlignments::readGAlignments, use.names = TRUE)
imgt_algn <- do.call(c, imgt_algn) %>%
  GenomicRanges::granges()

message("reading sequence data")
r1_seqs <- ShortRead::readFastq(r1)
r2_seqs <- ShortRead::readFastq(r2)

trb_aligns <- IRanges::subsetByOverlaps(alignments, imgt_algn)
reads_names <- names(trb_aligns)

# get indexes
get_idx <- function(sread, read_names) {

  ss <- as.character(ShortRead::id(sread))
  ss %<>% stringr::str_split(" ")
  ss %<>% unlist()

  ss <- ss[seq(1, length(ss), by = 2)]

  idx <- ss %in% read_names
  idx <- which(idx)

  idx

}

message("getting indexes")
idx1 <- get_idx(r1_seqs, reads_names)
idx2 <- get_idx(r2_seqs, reads_names)

# sanity check
identical(idx1, idx2)

out1 <- r1_seqs[idx1]
out2 <- r2_seqs[idx2]

message("saving results")
ShortRead::writeFastq(out1, snakemake@output[["r1"]])
ShortRead::writeFastq(out2, snakemake@output[["r2"]])
