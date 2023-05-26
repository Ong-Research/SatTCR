library(Biostrings)
library(GenomicAlignments)
library(magrittr)
library(stringr)

fasta_ref <- Biostrings::readDNAStringSet(snakemake@input[["fa"]])
alignment <- GenomicAlignments::readGAlignments(snakemake@input[["align"]],
  use.names = TRUE)

clean_ref <- as.character(fasta_ref) %>%
  str_remove_all("\\.") %>%
  str_remove_all("N")
names(clean_ref) <- names(fasta_ref)
clean_ref <- DNAStringSet(clean_ref)

common_seqs <- intersect(names(fasta_ref), names(alignment))

clean_ref <- clean_ref[common_seqs]
alignment <- alignment[common_seqs]

names(clean_ref) <-
  str_c(names(clean_ref),
    as.character(seqnames(alignment)),
    start(alignment),
    end(alignment),
    as.character(strand(alignment)),
    sep = " ")

Biostrings::writeXStringSet(clean_ref, snakemake@output[["ref"]])
