# save(snakemake, file = "zz.RData")

# load("zz.RData")

library(magrittr)
library(tidyverse)
library(vroom)
library(Biostrings)
library(GenomicAlignments)
library(fuzzyjoin)


conflicted::conflict_prefer_all("dplyr")
conflicted::conflicts_prefer(magrittr::set_names)
# conflicted::conflicts_prefer(Biostrings::pattern)
# conflicted::conflicts_prefer(purrr::set_names)

stopifnot(
  file.exists(snakemake@input[["imgt"]]),
  file.exists(snakemake@input[["align"]]),
  file.exists(snakemake@input[["gtf"]]))

imgt <- readDNAStringSet(snakemake@input[["imgt"]])
align <- readGAlignments(snakemake@input[["align"]], use.names = TRUE)

gtf <- vroom(snakemake@input[["gtf"]], skip = 5, delim = "\t",
  col_names = c("seqname", "source", "feature", "start", "end", "score",
    "strand", "frame", "attribute"))
gtf %<>%
  filter(feature == "gene")

get_gene_name <- function(attribute) {

  splits <- attribute %>%
    str_split("\\;")
  gene_names <- map(splits, str_subset, "gene_name")
  map_chr(gene_names,
    ~ ifelse(length(.) == 0, NA_character_, .)) %>%
    str_trim() %>%
    str_remove_all("\"") %>%
    str_remove("gene_name ")

}

gtf %<>%
  mutate(
    gene_name = get_gene_name(attribute)) %>%
  filter(!is.na(gene_name))

imgt_names <- names(imgt) %>%
  str_split(regex("[\\*|\\-]")) %>%
  map_chr(1) %>%
  unique()

gtf_repertoire <- map_df(imgt_names,
  ~ filter(gtf, str_detect(gene_name, .))) %>%
  distinct() %>%
  filter(str_detect(gene_name, regex("[V|D|J|C]")))

align %<>% granges()


align_genes <- tibble(gene_name = names(align)) %>%
  stringdist_left_join(
    gtf_repertoire, by = "gene_name", max_dist = 3) %>%
  rename(
    gene.algn = gene_name.x,
    gene.gtf = gene_name.y) %>%
  filter(str_sub(gene.algn, 1, 6) == str_sub(gene.gtf, 1, 6)) 

gtf_genes <- GRanges(
  seqnames = align_genes$seqname,
  ranges = IRanges(
    start = align_genes$start,
    end = align_genes$end), strand = align_genes$strand)
names(gtf_genes) <- align_genes$gene.algn

align <- align[!names(align) %in% names(gtf_genes)]
align <- c(align, gtf_genes)

imgt_ref <- imgt[names(imgt) %in% names(align)]
imgt_seqs <- imgt_ref %>%
  as.character() %>%
  set_names(NULL) %>%
  str_remove_all("\\.") %>%
  str_remove_all("N")

imgt_names <- tibble(gene = names(imgt_ref)) %>%
  inner_join(
    align %>%
      as.data.frame() %>%
      as_tibble(rownames = "gene"), by = "gene") %>%
  mutate(
    name = str_c(gene, seqnames, start, end, strand, sep = " ")
  )

imgt_ref <- DNAStringSet(imgt_seqs)
names(imgt_ref) <- imgt_names$name

writeXStringSet(imgt_ref, snakemake@output[["trust4_fa"]])
imgt_names %>%
  select(gene) %>%
  write_tsv(snakemake@output[["name_file"]], col_names = FALSE)



# imgt ref
# imgt aligned
# gtf file

# this creates a fasta file with entries of the form
# >IGKV4-1 chr2 88885397 88886153 +
# sequence


