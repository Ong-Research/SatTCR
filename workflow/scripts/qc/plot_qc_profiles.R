

library(magrittr)
library(tidyverse)
library(dada2)
library(cowplot)

r1_file <- snakemake@input[["r1"]]
r2_file <- snakemake@input[["r2"]]

qc_plot <- cowplot::plot_grid(
  dada2::plotQualityProfile(r1_file),
  dada2::plotQualityProfile(r2_file), nrow = 1)

ggsave(
  filename = snakemake@output[["plot"]],
  plot = qc_plot,
  width = snakemake@params[["width"]],
  height = snakemake@params[["height"]],
  units = snakemake@params[["units"]])
