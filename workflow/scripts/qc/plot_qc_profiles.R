#!/usr/bin/Rscript

 options(bspm.sudo = TRUE)

"Plots quality profiles using dada2

Usage:
plot_qc_profiles.R [<outfile>] [--end1=<r1_file> --end2=<r2_file>] [--width=<width> --height=<height>]
plot_qc_profiles.R (-h|--help)
plot_qc_profiles.R --version

Options:
-h --help    show this screen
--end1=<end1>    name of the R1 end fastq.gz file
--end2=<end2>    name of the R2 end fastq.gz file
--width=<width>    width of the plot in inches
--height=<height>    height of the plot in inches" -> doc

my_args <- commandArgs(trailingOnly = TRUE)

arguments <- docopt::docopt(
  doc, args = my_args, version = "plot quality profiles V1")

if (interactive()) {

  arguments$end1 <- "data/run1/K9013rcc-rxn3_S423_L003_R1_001.fastq.gz"
  arguments$end2 <- "data/run1/K9013rcc-rxn3_S423_L003_R2_001.fastq.gz"
  arguments$width <- 7
  arguments$height <- 3.2

}

library(magrittr)
library(tidyverse)
library(dada2)
library(cowplot)

r1_file <- arguments$end1
r2_file <- arguments$end2

qc_plot <- cowplot::plot_grid(
  dada2::plotQualityProfile(r1_file),
  dada2::plotQualityProfile(r2_file), nrow = 1)

ggsave(
  filename = arguments$outfile,
  plot = qc_plot,
  width = as.numeric(arguments$width),
  height = as.numeric(arguments$height))
