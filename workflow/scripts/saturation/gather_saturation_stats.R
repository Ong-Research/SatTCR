
library(magrittr)
library(tidyverse)
library(future)
library(furrr)

plan(multicore, workers = snakemake@threads)

conflicted::conflict_prefer_all("dplyr")
conflicted::conflicts_prefer(purrr::set_names)


out <- tibble(file = unlist(snakemake@input))

out %<>%
  mutate(
    summary = future_map(file, qs::qread),
    sample = basename(file),
    sample = str_remove(sample, "_stats.qs")) %>%
  select(sample, summary)

out %<>%
  separate(sample, into = c("run", "bio", "tech", "boots"),
    sep = "\\_", remove = FALSE)

out %>%
  qs::qsave(snakemake@output[["stats"]])
