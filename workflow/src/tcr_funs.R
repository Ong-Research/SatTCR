
# Clean an AIRR file by removing columns that are not
# used when assembling TCRB chains with MIXCR
clean_clono <- function(clono) {

  clono %<>%
    select(-contains("alignment"),
      -contains("germline"), -ends_with("cigar"),
      -starts_with("fwr"), -ends_with("end"), -ends_with("start"),
      -starts_with("cdr1"), -starts_with("cdr2"), -starts_with("np1"),
      -starts_with("np2"), -rev_comp, -complete_vdj)

  clono %>%
    select(-c_score) %>%
    rename(complete = productive, count = duplicate_count)

}

# Compute the D50 index
calculate_d50 <- function(data) {

  data %>%
    group_by(cdr3, v_call, j_call) %>%
    summarize(
      count = sum(count), .groups = "drop") %>%
    arrange(desc(count)) %>%
    mutate(
      id = row_number(),
      acc_count = cumsum(count) / sum(count)) %>%
    filter(acc_count >= 0.5) %>%
    pull(id) %>%
    min()

}

## Computes Shannon diversity
calculate_shannon <- function(data) {

  data %<>%
    group_by(cdr3, v_call, j_call) %>%
    summarize(
      count = sum(count), .groups = "drop") %>%
    arrange(desc(count)) %>%
    mutate(
      p = count / sum(count))

  p <- data[["p"]]
  - sum(p * log(p)) / log(nrow(data))

}


summarize_clonotypes <- function(data) {

  data %>%
    summarize(
      nclono = n(),
      shannon = calculate_shannon(.),
      d50 = calculate_d50(.)) %>%
    mutate(clonality = 1 - shannon)
}
