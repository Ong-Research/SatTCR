#------------------------------------------------------------------------------#
# 
#   Rules to gather the results from assembling the clonotypes
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule gather_trust4_report:
  """
  Gather the TRUST4 reports into a unique file
  """
  input:
    expand("output/clonotypes/trust4/{specie}/{sample}/TRUST_{sample}_R1_report.tsv", sample = sample_names, specie = config["species"])
  output:
    "output/gather/trust4/gathered_report.qs"
  params:
    method = "trust4"
  script:
    "../scripts/gather_results/gather_tsv_files.R"

rule gather_trust4_saturation:
  """
  Gather the TRUST4 reports into a unique file
  """
  input:
    expand("output/saturation/{seed}/{perc}/trust4/{specie}/{sample}/TRUST_{sample}_R1_report.tsv",
      seed = seeds, perc = saturation_frequencies,
      specie = config["species"], sample = sample_names)
  output:
    "output/gather/trust4/gathered_report_saturation.qs"
  params:
    method = "trust4"
  script:
    "../scripts/gather_results/gather_tsv_files.R"

rule gather_mixcr_report:
  """
  Gather the MIXCR reports into a unique file
  """
  input:
    expand("output/clonotypes/mixcr/{specie}/{sample}/{sample}.clones_TRB.tsv", sample = sample_names, specie = config["species"])
  output:
    "output/gather/mixcr/gathered_report.qs"
  params:
    method = "mixcr"
  script:
    "../scripts/gather_results/gather_tsv_files.R"

rule gather_mixcr_saturation:
  """
  Gather the MIXCR reports into a unique file
  """
  input:
    expand("output/saturation/{seed}/{perc}/mixcr/{specie}/{sample}/{sample}.clones_TRB.tsv",
      seed = seeds, perc = saturation_frequencies,
      specie = config["species"], sample = sample_names)
  output:
    "output/gather/mixcr/gathered_report_saturation.qs"
  params:
    method = "mixcr"
  script:
    "../scripts/gather_results/gather_tsv_files.R"
