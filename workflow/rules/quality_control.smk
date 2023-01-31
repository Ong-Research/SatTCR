#------------------------------------------------------------------------------#
# 
#   Quality control rules
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule fastqc:
  """
  Process the quality control of a paired end 
  """
  input:
    r1 = lambda wc: sample_dict[wc.sample]["end1"],
    r2 = lambda wc: sample_dict[wc.sample]["end2"],
  threads: 1
  output:
    multiext("output/qc/fastqc/{sample}_",
      "R1_fastqc.zip", "R2_fastqc.zip",
      "R1_fastqc.html", "R2_fastqc.html")
  log: "logs/qc/{sample}_fastqc.txt"
  shell:
    """fastqc -o output/qc/fastqc -t {threads} {input.r1} {input.r2}"""
    
rule multiqc:
  input:
    expand("output/qc/fastqc/{sample}_R1_fastqc.zip", sample = sample_names),
    expand("output/qc/fastqc/{sample}_R2_fastqc.zip", sample = sample_names),
  output:
    "output/qc/multiqc/multiqc_report.html"
  shell:
    """multiqc output/qc/fastqc -o output/qc/multiqc"""

rule dada2_qc_profiles:
  """
  Plots QC profiles with dada2
  """
  input:
    r1 = lambda wc: sample_dict[wc.sample]["end1"],
    r2 = lambda wc: sample_dict[wc.sample]["end2"],
  output:
    plot = "output/qc/figs/{sample}_qc_profile.png"
  params:
    width = 7,
    height = 3.2,
    units = "in"
  threads: 1
  script:
    """../scripts/qc/plot_qc_profiles.R"""

rule trimmomatic_pe:
    input:
      r1 = lambda wc: sample_dict[wc.sample]["end1"],
      r2 = lambda wc: sample_dict[wc.sample]["end2"],
    output:
      r1 = "output/trimmed/{sample}_R1.fastq.gz",
      r2 = "output/trimmed/{sample}_R2.fastq.gz",
      # reads where trimming entirely removed the mate
      r1_unpaired = "output/trimmed/{sample}_R1.unpaired.fastq.gz",
      r2_unpaired = "output/trimmed/{sample}_R2.unpaired.fastq.gz"
    log: "logs/trimmomatic/{sample}.log"
    params:
      # list of trimmers (see manual)
      trimmer=config["qc"]["trimmer"],
      # optional parameters
      extra="",
      compression_level="-9"
    threads: config["threads"]
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
      mem_mb = 1024
    wrapper:
      "v1.21.6/bio/trimmomatic/pe"

rule count_sequences:
  """
  Count the number of sequences in the main and trimmed sequenced file
  """
  input:
    original = lambda wc: sample_dict[wc.sample]["end1"],
    trimmed = "output/trimmed/{sample}_R1.fastq.gz"
  output:
    summary = "output/qc/{sample}_trimmed.tsv"
  params:
    sample = "{sample}"
  script:
    """../scripts/qc/count_sequences.R"""

