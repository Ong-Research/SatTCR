#------------------------------------------------------------------------------#
# 
#   Quality control rules
#
#     Author: Rene Welch rwelch2@wisc.edu
#     Date: 2024-02-09
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
  params:
    sample = "{sample}",
    prefix1 = lambda wc: os.path.basename(sample_dict[wc.sample]["end1"].replace(".fastq.gz", "")),
    prefix2 = lambda wc: os.path.basename(sample_dict[wc.sample]["end2"].replace(".fastq.gz", "")),
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["fastqc"]
  log: "logs/qc/{sample}_fastqc.txt"
  shell:
    """{params.docker_run} {params.image} \
    fastqc -o output/qc/fastqc -t {threads} {input.r1} {input.r2} > {log}
    mv output/qc/fastqc/{params.prefix1}_fastqc.zip output/qc/fastqc/{params.sample}_R1_fastqc.zip
    mv output/qc/fastqc/{params.prefix2}_fastqc.zip output/qc/fastqc/{params.sample}_R2_fastqc.zip
    mv output/qc/fastqc/{params.prefix1}_fastqc.html output/qc/fastqc/{params.sample}_R1_fastqc.html
    mv output/qc/fastqc/{params.prefix2}_fastqc.html output/qc/fastqc/{params.sample}_R2_fastqc.html
    """

rule multiqc:
  """
  Builds a MultiQC report based on the FASTQC ones
  """  
  input:
    expand("output/qc/fastqc/{sample}_R1_fastqc.zip", sample = sample_names),
    expand("output/qc/fastqc/{sample}_R2_fastqc.zip", sample = sample_names),
  output:
    "output/qc/multiqc/multiqc_report.html"
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["multiqc"]
  shell:
    """{params.docker_run} {params.image} \
      multiqc -f output/qc/fastqc -o output/qc/multiqc/"""

rule dada2_qc_profiles:
  """
  Plots QC profiles with dada2
  """
  input:
    r1 = lambda wc: sample_dict[wc.sample]["end1"],
    r2 = lambda wc: sample_dict[wc.sample]["end2"],
  output:
    "output/qc/figs/{sample}_qc_profile.png"
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
    width = 7,
    height = 3.2,
    units = "in"
  threads: 1
  shell:
    """{params.docker_run} {params.image} \
      Rscript workflow/scripts/qc/plot_qc_profiles.R {output} \
        --end1={input.r1} --end2={input.r2} \
        --width={params.width} --height={params.height} > {log}"""

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
      # docker parameters
      docker_run = config["docker"]["run_line"],
      image = config["docker"]["trimmomatic"],
      # list of trimmers (see manual)
      trimmer=" ".join(config["trim"]["trimmer"]),
      compression_level="-9"
    threads: config["threads"]
    resources:
      mem_mb = 1024
    shell:
      """{params.docker_run} {params.image} \
      trimmomatic PE {input.r1} {input.r2} \
      {output.r1} {output.r1_unpaired} \
      {output.r2} {output.r2_unpaired} \
      {params.trimmer}"""


rule count_sequences:
  """
  Count the number of sequences in the raw and trimmed sequenced file
  """
  input:
    raw = lambda wc: sample_dict[wc.sample]["end1"],
    trimmed = "output/trimmed/{sample}_R1.fastq.gz"
  output:
    summary = "output/qc/{sample}_trimmed.tsv"
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
    sample = "{sample}"
  shell:
    """{params.docker_run} {params.image} \
      Rscript workflow/scripts/qc/count_sequences.R {output} \
        --input={input.raw} --trimmed={input.trimmed} \
        --sample={params.sample}"""
