#-------------------------------------------------------------------------------
# 
#   Snakefile rule to assemble clonotypes from repertoire data / RNA-seq using
#    MIXCR
#
#   Note:
#     MIXCR requires a license from
#     https://docs.milaboratories.com/mixcr/getting-started/license/
#     to work
#
#     Author: Rene Welch rwelch2@wisc.edu
#     Date: 2024-02-12
#
#-------------------------------------------------------------------------------


rule process_sample_mixcr:
  """
  Process a sample with mixcr
  """
  input:
    r1 = "output/trimmed/{sample}_R1.fastq.gz",
    r2 = "output/trimmed/{sample}_R2.fastq.gz",
  output:
    outdir=directory("output/clonotypes/mixcr/{sample}/"),
    clono = "output/clonotypes/mixcr/{sample}/{sample}.contigs.clns",
    airr = "output/clonotypes/mixcr/{sample}/{sample}_airr.tsv"
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["mixcr"],
    mixcr_license_file = config["mixcr"]["license_file"],
    mixcr=config["mixcr"]["params"],
    sample="{sample}"
  threads: 12
  resources:
    mem_mb = lambda wildcards, threads: 1200 * threads
  log: "logs/mixcr/{sample}_run.log"
  shell:
    """
    {params.docker_run} \
    -v ./license_mixcr:/opt/mixcr/mi.license:ro \
    {params.image} \
    mixcr analyze {params.mixcr} \
      -t {threads} -f \
      {input.r1} {input.r2} \
      {output.outdir}/{params.sample}
    {params.docker_run} \
      -v ./license_mixcr:/opt/mixcr/mi.license:ro \
      {params.image} \
    mixcr exportAirr -f {output.clono} {output.airr}
    """



rule process_saturation_mixcr:
  """
  Process a saturation sample with mixcr
  """
  input:
    r1 = "output/seq_bootstrap/{seed}/{sample}/{subsample}_R1.fastq.gz",
    r2 = "output/seq_bootstrap/{seed}/{sample}/{subsample}_R2.fastq.gz"
  output:
    outdir = directory("output/seq_bootstrap/mixcr/{seed}/{sample}/{subsample}"),
    clns = "output/seq_bootstrap/mixcr/{seed}/{sample}/{subsample}/{subsample}.contigs.clns",
    airr = "output/seq_bootstrap/mixcr/{seed}/{sample}/{subsample}/{subsample}_airr.tsv"
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["mixcr"],
    mixcr_license_file = config["mixcr"]["license_file"],
    mixcr=config["mixcr"]["params"],
    subsample="{subsample}"
  threads: 4
  resources:
    mem_mb = lambda wildcards, threads: 1200 * threads
  shell:
    """
    {params.docker_run} \
    -v ./license_mixcr:/opt/mixcr/mi.license:ro \
    {params.image} \
    mixcr analyze {params.mixcr} \
      -t {threads} \
      {input.r1} {input.r2} \
      {output.outdir}/{params.subsample}
    {params.docker_run} \
      -v ./license_mixcr:/opt/mixcr/mi.license:ro \
      {params.image} \
    mixcr exportAirr {output.clns} {output.airr}
    """


rule gather_results_saturation_mixcr:
  input:
    airr = "output/clonotypes/mixcr/{species}/{sample}/{sample}_airr.tsv",
    summary = "output/seq_bootstrap/{seed}/{sample}/bootstrap_summary.qs",
  output:
    out = "output/seq_bootstrap/results/mixcr/{species}/{seed}/{sample}.qs"
  params:
    sample = "{sample}",
    bootdir = "output/seq_bootstrap/mixcr/{seed}/{sample}"
  threads: 24
  script:
    "../scripts/gather_results/gather_mixcr_airr.R"

rule gather_results_saturation_mixcr_perc:
  input:
    airr = "output/clonotypes/mixcr/{species}/{sample}/{sample}_airr.tsv",
    summary = "output/seq_bootstrap_perc/{seed}/{sample}/bootstrap_summary.qs",
  output:
    out = "output/seq_bootstrap_perc/results/mixcr/{species}/{seed}/{sample}.qs"
  params:
    sample = "{sample}",
    bootdir = "output/seq_bootstrap_perc/mixcr/{seed}/{sample}"
  threads: 24
  script:
    "../scripts/gather_results/gather_mixcr_airr.R"
