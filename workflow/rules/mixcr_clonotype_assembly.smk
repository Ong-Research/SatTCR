#-------------------------------------------------------------------------------
# 
#   Snakefile to assemble clonotypes from repertoire data / RNA-seq using
#    MIXCR
#
#   Note:
#     MIXCR requires a license from
#     https://docs.milaboratories.com/mixcr/getting-started/license/
#     to work
#
#     Author: Rene Welch rwelch2@wisc.edu
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
    outdir=directory("output/clonotypes/mixcr/{species}/{sample}/"),
    files = multiext("output/clonotypes/mixcr/{species}/{sample}/{sample}",
      ".clones_TRB.tsv",
      ".assemble.report.json",
      ".assemble.report.txt",
      ".clns",
      ".extended.vdjca",
      ".extend.report.json",
      ".extend.report.txt",
      ".passembled.2.vdjca",
      ".assemblePartial.report.json",
      ".assemblePartial.report.txt",
      ".passembled.1.vdjca",
      ".vdjca",
      ".align.report.json",
      ".align.report.txt")
  params:
    mixcr=config["mixcr"]["params"],
    sample="{sample}"
  threads: 12
  resources:
    mem_mb = lambda wildcards, threads: 1200 * threads
  log: "logs/mixcr/{species}/{sample}_run.log"
  shell:
    """
    mixcr analyze {params.mixcr} \
      -t {threads} \
      {input.r1} {input.r2} \
      {output.outdir}/{params.sample}
    """


rule process_saturation_mixcr:
  """
  Process a sample with mixcr
  """
  input:
    r1="output/saturation/{seed}/{perc}/fastq/{sample}_R1.fastq.gz",
    r2="output/saturation/{seed}/{perc}/fastq/{sample}_R2.fastq.gz",
  output:
    outdir=directory("output/saturation/{seed}/{perc}/mixcr/{species}/{sample}/"),
    file = "output/saturation/{seed}/{perc}/mixcr/{species}/{sample}/{sample}.assemble.report.json"
  params:
    mixcr=config["mixcr"]["params"],
    sample="{sample}"
  threads: 12
  resources:
    mem_bm = lambda wildcards, threads: 200 * threads
  log: "logs/saturation/mixcr/{species}/{seed}/{perc}/{sample}_run.log"
  shell:
    """
    mixcr analyze {params.mixcr} \
      -t {threads} \
      {input.r1} {input.r2} \
      {output.outdir}/{params.sample}
    """


