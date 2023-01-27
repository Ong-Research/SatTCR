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
    r1= lambda wc: sample_dict[wc.sample]["end1"],
    r2= lambda wc: sample_dict[wc.sample]["end2"],
  output:
    outdir=directory("output/clonotypes/mixcr/{species}/{sample}/"),
    file = "output/clonotypes/mixcr/{species}/{sample}/{sample}.assemble.report.json"
  params:
    ref=config["mixcr"]["ref"],
    sample="{sample}"
  threads: 1
  log: "logs/mixcr/{species}/{sample}_run.log"
  shell:
    """
    mixcr analyze {params.ref} \
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
    ref=config["mixcr"]["ref"],
    sample="{sample}"
  threads: 1
  log: "logs/saturation/mixcr/{species}/{seed}/{perc}/{sample}_run.log"
  shell:
    """
    mixcr analyze {params.ref} \
      {input.r1} {input.r2} \
      {output.outdir}/{params.sample}
    """


