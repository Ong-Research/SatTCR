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

rule tcr_stats:
  """
  Calculate TCR statistics for MIXCR results
  """
  input:
    airr = "output/clonotypes/mixcr/{sample}/{sample}_airr.tsv"
  output:
    stats = "output/clonotypes/mixcr/{sample}/{sample}_summary.tsv"
  threads: 1
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
    min_count = config["summary"]["min_count"]
  shell:
    """
    {params.docker_run} {params.image} \
      Rscript workflow/scripts/repertoire/calculate_tcr_stats.R \
         {output.stats} --airr={input.airr} --min_count={params.min_count}
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

rule calculate_tcr_stats_saturation:
  """
  Calculate TCR statistics for MIXCR saturation results, and calculates the 
  overlap between the clonotypes assembled with the whole 
  sample vs the ones assembled with a subset of all the reads
  """
  input:
    airr = "output/seq_bootstrap/mixcr/{seed}/{sample}/{subsample}/{subsample}_airr.tsv",
    airr_full = "output/clonotypes/mixcr/{sample}/{sample}_airr.tsv"
  output:
    stats = "output/seq_bootstrap/mixcr/{seed}/{sample}/{subsample}/{subsample}_summary.tsv",
    overlap = "output/seq_bootstrap/mixcr/{seed}/{sample}/{subsample}/{subsample}_overlap.tsv"    
  threads: 1
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
    min_count = config["summary"]["min_count"]
  shell:
    """
    {params.docker_run} {params.image} \
      Rscript workflow/scripts/repertoire/calculate_tcr_stats.R \
         {output.stats} --airr={input.airr} --min_count={params.min_count}

    {params.docker_run} {params.image} \
      Rscript workflow/scripts/repertoire/calculate_tcr_overlap.R \
        {output.overlap} \
        --airr1={input.airr_full} --airr2={input.airr} \
        --min_count={params.min_count}
    """
