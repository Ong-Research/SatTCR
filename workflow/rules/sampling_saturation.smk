#-------------------------------------------------------------------------------
# 
#   Snakefile to perform saturation / rarefaction of fastq sequence files
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule sequential_saturation:
  """
  Samples sequences from the trimmed fastq files that were subsetted to overlap the
  TRB loci
  """
  input:
    r1 = "output/trimmed/{sample}_R1.fastq.gz",
    r2 = "output/trimmed/{sample}_R2.fastq.gz"
  output:
    outdir = directory("output/seq_bootstrap/{seed}/{sample}/"),
    summary_file = "output/seq_bootstrap/{seed}/{sample}/bootstrap_summary.qs"
  threads: 1
  log: "logs/seq_bootstrap/{sample}_{seed}.log"
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
    nreps = config["saturation"]["bootstrap_replicates"],
    nblocks = config["saturation"]["nblocks"],
    block_size = config["saturation"]["block_size"],
    seed = "{seed}",
    sample = "{sample}"
  shell:
    """{params.docker_run} {params.image} \
      Rscript workflow/scripts/sequences/sample_sequences.R \
      {output.outdir} --summary={output.summary_file} \
      --r1={input.r1} --r2={input.r2} \
      --block_size={params.block_size} --nblocks={params.nblocks} \
      --nreplicates={params.nreps} --seed={params.seed} --sample={params.sample}
    """

rule number_of_sequences_seq_saturation:
  """
  Count the number of sequences in the raw and trimmed sequenced file
  """
  input:
    raw = "output/seq_bootstrap/{seed}/{sample}/{subsample}_R1.fastq.gz",
  output:
    summary = "output/seq_bootstrap/nseqs/{seed}/{sample}/{subsample}_nseqs.tsv"
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],  
    sample = "{sample}",
  threads: 1
  shell:
    """{params.docker_run} {params.image} \
      Rscript workflow/scripts/qc/count_sequences.R {output} \
        --input={input.raw} \
        --sample={params.sample}"""
