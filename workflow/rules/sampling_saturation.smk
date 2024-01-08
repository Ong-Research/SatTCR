#-------------------------------------------------------------------------------
# 
#   Snakefile to perform saturation / rarefaction of fastq sequence files
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule sample_sequences:
  """
  Sample sequences from the fastq file that were subsetted to overlapping the
  TRB loci
  """
  input:
    r1 = "output/trimmed/{sample}_R1.fastq.gz",
    r2 = "output/trimmed/{sample}_R2.fastq.gz"
  output:
    r1 = "output/saturation/{seed}/{perc}/fastq/{sample}_R1.fastq.gz",
    r2 = "output/saturation/{seed}/{perc}/fastq/{sample}_R2.fastq.gz"
  log: "logs/saturation/subsample/{sample}_{perc}_{seed}.log"
  params:
    sample_perc="{perc}",
    seed="{seed}"
  script:
    "../scripts/sequences/sample_sequences.R"

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
    seed = "{seed}",
    sample = "{sample}"
  script:
    "../scripts/sequences/sample_saturation.R"

rule number_of_sequences_seq_saturation:
  """
  Gets the number of sequences in the experiment
  """
  input:
    original = "output/seq_bootstrap/{seed}/{sample}/{subsample}_R1.fastq.gz",
  output:
    summary = "output/seq_bootstrap/nseqs/{seed}/{sample}/{subsample}_nseqs.tsv"
  params:
    sample = "{sample}"
  threads: 1
  script:
    """../scripts/qc/count_sequences.R"""

rule process_sample_trust4_seq_saturation:
  """
  Process a sample with trust4
  """
  input:
    r1="output/seq_bootstrap/{seed}/{sample}/{subsample}_R1.fastq.gz",
    r2="output/seq_bootstrap/{seed}/{sample}/{subsample}_R2.fastq.gz",
    imgt="output/trust4/imgt_annots/{species}/imgt.fa"
  output:    
    files = multiext(
      "output/seq_bootstrap/{seed}/clonotypes/trust4/{species}/{sample}/TRUST_{subsample}_R1_",
      "airr_align.tsv", "final.out", "airr.tsv", "raw.out", "annot.fa",
      "report.tsv", "assembled_reads.fa", "toassemble_1.fq",
      "barcode_airr.tsv", "toassemble_2.fq", "barcode_report.tsv",
      "toassemble_bc.fa", "cdr3.out")
  params:
    outdir = "output/seq_bootstrap/{seed}/clonotypes/trust4/{species}/{sample}",
    barcode=config["trust4"]["barcode_range"],
    r1=config["trust4"]["read1_range"],
    r2=config["trust4"]["read2_range"]
  threads: 16
  shell:
    """
    run-trust4 --ref {input.imgt} -f {input.imgt} \
      -1 {input.r1} -2 {input.r2} --barcode {input.r1} \
      --barcodeRange 0 {params.barcode} + \
      --read1Range {params.r1} -1 \
      --read2Range {params.r2} -1 \
      -t {threads} --od {params.outdir}
    """

rule filter_simple_repertoire_seq_saturation:
  """
  Formats a simple report based on TRUST4's `trust-simplerep.pl` function
  """
  input:
    cdr3 = "output/seq_bootstrap/{seed}/clonotypes/trust4/{species}/{sample}/TRUST_{subsample}_R1_cdr3.out"
  output:
    report = "output/seq_bootstrap/{seed}/clonotypes/trust4/{species}/{sample}/TRUST_{subsample}_simple_report.tsv"
  params:
    tcr_error = config["trust4"]["filter_tcr_error"]
  threads: 1
  shell:
    """
    trust-simplerep.pl {input.cdr3} \
      --filterTcrError {params.tcr_error} > {output.report}
    """

rule compute_stats_repertoire_seq_saturation:
  """
  Computes the sample specific stats for each saturation sample
  """
  input:
    simple_report = "output/seq_bootstrap/{seed}/clonotypes/trust4/{species}/{sample}/TRUST_{subsample}_simple_report.tsv",
    nseqs = "output/seq_bootstrap/nseqs/{seed}/{sample}/{subsample}_nseqs.tsv"
  output:
    stats = "output/seq_bootstrap/stats/{species}/{seed}/{sample}/{subsample}_stats.qs"
  threads: 4    
  script:
    """../scripts/saturation/compute_saturation_stats.R"""

seed, samples, subsamples = glob_wildcards("output/seq_bootstrap/{seed}/{sample}/{subsample}_R1.fastq.gz")
all_reports = []
for zip_words in zip(seed, samples, subsamples):
  zip_words = list(zip_words)
  seed = zip_words[0]
  sample = zip_words[1]
  subsample = zip_words[2]
  all_reports.append(
    "output/seq_bootstrap/stats/" + config["species"] + "/" +
      seed + "/" + sample + "/" + subsample + "_stats.qs")

rule gather_seq_saturation:
  """
  Gathers the saturation stats
  """
  input:
    all_reports 
  output:
    stats = "output/seq_bootstrap/summary_stats.qs"
  threads: 4
  script:
    """../scripts/saturation/gather_saturation_stats.R"""
