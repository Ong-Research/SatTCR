#-------------------------------------------------------------------------------
# 
#   Snakefile to assemble clonotypes from repertoire data / RNA-seq using
#    TRUST4
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#-------------------------------------------------------------------------------

rule get_imgt_annot:
  """
  Gets a species reference based on IMGT, the reference contains all
  IG and TR references
  """
  output:
    "output/trust4/imgt_annots/{species}/imgt.fa"
  params:
    species = config["species"]
  log: "logs/trust4/get_imgt_annot_{species}.log"
  shell:
    """BuildImgtAnnot.pl {params.species} > {output}"""

rule align_imgt:
  """
  Aligns the IMGT index to the reference genome using bowtie2
  """
  input:
    fa = "output/trust4/imgt_annots/{species}/imgt.fa"
  output:
    align = "output/trust4/imgt_annots/{species}.bam",
    index = "output/trust4/imgt_annots/{species}.bam.bai"
  params:
    idx = config["index"],
  threads: 16
  log: "logs/trust4/align_imgt_annot_{species}.log"
  shell:
    """
    bowtie2 --no-unal --very-fast --very-fast-local -p {threads} -x {params.idx} \
      -f {input.fa} | samtools view -bS - |
      samtools sort - -o {output.align}
    samtools index {output.align}
    """


rule generate_fa_ref:
  """
  generates the references for 
  """
  input:
    fa = "output/trust4/imgt_annots/{species}/imgt.fa",
    align = "output/trust4/imgt_annots/{species}.bam"
  output:
    ref = "output/trust4/imgt_annots/{species}/imgt_ref.fa"
  params:
    idx = config["index"]
  threads: 1
  log: "logs/trust4/create_imgt_ref_{species}.log"
  script:
    """../scripts/sequences/generate_fasta_ref.R"""

rule subset_fastq:
  """
  Subset fastq files to the reference alignments
  """
  input:
    sort_bam="output/aligned/{sample}.bam",
    r1 = "output/trimmed/{sample}_R1.fastq.gz" if config["trim_sequences"] else lambda wc: sample_dict[wc.sample]["end1"],
    r2 = "output/trimmed/{sample}_R2.fastq.gz" if config["trim_sequences"] else lambda wc: sample_dict[wc.sample]["end2"],
    imgt=expand("output/trust4/imgt_annots/{species}.bam",
      species = config["species"])
  output:
    r1="output/fastq_corrected/{sample}_R1.fastq.gz",
    r2="output/fastq_corrected/{sample}_R2.fastq.gz"
  params:
    window = config["window"]
  threads: 1
  script:
    """../scripts/sequences/subset_fastq.R"""

rule process_sample_trust4:
  """
  Process a sample with trust4
  """
  input:
    r1="output/fastq_corrected/{sample}_R1.fastq.gz",
    r2="output/fastq_corrected/{sample}_R2.fastq.gz",
    imgt="output/trust4/imgt_annots/{species}/imgt.fa"
  output:
    outdir=directory("output/clonotypes/trust4/{species}/{sample}/"),
    files = multiext(
      "output/clonotypes/trust4/{species}/{sample}/TRUST_{sample}_R1_",
      "airr_align.tsv", "final.out", "airr.tsv", "raw.out", "annot.fa",
      "report.tsv", "assembled_reads.fa", "toassemble_1.fq",
      "barcode_airr.tsv", "toassemble_2.fq", "barcode_report.tsv",
      "toassemble_bc.fa", "cdr3.out")
  params:
    barcode=config["trust4"]["barcode_range"],
    r1=config["trust4"]["read1_range"],
    r2=config["trust4"]["read2_range"]
  threads: 32
  shell:
    """
    run-trust4 --ref {input.imgt} -f {input.imgt} \
      -1 {input.r1} -2 {input.r2} --barcode {input.r1} \
      --barcodeRange 0 {params.barcode} + \
      --read1Range {params.r1} -1 \
      --read2Range {params.r2} -1 \
      -t {threads} --od {output.outdir}
    """

rule filter_simple_repertoire:
  """
  Formats a simple report based on TRUST4's `trust-simplerep.pl` function
  """
  input:
    cdr3 = "output/clonotypes/trust4/{species}/{sample}/TRUST_{sample}_R1_cdr3.out"
  output:
    report = "output/clonotypes/trust4/{species}/TRUST_{sample}_simple_report.tsv"
  params:
    tcr_error = config["trust4"]["filter_tcr_error"]
  threads: 1
  shell:
    """
    trust-simplerep.pl {input.cdr3} \
      --filterTcrError {params.tcr_error} > {output.report}
    """