#-------------------------------------------------------------------------------
# 
#   Snakefile to assemble clonotypes from repertoire data / RNA-seq using
#    TRUST4 and bam files
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#-------------------------------------------------------------------------------

rule align_imgt_bowtie:
  """
  Aligns the IMGT index to the reference genome using bowtie2
  """
  input:
    fa = "output/trust4/imgt_annots/{species}/imgt.fa"
  output:
    align = "output/trust4_bam/imgt_annots/{species}.bam",
    index = "output/trust4_bam/imgt_annots/{species}.bam.bai"
  params:
    idx = config["index"],
  threads: 16
  log: "logs/trust4_bam/align_imgt_annot_{species}.log"
  shell:
    """
    bowtie2 --no-unal --very-sensitive --very-sensitive-local -p {threads} -x {params.idx} \
      -f {input.fa} | samtools view -bS - |
      samtools sort - -o {output.align}
    samtools index {output.align}
    """

rule build_trust4_db:
  """
  Builds a fasta reference with GenomicAlignments to be used with TRUST4
  """
  input:
    imgt = "output/trust4/imgt_annots/{species}/imgt.fa",
    align = "output/trust4_bam/imgt_annots/{species}.bam",
    index = "output/trust4_bam/imgt_annots/{species}.bam.bai",
    gtf = "data/refs/ensembl/genes/Canis_lupus_familiaris.ROS_Cfam_1.0.109.gtf.gz"
  output:
    trust4_fa = "output/trust4_bam/imgt_annots/{species}_bcrtcr.fa",
    name_file = "output/trust4_bam/imgt_annots/{species}_names.list"
  threads: 16
  script:
    """../scripts/optional_files/build_db_custom.R"""

rule align_samples_ownref:
  """
  Aligns the trust4 samples to our new reference
  """
  input:
    r1 = "output/trimmed/{sample}_R1.fastq.gz" if config["trim_sequences"] else lambda wc: sample_dict[wc.sample]["end1"],
    r2 = "output/trimmed/{sample}_R2.fastq.gz" if config["trim_sequences"] else lambda wc: sample_dict[wc.sample]["end2"],
  output:
    align = "output/trust4_bam/aligned/{sample}.bam",
    index = "output/trust4_bam/aligned/{sample}.bam.bai"
  params:
    idx = config["index"]
  threads: 16
  log: "logs/alignments/align_{sample}.log"
  shell:
    """
    bowtie2 --sensitive-local -p {threads} -x {params.idx} \
      -1 {input.r1} -2 {input.r2} | samtools view -bS - |
      samtools sort - -o {output.align} 2> {log}
    samtools index {output.align}
    """

rule process_sample_trust4_bam:
  """
  Process a sample with trust4
  """
  input:
    bam = "output/trust4_bam/aligned/{sample}.bam",
    index = "output/trust4_bam/aligned/{sample}.bam.bai",
    bcrtcr = "output/trust4_bam/imgt_annots/{species}_bcrtcr.fa",
    names = "output/trust4_bam/imgt_annots/{species}_names.list",
    imgt ="output/trust4/imgt_annots/{species}/imgt.fa"
  output:
    outdir=directory("output/clonotypes/trust4_bam/{species}/{sample}/"),
    files = multiext(
      "output/clonotypes/trust4_bam/{species}/{sample}/TRUST_{sample}_R1_",
      "airr_align.tsv", "final.out", "airr.tsv", "raw.out", "annot.fa",
      "report.tsv", "assembled_reads.fa", "toassemble_1.fq",
      "barcode_airr.tsv", "toassemble_2.fq", "barcode_report.tsv",
      "toassemble_bc.fa", "cdr3.out")
  params:
    barcode=config["trust4"]["barcode_range"],
    r1=config["trust4"]["read1_range"],
    r2=config["trust4"]["read2_range"]
  threads: config["threads"]
  log: "logs/trust4_bam/{species}/{sample}_run.log"
  shell:
    """
    run-trust4 \
      --ref {input.imgt} -f {input.bcrtcr} \
      -b {input.bam} -t {threads} --od {output.outdir}
    """
