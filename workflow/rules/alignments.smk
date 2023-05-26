#------------------------------------------------------------------------------#
# 
#   Optional rules to pre-align the TCR data to a genome of reference
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule align_sample:
  """
  This rule aligns either a trimmed sequence file or the raw sequence file to a 
  genome of reference
  """
  input:
    r1 = "output/trimmed/{sample}_R1.fastq.gz" if config["trim_sequences"] else lambda wc: sample_dict[wc.sample]["end1"],
    r2 = "output/trimmed/{sample}_R2.fastq.gz" if config["trim_sequences"] else lambda wc: sample_dict[wc.sample]["end2"],
  output:
    align = "output/aligned/{sample}.bam",
    index = "output/aligned/{sample}.bam.bai"
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

rule get_number_aligned_reads:
  """
  Gets the number of aligned reads
  """
  input:
    align = "output/aligned/{sample}.bam",
    index = "output/aligned/{sample}.bam.bai"
  output:
    "output/aligned_reads/{sample}.tsv"
  threads: 1
  shell:
    """
    samtools view -c -F 260 {input.align} > {output}
    """
