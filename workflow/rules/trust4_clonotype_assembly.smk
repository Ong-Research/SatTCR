#-------------------------------------------------------------------------------
# 
#   Snakefile to assemble clonotypes from repertoire data / RNA-seq using
#    TRUST4
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule get_imgt_annot:
  """
  Gets a species reference based on IMGT
  """
  output:
    "output/trust4/imgt_annots/{species}/imgt.fa"
  params:
    species = config["trust4"]["species"]
  shell:
    """BuildImgtAnnot.pl {params.species} > {output}"""

rule process_sample_trust4:
  """
  Process a sample with trust4
  """
  input:
    r1="output/fastq/{perc}/{sample}_R1_001.fastq.gz",
    r2="output/fastq/{perc}/{sample}_R2_001.fastq.gz",
    imgt="output/trust4/imgt_annots/{species}/imgt.fa"
  output:
    outdir=directory("output/trust4/{species}/{perc}/{sample}/")
  params:
    barcode=config["trust4"]["barcode_range"],
    r1=config["trust4"]["read1_range"],
    r2=config["trust4"]["read2_range"]
  threads: config["threads"]
  shell:
    """
    run-trust4 --ref {input.imgt} -f {input.imgt} \
      -1 {input.r1} -2 {input.r2} --barcode {input.r1} \
      --barcodeRange 0 {params.barcode} + \
      --read1Range {params.r1} -1 \
      --read2Range {params.r2} -1 \
      -t {threads} --od {output.outdir}
    """
