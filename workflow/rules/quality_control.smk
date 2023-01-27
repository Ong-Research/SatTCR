#------------------------------------------------------------------------------#
# 
#   Quality control rules
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule fastqc:
  """
  Process the quality control of a paired end 
  """
  input:
    r1= lambda wc: sample_dict[wc.sample]["end1"],
    r2= lambda wc: sample_dict[wc.sample]["end2"],
  threads: 1
  output:
    multiext("output/qc/fastqc/{sample}_",
      "R1_fastqc.zip", "R2_fastqc.zip",
      "R1_fastqc.html", "R2_fastqc.html")
  log: "logs/qc/{sample}_fastqc.txt"
  shell:
    """fastqc -o output/qc/fastqc -t {threads} {input.r1} {input.r2}"""
    
rule multiqc:
  input:
    expand("output/qc/fastqc/{sample}_R1_fastqc.zip", sample = sample_names),
    expand("output/qc/fastqc/{sample}_R2_fastqc.zip", sample = sample_names),
  output:
    "output/qc/multiqc/multiqc_report.html"
  shell:
    """multiqc output/qc/fastqc -o output/qc/multiqc"""