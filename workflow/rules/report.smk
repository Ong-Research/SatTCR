#-------------------------------------------------------------------------------
# 
#   Snakefile with rules to generate a quarto bookdown report
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

quarto_vec = ["index", "qc", "repertoire", "saturation"] 
quarto_template_vec = ["template_upset", "template_qc_profile", "template_chord",
  "appendix"]

rule report:
  """
  Generates a report with summarized results of the analysis
  """
  input:
    expand("workflow/templates/{file}.qmd",
      file = quarto_vec)
  output:
    expand("output/report/{file}.html",
      file = quarto_vec)
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
  threads: config["threads"]
  shell:
    """
    cp output/qc/multiqc/multiqc_report.html output/report/
    cp -r output/qc/figs output/report/
    {params.docker_run} {params.image} \
    quarto render workflow/templates --no-cache
    """
  
qc_vec = ["index", "qc"]
rule qc_report:
  """
  Generates a report with the QC of the analysis
  """
  input:
    expand("workflow/templates/{file}.qmd",
      file = qc_vec)
  output:
    expand("output/report/{file}.html",
      file = qc_vec)
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
  threads: config["threads"]
  shell:
    """
    cp output/qc/multiqc/multiqc_report.html output/report/
    cp -r output/qc/figs output/report/
    {params.docker_run} {params.image} \
    quarto render workflow/templates/index.qmd --no-cache
    {params.docker_run} {params.image} \
    quarto render workflow/templates/qc.qmd --no-cache
    """  

quarto_vec = ["index", "qc", "repertoire"]
rule repertoire_report:
  """
  Generates a report with the QC of the analysis
  """
  input:
    expand("workflow/templates/{file}.qmd",
      file = qc_vec)
  output:
    expand("output/report/{file}.html",
      file = qc_vec)
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
  threads: config["threads"]
  shell:
    """
    cp output/qc/multiqc/multiqc_report.html output/report/
    cp -r output/qc/figs output/report/
    {params.docker_run} {params.image} \
      quarto render workflow/templates/index.qmd --no-cache
    {params.docker_run} {params.image} \
      quarto render workflow/templates/qc.qmd --no-cache
    {params.docker_run} {params.image} \
    quarto render workflow/templates/repertoire.qmd --no-cache
    """  