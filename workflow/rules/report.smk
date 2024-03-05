#-------------------------------------------------------------------------------
# 
#   Snakefile with rules to generate a quarto bookdown report
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule report:
  """
  Generates a report with summarized results of the analysis
  """
  input:
    expand("workflow/templates/{file}.qmd",
      file = ["index", "qc", "repertoire", "saturation", "template_upset"])
  output:
    expand("output/report/{file}.html",
      file = ["index", "qc", "repertoire", "saturation"])
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
  

  
