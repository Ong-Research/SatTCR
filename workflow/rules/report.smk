#-------------------------------------------------------------------------------
# 
#   Snakefile with rules to generate a quarto bookdown report
#
#     Author: Rene Welch rwelch2@wisc.edu
#
#------------------------------------------------------------------------------#

rule quarto_report:
  """
  Generates a report with summarized results of the analysis
  """
  input:
    expand("workflow/templates/{file}.qmd",
      file = ["index", "qc", "repertoire", "saturation"])
  output:
    expand("output/report/{file}.html",
      file = ["index", "qc", "repertoire", "saturation"])
  params:
    docker_run = config["docker"]["run_line"],
    image = config["docker"]["rquarto"],
  threads: config["threads"]
  shell:
    """{params.docker_run} {params.image} \
    quarto render --no-cache workflow/templates"""
  