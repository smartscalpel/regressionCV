source("R/functions.R") # defines the create_plot() function
source("R/plan.R")      # creates the drake plan
library(drake)

# Create the template file. You may have to modify it.
drake_hpc_template_file("slurm_clustermq.tmpl")

# Configure clustermq.
options(clustermq.scheduler = "slurm", template = "slurm_clustermq.tmpl")

make(
  plan, max_expand = 8,# defined in R/plan.R
  verbose = 2,cache_log_file = TRUE 
  # ,parallelism = "clustermq",
  # jobs = 8,
  #   log_make = "drake.log",
  # template = list(   memory = 12288,
  #                    log_file = "logs/slurm_%A_%a.log")
  # 

)
