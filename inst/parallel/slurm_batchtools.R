cluster.functions <-
  makeClusterFunctionsSlurm(system.file("parallel/slurm_script.tmpl"
                                      , package = "RMINC"))

default.resources <-
  list(memory = "8G"
     , walltime = 120
     , ntasks = 1
     , ncpus = 1)
