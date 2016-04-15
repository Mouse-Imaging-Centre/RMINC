cluster.functions <- makeClusterFunctionsTorque(system.file("parallel/pbs_script.tmpl", package = "RMINC"))
default.resources <-
  list(nodes = 1,
       vmem = "2G",
       walltime = "01:00:00",
       modules = NULL,
       module_paths = NULL)

