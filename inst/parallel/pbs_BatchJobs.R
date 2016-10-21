cluster.functions <- makeClusterFunctionsTorque(system.file("parallel/pbs_script.tmpl", package = "RMINC"))
default.resources <-
  list(nodes = 1,
       vmem = "8G",
       walltime = "01:00:00")

