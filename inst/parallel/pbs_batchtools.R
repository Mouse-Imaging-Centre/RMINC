cluster.functions <- makeClusterFunctionsTORQUE(system.file("parallel/pbs_script.tmpl", package = "RMINC"))
default.resources <-
  list(nodes = 1,
       memory = "8G",
       walltime = "01:00:00")

