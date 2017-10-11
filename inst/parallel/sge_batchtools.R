cluster.functions <- makeClusterFunctionsSGE(system.file("parallel/sge_script.tmpl", package = "RMINC"))
default.resources <-
  list(nodes = 1,
       vmem = "8G",
       walltime = "01:00:00")

