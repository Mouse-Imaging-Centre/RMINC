cluster.functions <- makeClusterFunctionsTorque(system.file("parallel/sge_script.tmpl", package = "RMINC"))
default.resources <-
  list(nodes = 1,
       memory = "8G",
       walltime = "01:00:00",
       modules = NULL,
       module_paths = NULL)
db.options <- list(pragmas = c("busy_timeout=5000", "journal_mode=WAL"))
fs.timeout <- 10

