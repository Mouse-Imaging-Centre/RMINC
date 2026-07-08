# Parallel MincApply

Apply an arbitrary R function across a collection of minc files,
distributing the computation to multiple cores or workers on a grid
computing environment

## Usage

``` r
pMincApply(
  filenames,
  fun,
  ...,
  mask = NULL,
  tinyMask = FALSE,
  batches = 4,
  slab_sizes = NULL,
  method = NULL,
  local = FALSE,
  cores = NULL,
  resources = list(),
  packages = NULL,
  vmem = NULL,
  walltime = NULL,
  workers = batches,
  temp_dir = getwd(),
  cleanup = TRUE,
  collate = simplify2minc,
  conf_file = getOption("RMINC_BATCH_CONF"),
  registry_name = new_file("pMincApply_registry"),
  registry_dir = getwd()
)
```

## Arguments

- filenames:

  Paths to the minc files to be applied accross

- fun:

  The function to apply

- ...:

  Additional arguments to fun through
  [mincApplyRCPP](https://mouse-imaging-centre.github.io/RMINC/reference/mincApplyRCPP.md)
  see notes for a warnings

- mask:

  The path to a mask for the minc files

- tinyMask:

  whether to use a small subset of voxels to test the computation

- batches:

  The number of jobs to break the computation into, ignored for
  snowfall/local mode

- slab_sizes:

  A 3 element vector indicating how large a chunk of data to read from
  each minc file at a time defaults to one slice along the first
  dimension.

- method:

  A deprecated argument formerly used to configure how to parallelize
  the jobs now this is handled with `conf_file`

- local:

  boolean whether to run the jobs locally (with the parallel package) or
  with batchtools

- cores:

  defaults to 1 or
  `max(getOption("mc.cores"), parallel::detectCores() - 1)` if running
  locally see
  [qMincApply](https://mouse-imaging-centre.github.io/RMINC/reference/qMincApply.md)
  for details.

- resources:

  A list of resources to request from the queueing system common
  examples including memory, walltime, and nodes see
  `system.file("parallel/pbs_script.tmpl", package = "RMINC")` and
  `system.file("parallel/sge_script.tmpl", package = "RMINC")` for more
  examples

- packages:

  Character vector of packages to load for all jobs

- vmem:

  The number of gigabytes of memory to request for each batched job. It
  is a compatibility argument and will overload `vmem` set in the
  resource list (if it is defined)

- walltime:

  The amount of walltime to request for each batched job. It is a
  compatibility argument and will overwrite `walltime` set in the
  resource list (if it is defined)

- workers:

  The number of workers to use. It is a compatibility option and will
  overwrite `batches` if it is supplied.

- temp_dir:

  A path to a temporary directory to hold the job registry created when
  using a true queuing system and for writing temporary mask files. This
  must be a location read/writable by all nodes when using a true
  queuing system (so /tmp will not work).

- cleanup:

  Whether to clean up registry after a queue parallelization job

- collate:

  A function to be applied to collapse the results of the the
  pMincApply. Defaults to
  [simplify2minc](https://mouse-imaging-centre.github.io/RMINC/reference/simplify2minc.md).

- conf_file:

  A batchtools config file, defaults to `option("RMINC_BATCH_CONF")`

- registry_name:

  a name for the registry

- registry_dir:

  where batchtools should create the registry

## Value

The results of applying `fun` to each voxel accross `filenames` after
collation with `collate`

## Details

This is a convenience wrapper for two underlying functions
[qMincApply](https://mouse-imaging-centre.github.io/RMINC/reference/qMincApply.md)
and
[mcMincApply](https://mouse-imaging-centre.github.io/RMINC/reference/mcMincApply.md)
for queueing and multicore processing respectively. Each of these
functions divides all of the voxels that are masked by `mask` into
`batches`. Batches are processed in parallel, with calling
[mincApplyRCPP](https://mouse-imaging-centre.github.io/RMINC/reference/mincApplyRCPP.md)
to process the voxels in the batch. Arguments passed in through `...`
will be bound by
[mincApplyRCPP](https://mouse-imaging-centre.github.io/RMINC/reference/mincApplyRCPP.md)
before `fun`, so be wary of potential partial matches. When in doubt,
partially apply your function before hand, and do not rely on positional
matching.
