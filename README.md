# RMINC

### Statistical Tools for Medical Imaging NetCDF (MINC) Files

[![R Tests](https://github.com/Mouse-Imaging-Centre/RMINC/actions/workflows/R-tests.yml/badge.svg?branch=develop)](https://github.com/Mouse-Imaging-Centre/RMINC/actions/workflows/R-tests.yml)
[![License: BSD-3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-blue.svg)](https://github.com/Mouse-Imaging-Centre/RMINC/blob/develop/LICENSE)
[![R >= 4.0](https://img.shields.io/badge/R-%3E%3D%204.0-blue.svg)](https://www.r-project.org/)
[![Unix](https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey.svg)](#installation)

Tools for reading, writing, analyzing, and visualizing MINC neuroimaging files
with R. Supports voxel-wise linear models, mixed-effects models, correlations,
ANOVA, FDR correction, TFCE, ROI-based anatomy analysis, surface-based vertex
statistics, and parallel execution.

> **Unix only** (Linux, macOS). Requires **R >= 4.0**. No native Windows support
> (WSL may work — see [INSTALL](https://github.com/Mouse-Imaging-Centre/RMINC/blob/develop/INSTALL)).

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
  - [System Requirements](#system-requirements)
  - [Troubleshooting](#troubleshooting)
- [Quick Example](#quick-example)
- [Dockerized Testing](#dockerized-testing)
- [Tutorials](#tutorials)
- [Getting Help](#getting-help)
- [Development](#development)
- [License](#license)

## Features

| Capability | Key functions |
| --- | --- |
| **Voxel-wise linear models** | `mincLm`, `mincAnova`, `mincTtest`, `mincWilcoxon` |
| **Correlations & summaries** | `mincCorrelation`, `mincMean`, `mincSummary` |
| **Mixed-effects models** | `mincLmer`, `mincLmerEstimateDF`, `mincLogLikRatio` |
| **Surface / vertex statistics** | `vertexLm`, `vertexLmer`, `vertexTFCE` |
| **ROI / anatomy analysis** | `anatLm`, `anatLmer`, `hanatLm` (hierarchical atlases) |
| **Multiple comparisons** | `mincFDR`, `mincTFCE`, `mincRandomize` |
| **I/O** | `mincGetVolume`, `mincWriteVolume`, `mincArray` |
| **Visualization** | `mincPlotSliceSeries`, `mincTriplanarSlicePlot` (2D), `rgl` (3D), `hanatView` |
| **Parallel execution** | `pMincApply`, `qMincApply`, `mincLm(parallel = ...)` |

## Installation

The easiest way to install is with devtools from the `develop` branch:

```r
devtools::install_github(
  "Mouse-Imaging-Centre/RMINC@develop",
  upgrade_dependencies = FALSE
)
```

### System Requirements

RMINC depends on the MINC C library. You have two options:

1. **minc-toolkit-v2** (recommended): Download binary installers from
   <http://bic-mni.github.io/> or build from source at
   <https://github.com/BIC-MNI/minc-toolkit-v2>. This provides the MINC library
   plus command-line tools.

2. **libminc from source**: If the toolkit is not found, RMINC's `configure`
   will attempt to build libminc for you. This requires CMake (> 2.6), git, and
   the HDF5 development headers.

**Debian / Ubuntu:**

```bash
sudo apt-get install libhdf5-dev cmake git lsof autoconf automake libtool
```

**macOS (Homebrew):**

```bash
brew install hdf5 cmake autoconf automake libtool
```

> **macOS note:** RMINC's compiled backend needs a Fortran compiler. Apple's
> default clang does not include one — install `gcc` via Homebrew
> (`brew install gcc libssh2 libgit2 libomp zlib`) or
> [gfortran](https://gcc.gnu.org/wiki/GFortranBinaries) directly, and point R at
> it via `~/.R/Makevars`. See the
> [INSTALL](https://github.com/Mouse-Imaging-Centre/RMINC/blob/develop/INSTALL)
> file for the exact `FLIBS`/`LDFLAGS` setup and detailed instructions.

### Troubleshooting

If `configure` cannot find MINC, set the install prefix via an environment
variable:

```bash
export MINC_TOOLKIT=/opt/minc/1.9.18
```

Or pass the path at install time from R:

```r
devtools::install_github(
  "Mouse-Imaging-Centre/RMINC@develop",
  args = "--configure-args='--with-build-path=/opt/minc/1.9.18'",
  upgrade_dependencies = FALSE
)
```

## Quick Example

Fit a voxel-wise linear model across a set of MINC volumes. The response column
in `data` holds MINC file paths; the model is fit independently at every voxel:

```r
library(RMINC)

# `gf` is a data frame whose `volumes` column holds MINC file paths
fit <- mincLm(volumes ~ group, data = gf)

# Output columns: "F-statistic", "R-squared", beta coefficients,
# t-statistics, "logLik". Write the whole-model F-statistic (column 1):
mincWriteVolume(fit, output.filename = "group_Fstat.mnc", column = 1)
```

Function-level documentation is available in R via `?mincLm`, `?mincLmer`,
`?vertexLm`, `?anatLm`, etc.

## Dockerized Testing

Build the test image:

```bash
docker build -t rminc-test .
```

Run checks and tests (mount the source directory):

```bash
docker run -v $(pwd):/RMINC rminc-test              # check + tests (default)
docker run -v $(pwd):/RMINC rminc-test check          # R CMD check --as-cran only
docker run -v $(pwd):/RMINC rminc-test test            # tests only
docker run -v $(pwd):/RMINC rminc-test test mincLm     # filtered tests
docker run -v $(pwd):/RMINC rminc-test shell           # interactive bash shell
```

## Tutorials

- [Voxel-wise Statistics](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/VBMstats.pdf) (somewhat out of date)
- [How RMINC Parallel Works](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/RMINC_Parallelism.html)
- [RMINC 2D Visualization](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/visualizationTutorial.html)
- [Visualizing 3D Objects with RMINC](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/RMINC_rgl.html)
- [Analyzing Anatomy with Hierarchical Atlases](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/hierarchiesTutorial.html)

## Getting Help

Function documentation is available through R's `?` interface (e.g. `?mincLm`).
The [MICePUB wiki](https://wiki.mouseimaging.ca/display/MICePub/RMINC) has
additional background and usage notes. Report bugs or request features at the
[issues page](https://github.com/Mouse-Imaging-Centre/RMINC/issues).

## Development

The `develop` branch is where active development happens; the `master` branch
may lag behind. Contributions via pull requests to `develop` are welcome.

The C/C++ backend is built with Rcpp and autoconf — `src/Makevars` is generated
from `src/Makevars.in`, and `RcppExports.{R,cpp}` are auto-generated by Rcpp.
Do not edit these files directly.

After editing R source containing roxygen comments, regenerate documentation
and the `NAMESPACE` with:

```r
roxygen2::roxygenise()
```

## License

RMINC is released under the [BSD 3-Clause License](LICENSE).
