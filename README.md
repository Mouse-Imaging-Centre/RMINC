# RMINC

[![R Tests](https://github.com/Mouse-Imaging-Centre/RMINC/actions/workflows/R-tests.yml/badge.svg?branch=develop)](https://github.com/Mouse-Imaging-Centre/RMINC/actions/workflows/R-tests.yml)
[![License: BSD 3-Clause](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://github.com/Mouse-Imaging-Centre/RMINC/blob/develop/LICENSE)

Statistical Tools for Medical Imaging NetCDF (MINC) Files.

**Unix only** (Linux, macOS). Requires **R >= 4.0**.

## What is MINC?

MINC (Medical Imaging NetCDF) is a medical imaging format developed at the
McConnell Brain Imaging Centre (BIC) at McGill University. It is built on top of
HDF5 and designed for efficient storage and processing of multi-dimensional
neuroimaging data. RMINC provides a native R interface for working with these
files.

## Features

- **Voxel-wise statistics** — linear models (`mincLm`), mixed-effects models
  (`mincLmer`), ANOVA (`mincAnova`), t-tests, Wilcoxon tests, and flexible
  models with voxel-varying covariates (`flexLm`)
- **Surface/vertex statistics** — vertex-wise linear and mixed-effects models,
  vertex TFCE, and vertex-level FDR correction
- **ROI anatomy analysis** — atlas-based volume extraction (`anatGetAll`),
  structure-wise statistics (`anatLm`, `anatLmer`), and hierarchical atlas
  summarization
- **Multiple comparisons correction** — FDR (Benjamini-Hochberg andStorey's
  q-value) and Threshold-Free Cluster Enhancement (TFCE)
- **Effect sizes** — Hedges' G for lm/lmer models
- **Model comparison** — AIC, BIC, and AICc with log-likelihood reporting
- **Visualization** — 2D slice series plots, interactive 3D surface rendering
  via `rgl`, colour lookup tables
- **Interactive exploration** — Shiny app (`launch_shinyRMINC`) for browsing
  anatomy and statistics
- **Parallel execution** — distributed computing via `batchtools` for local
  multicore and HPC cluster backends
- **CIVET pipeline support** — reading CIVET 1.1.12/2.0/2.1 outputs, QC
  data, and brain-view ROI generation
- **Large data support** — file-backed matrices via `mincTable` for datasets
  exceeding available memory
- **Peak finding** — local optima detection in MINC volumes

## Installation

The easiest way to install is with devtools:

```r
devtools::install_github("Mouse-Imaging-Centre/RMINC@develop"
                       , upgrade_dependencies = FALSE)
```

### System Requirements

RMINC depends on the MINC C library. You have two options:

1. **minc-toolkit-v2** (recommended): Download binary installers from
   http://bic-mni.github.io/ or build from source at
   https://github.com/BIC-MNI/minc-toolkit-v2. This provides the MINC
   library plus command-line tools.

2. **libminc from source**: If the toolkit is not found, RMINC's `configure`
   will attempt to build libminc for you. This requires CMake (> 2.6), git,
   and HDF5 development headers.

On Debian/Ubuntu:

```
sudo apt-get install libhdf5-dev cmake git autoconf automake libtool
```

On macOS (Homebrew):

```
brew install hdf5 cmake autoconf automake libtool
```

See the [INSTALL](https://github.com/Mouse-Imaging-Centre/RMINC/blob/develop/INSTALL)
file for detailed instructions, including macOS-specific setup (gfortran,
libz linking).

If `configure` cannot find MINC, set the install prefix:

```bash
export MINC_TOOLKIT=/opt/minc/1.9.18
```

Or pass the path at install time from R:

```r
devtools::install_github("Mouse-Imaging-Centre/RMINC@develop"
  , args = "--configure-args='--with-build-path=/opt/minc/1.9.18'"
  , upgrade_dependencies = FALSE)
```

## Quick Example

```r
library(RMINC)

# Define a data frame mapping MINC files to covariates
gf <- data.frame(
  file = c("subject1.mnc", "subject2.mnc", "subject3.mnc", "subject4.mnc"),
  group = factor(c("control", "treatment", "control", "treatment"))
)

# Run a voxel-wise linear model testing group differences
result <- mincLm(file ~ group, data = gf, mask = "mask.mnc")

# Apply FDR correction
qvals <- mincFDR(result)

# Visualize a slice series
mincPlotSliceSeries(result, column = "groupcontrol", threshold = 2.5)
```

## Dockerized Testing

Build the test image:

```
docker build -t rminc-test .
```

Run checks and tests (mount the source directory):

```
docker run -v $(pwd):/RMINC rminc-test              # check + tests (default)
docker run -v $(pwd):/RMINC rminc-test check          # R CMD check --as-cran only
docker run -v $(pwd):/RMINC rminc-test test            # tests only
docker run -v $(pwd):/RMINC rminc-test test mincLm     # filtered tests
docker run -v $(pwd):/RMINC rminc-test shell           # interactive bash shell
```

## Tutorials

- [Voxel-wise Statistics](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/VBMstats.pdf)
- [How RMINC Parallel Works](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/RMINC_Parallelism.html)
- [RMINC 2D Visualization](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/visualizationTutorial.html)
- [Visualizing 3D Objects with RMINC](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/RMINC_rgl.html)
- [Analyzing Anatomy with Hierarchical Atlases](https://raw.githubusercontent.com/Mouse-Imaging-Centre/RMINC/develop/inst/documentation/hierarchiesTutorial.html)

## Getting Help

Function documentation is available through R's `?` interface (e.g. `?mincLm`).
Report bugs or request features at the
[issues page](https://github.com/Mouse-Imaging-Centre/RMINC/issues).

## Contributing

Active development happens on the `develop` branch. The `master` branch may lag
behind. Contributions via pull requests to `develop` are welcome.

## Citation

To cite RMINC in publications, use the citation information returned by:

```r
citation("RMINC")
```

## License

RMINC is released under the [BSD 3-Clause License](LICENSE).

Copyright (c) 2003, Jason Lerch, The Hospital for Sick Children.

## Authors

- **Jason Lerch** — author
- **Chris Hammill** — author, maintainer
- **Matthijs van Eede** — author
- **Daniel Cassel** — author
- Jim Nikelski, Yohan Yee, Miriam Friedel, Nick Wang, Gabriel Devenyi,
  Benjamin Darwin, Bilal Syed — contributors
