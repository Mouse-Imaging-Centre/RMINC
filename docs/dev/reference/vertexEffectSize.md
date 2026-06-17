# Effect Sizes

Takes the output of a minc modelling function and computes the unbiased
hedges g\* and variance of hedges g\*

## Usage

``` r
vertexEffectSize(buffer, predictors = NULL)

mincEffectSize(buffer, predictors = NULL)

anatEffectSize(buffer, predictors = NULL)
```

## Arguments

- buffer:

  The results of a vertex/anat/mincLm run

- predictors:

  A vector of factor predictor names. By default the effect size be
  computed for all treatment-coded factor columns.

## Value

An object with columns of hedgesg-\<factorlevel\> and
hedgesg_var-\<factorlevel\> for each factor predictor in the GLM. The
class and attributes of the input are preserved:

- `mincEffectSize`/`mincMultiDim` for voxel-wise inputs

- `vertexEffectSize`/`vertexMultiDim` for vertex-wise inputs

- `anatEffectSize`/`anatModel` for anatomy-wise inputs

Attributes `likeVolume`, `filenames`, `model`, `data`, `call`, `df`,
`atlas`, `definitions`, and `stat-type` are carried over from the input.

## Details

This code implements the methods from Nakagawa, S., Cuthill, I.C., 2007.
Effect size, confidence interval and statistical significance: a
practical guide for biologists. Biol. Rev. Camb. Philos. Soc. 82,
591-05. https://doi.org/10.1111/j.1469-185X.2007.00027.x for computing
effect size of group comparisons from a GLM.

For now, interactions are explicitly excluded from being predictors. To
get effect size for interactions, use the interaction() function to
create a new treatment coded factor to use as a predictor.

## Functions

- `mincEffectSize()`: mincEffectSize

- `anatEffectSize()`: anatEffectSize

## Examples

``` r
if (FALSE) { # \dontrun{
getRMINCTestData()
# read the text file describing the dataset
gf <- read.csv("/tmp/rminctestdata/test_data_set.csv")
# run a linear model relating the data in all voxels to Genotype
vs <- mincLm(jacobians_fixed_2 ~ Sex, gf)
effectsize <- mincEffectSize(vs)
} # }
```
