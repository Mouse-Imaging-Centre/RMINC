# This file serves three purposes - a quick and dirty way to document
# what RMINC is capable of, a way of prototyping what RMINC should be
# capable of in the future, and simultaneously runnable code to make
# sure that it actually is capable of what it claims to be capable of.

version <- 0.4

# import the data
gf <- read.csv("sample-in.csv")

# everything below here should work for version 0.4


####
# Linear Models
####

# basic linear model
# mincLm - takes a formula as an argument, the first part of which specifies
#          the filenames, the other sides the covariates.
vs <- mincLm(jacobians ~ genotype, gf)

# just typing the output should print some information about the model
vs
# example output here

# mincLm with subsetting - only parts of the data will be used
vs <- mincLm(jacobians ~ genotype, subset=coil==1, gf)

# mincLm with masking - linear model only computed inside the mask
vs <- mincLm(jacobians ~ genotype, mask="mask.mnc", gf)


# missing values: by default missing values should fail
gf$genotype[3] <- NA
vs <- mincLm(jacobians ~ genotype, gf)

# but with an option they should be removed
vs <- mincLm(jacobians ~ genotype, gf, na.rm=TRUE)

# t-tests should be just as simple:
t <- mincTtest(jacobians ~ genotype, gf)
# same deal but use two columns of filenames rather than a formula
t <- mincTtest(jacobians[genotype=="+"], jacobians[genotype=="-"])
# masking should work in the same way
t <- mincTtest(jacobians ~ genotype, mask="mask.mnc", gf)
# and subsetting ought to work in the same way as for mincLm
t <- mincTtest(jacobians ~ genotype, gf, subset=coil==1)

# correlations need the filenames and the correlate.
r <- mincCor(jacobians, age, gf)
# subsetting is the same as above.
r <- mincCor(jacobians, age, subset=coil==1, gf)

# everything below here is to be implemented for version 0.5
if (version >= 0.5) {

  # paired t-tests
  paired <- mincTtest(jacobians[genotype=="+"], jacobians[genotype=="-"],
                      paired=TRUE, gf)
  
  # mask should be either a filename or a vector of equal length to the
  # input files
  mask <- mincGetVolume("mask.mnc")
  vs <- mincLm(jacobians ~ genotype, mask=mask, gf)
  
  # linear mixed effect models
  vs <- mincLme(jacobians ~ genotype, random=~1|coil, gf)
  vs
}

# everything below here is for the distant version 0.6
if (version >= 0.6) {

  # hotelling T2
  vs <- mincHt2(displacements ~ genotype, gf)
  vs
}
