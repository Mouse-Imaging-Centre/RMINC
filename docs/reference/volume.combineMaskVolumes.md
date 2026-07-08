# Combine Multiple Mask Volumes into a Single Mask Volume

Given a list containing more than one label volume, combine those
volumes, creating an aggregate mask volume.

## Usage

``` r
volume.combineMaskVolumes(vol_list)
```

## Arguments

- vol_list:

  A list containing more than one mask volume. Note that all volumes
  must reflect the same sampling.

## Value

A single aggregate MincVolumeIO volume is returned.

## Author

Jim Nikelski <nikelski@bic.mni.mcgill.ca>
