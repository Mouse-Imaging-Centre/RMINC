# Voxel Atlas Definitions

Voxel based analyses often require an atlas to differentiate voxels by
structure membership. In conjunction with the atlas, a set of atlas
definitions allow you to calculate summaries by region for measures and
statistics of interest. For example the `defs` argument of
[anatGetAll](https://mouse-imaging-centre.github.io/RMINC/reference/anatGetAll.md)
is a character vector pointing to atlas definitions.

## Details

the basic format for voxel atlas definitions is a 3-column csv file. The
file should be formatted like:\

|                                   |             |            |
|-----------------------------------|-------------|------------|
| Structure                         | right.label | left.label |
| amygdala                          | 51          | 151        |
| anterior commisure: pars anterior | 115         | 215        |

The csv file may contain additional columns, but the three columns above
must by present with names spelling/case sensitive. Additional columns
may include a tissue.type column which is used in some functions.
