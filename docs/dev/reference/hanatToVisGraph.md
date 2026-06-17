# Turns an anatomical tree into a nodes and edges for visNetwork plotting

Turns an anatomical tree into a nodes and edges for visNetwork plotting

## Usage

``` r
hanatToVisGraph(
  hanatTree,
  colourVariable = "color_hex_triplet",
  colourScale = colorRampPalette(c("red", "yellow"))(255),
  rColourScale = colorRampPalette(c("blue", "turquoise1"))(255),
  low = NULL,
  high = NULL,
  symmetric = F,
  transparent = "#FDFDFD",
  edgeColourFromABI = F
)

hanatView(..., fontsize = 14, levelSeparation = 500)
```

## Arguments

- hanatTree:

  The input anatomical tree

- colourVariable:

  String containing the variable to colour nodes with

- colourScale:

  Colour scale to use

- rColourScale:

  Reverse colour scale to use if symmetric colours enabled

- low:

  Lower cut-off for colour scale

- high:

  Upper cut-off for colour scale

- symmetric:

  Boolean for whether colour scale is symmetric

- transparent:

  Colour to use for transparent nodes

- edgeColourFromABI:

  Whether to set edge colours based on ABI atlas

- ...:

  Extra arguments to `hanatToVisGraph` from `hanatView`

- fontsize:

  The font size for `hanatView`

- levelSeparation:

  Spacing between hierarchical layers for `hanatView`

## Value

list with nodes and edges data frames

## Functions

- `hanatView()`: Directly view the graph
