library(colorspace)
library(viridisLite)

viridis_hex <- substr(viridisLite::viridis(256), 1, 7)

generate_colors <- function(n) {
  # Paul Tol's qualitative palette (colorblind-friendly, good on white bg)
  tol_colors <- c(
    "#332288", # dark blue
    "#88CCEE", # cyan
    "#44AA99", # teal
    "#117733", # green
    "#999933", # olive
    "#DDCC77", # sand
    "#CC6677", # rose
    "#882255", # wine
    "#AA4499", # purple
    "#DDDDDD"  # light gray
  )
  
  if (n <= length(tol_colors)) {
    return(tol_colors[1:n])
  } else {
    # For larger sets, interpolate using HCL colors (colorblind-friendly)
    # return(grDevices::hcl.colors(n, palette = "Dynamic"))
    # Use highly distinct hues via HCL (even spacing in hue)
    hues <- seq(15, 375, length.out = n + 1)[-1]
    return(hcl(h = hues, c = 80, l = 60))
  }
}

darken_and_saturate_color <- function(hex, darken = 0.35, saturate = 2) {
  col <- hex2RGB(hex)
  hsl <- as(col, "HLS")
  h <- hsl@coords[1,1]
  l <- hsl@coords[1,2]
  s <- hsl@coords[1,3]
  hsl_new <- HLS(H = h, L = pmax(pmin(l * darken, 1), 0), S = pmax(pmin(s * saturate, 1), 0))
  rgb_new <- as(hsl_new, "RGB")
  hex(rgb_new)
}
