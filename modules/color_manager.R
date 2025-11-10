# ==============================================================================
# COLOR MANAGEMENT MODULE
# ==============================================================================
# Centralized color management to ensure consistent colors across all plots
# (heatmaps, scatter plots, violin plots, QC plots, etc.)

#' R6 Class: ColorManager
#' 
#' Centralized color management system for consistent coloring across all plots
#' Maintains a cache of color assignments for categorical variables
#' 
#' @export
ColorManager <- R6::R6Class(
  "ColorManager",
  public = list(
    #' @field color_cache Internal cache of color mappings
    color_cache = list(),
    
    #' @field palette_type Type of palette ("Set1", "Set2", "Set3", "Paired", etc.)
    palette_type = "Set1",
    
    #' @field max_colors_per_palette Maximum colors in base palette
    max_colors_per_palette = 9,
    
    #' Initialize ColorManager
    #' @param palette_type Default palette type
    initialize = function(palette_type = "Set1") {
      self$palette_type <- palette_type
    },
    
    #' Get consistent colors for categories
    #' @param categories Vector of categories
    #' @param force_regenerate Force regeneration even if cached
    #' @return Named vector of colors
    get_colors = function(categories, force_regenerate = FALSE) {
      unique_cats <- sort(unique(as.character(categories)))
      cat_key <- paste(unique_cats, collapse = "|")
      
      # Return cached if available
      if (!force_regenerate && cat_key %in% names(self$color_cache)) {
        return(self$color_cache[[cat_key]])
      }
      
      # Generate new palette
      n_cats <- length(unique_cats)
      colors <- self$generate_palette(n_cats)
      colors <- setNames(colors, unique_cats)
      
      # Cache and return
      self$color_cache[[cat_key]] <- colors
      return(colors)
    },
    
    #' Generate color palette for n categories
    #' @param n Number of colors needed
    #' @return Vector of hex colors
    generate_palette = function(n) {
      library(RColorBrewer)
      
      if (n <= self$max_colors_per_palette) {
        brewer.pal(self$max_colors_per_palette, self$palette_type)[1:n]
      } else {
        colorRampPalette(brewer.pal(self$max_colors_per_palette, self$palette_type))(n)
      }
    },
    
    #' Map values to colors using existing color scheme
    #' @param values Vector of values to map
    #' @param colors_named Named vector of colors (from get_colors)
    #' @return Vector of color hex codes
    map_to_colors = function(values, colors_named) {
      as.character(colors_named[as.character(values)])
    },
    
    #' Get RGB values for a color (useful for transparency)
    #' @param color Hex color
    #' @param alpha Alpha value (0-1)
    #' @return rgba string for web use
    to_rgba = function(color, alpha = 1) {
      rgb_vals <- col2rgb(color) / 255
      sprintf("rgba(%d, %d, %d, %f)", 
              as.integer(rgb_vals[1] * 255),
              as.integer(rgb_vals[2] * 255),
              as.integer(rgb_vals[3] * 255),
              alpha)
    },
    
    #' Clear cache (e.g., when dataset changes)
    clear_cache = function() {
      self$color_cache <- list()
    },
    
    #' Get palette info
    #' @return List with palette info
    get_info = function() {
      list(
        palette_type = self$palette_type,
        cached_schemes = length(self$color_cache),
        n_unique_colors = length(unique(unlist(self$color_cache)))
      )
    }
  )
)

# ==============================================================================
# GLOBAL COLOR MANAGER INSTANCE
# ==============================================================================

#' Initialize global color manager
#' @param palette_type Palette type to use
#' @return ColorManager instance
init_global_color_manager <- function(palette_type = "Set1") {
  ColorManager$new(palette_type = palette_type)
}

# ==============================================================================
# COLOR UTILITIES
# ==============================================================================

#' Create consistent color mapping for annotation
#' 
#' @param annotation_data Vector of annotation values
#' @param color_manager ColorManager instance
#' 
#' @return List with colors, levels, and mappings
create_annotation_color_map <- function(annotation_data, color_manager) {
  unique_vals <- sort(unique(as.character(annotation_data)))
  colors <- color_manager$get_colors(unique_vals)
  
  list(
    colors = colors,
    levels = unique_vals,
    n_categories = length(unique_vals),
    color_for_value = function(val) colors[as.character(val)]
  )
}

#' Encode annotation data with consistent colors for JavaScript
#' 
#' @param annotation_data Vector of annotation values
#' @param color_manager ColorManager instance
#' 
#' @return List with encoded data and color info
encode_annotation_for_js = function(annotation_data, color_manager) {
  color_map <- create_annotation_color_map(annotation_data, color_manager)
  
  # Convert to numeric (0-based)
  numeric_data <- as.numeric(factor(annotation_data, levels = color_map$levels)) - 1
  
  # Encode as binary
  binary_data <- writeBin(as.numeric(numeric_data), raw(), size = 4, endian = "little")
  encoded_data <- base64enc::base64encode(binary_data)
  
  list(
    encoded = encoded_data,
    length = length(numeric_data),
    categories = color_map$levels,
    colors = color_map$colors,
    hex_colors = as.character(color_map$colors)
  )
}

#' Generate color palette for continuous data (gene expression)
#' 
#' @param n_colors Number of colors in palette
#' @param palette_name Type: "viridis", "plasma", "inferno", "cool", "warm"
#' 
#' @return Vector of hex colors
generate_continuous_palette <- function(n_colors = 256, palette_name = "viridis") {
  library(viridisLite)
  
  palette_func <- get(palette_name)
  palette_func(n_colors)
}

#' Create diverging color scale (for log-fold change, etc.)
#' 
#' @param values Numeric vector
#' @param center Center value (usually 0)
#' @param palette_name "RdBu", "PiYG", "PRGn", etc.
#' 
#' @return Named vector with colors
create_diverging_colors <- function(values, center = 0, palette_name = "RdBu") {
  library(RColorBrewer)
  
  n_colors <- 101
  colors <- brewer.pal(11, palette_name)
  color_func <- colorRampPalette(colors)(n_colors)
  
  # Normalize values to 0-1 range centered at center
  min_val <- min(values, na.rm = TRUE)
  max_val <- max(values, na.rm = TRUE)
  
  # Ensure center is in range
  range <- max(abs(min_val - center), abs(max_val - center))
  normalized <- (values - center) / range
  normalized <- (normalized + 1) / 2  # Scale to 0-1
  normalized <- pmax(pmin(normalized, 1), 0)  # Clip to 0-1
  
  color_indices <- round(normalized * (n_colors - 1)) + 1
  colors_mapped <- color_func[color_indices]
  
  names(colors_mapped) <- names(values)
  return(colors_mapped)
}

# ==============================================================================
# TESTING/VALIDATION
# ==============================================================================

#' Test color manager functionality
test_color_manager <- function() {
  # Initialize
  cm <- ColorManager$new("Set1")
  
  # Test 1: Basic color assignment
  cats <- c("A", "B", "C", "A", "B")
  colors <- cm$get_colors(cats)
  cat("✓ Test 1 - Basic colors:", paste(names(colors), collapse = ", "), "\n")
  
  # Test 2: Caching
  colors2 <- cm$get_colors(cats)
  identical(colors, colors2)
  cat("✓ Test 2 - Caching works:", identical(colors, colors2), "\n")
  
  # Test 3: Color mapping
  value_colors <- cm$map_to_colors(cats, colors)
  cat("✓ Test 3 - Value mapping works:", length(value_colors) == length(cats), "\n")
  
  # Test 4: Large palette
  large_cats <- paste0("C", 1:50)
  large_colors <- cm$get_colors(large_cats)
  cat("✓ Test 4 - Large palette (50 cats):", length(large_colors) == 50, "\n")
  
  # Test 5: Cache clearing
  cm$clear_cache()
  cat("✓ Test 5 - Cache cleared:", length(cm$color_cache) == 0, "\n")
  
  cat("\n✅ All color manager tests passed!\n")
}

# Run tests
# test_color_manager()