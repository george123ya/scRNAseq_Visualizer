# ==============================================================================
# VISUALIZATION MODULE
# ==============================================================================
# Utilities for creating consistent, color-linked visualizations

#' Calculate dynamic heatmap height
#' 
#' @param n_rows Number of rows in heatmap
#' @param row_cm Height per row in cm
#' @param row_px Height per row in pixels
#' @param min_cm Minimum height in cm
#' @param max_cm Maximum height in cm
#' @param min_px Minimum height in pixels
#' @param max_px Maximum height in pixels
#' 
#' @return List with cm and px values
#' @export
calc_heatmap_height <- function(n_rows, row_cm = 0.5, row_px = 20,
                               min_cm = 6, max_cm = 40,
                               min_px = 300, max_px = 2000) {
  list(
    cm = max(min(n_rows * row_cm, max_cm), min_cm),
    px = max(min(n_rows * row_px, max_px), min_px)
  )
}

#' Calculate heatmap dimensions
#' 
#' @param n_rows Number of rows
#' @param n_cols Number of columns
#' @param row_cm Height per row in cm
#' @param col_cm Width per column in cm
#' @param min_cm Minimum dimension
#' @param max_cm Maximum dimension
#' @param min_w Minimum width
#' @param max_w Maximum width
#' 
#' @return List with height_cm and width_cm
#' @export
calc_heatmap_size <- function(n_rows, n_cols, row_cm = 0.5, col_cm = 0.3,
                             min_cm = 6, max_cm = 40, min_w = 6, max_w = 50) {
  list(
    height_cm = max(min(n_rows * row_cm, max_cm), min_cm),
    width_cm  = max(min(n_cols * col_cm, max_w), min_w)
  )
}

#' Render ComplexHeatmap with group annotations using consistent colors
#' 
#' @param data List with $mat (matrix) and $group_info (optional grouping)
#' @param title Heatmap title
#' @param color_manager ColorManager instance for consistent colors
#' @param show_row_names Show row names?
#' @param show_col_names Show column names?
#' 
#' @return ComplexHeatmap object
#' @export
render_heatmap_with_colors <- function(data, title = NULL, color_manager,
                                       show_row_names = TRUE, show_col_names = FALSE) {
  
  library(ComplexHeatmap)
  library(grid)
  
  mat <- data$mat
  n_rows <- nrow(mat)
  n_cols <- ncol(mat)
  sz <- calc_heatmap_size(n_rows, n_cols)
  
  # Create heatmap base
  ht <- if (!is.null(data$group_info)) {
    # With group annotations
    groups <- levels(data$group_info)
    colors <- color_manager$get_colors(groups)
    
    ha <- HeatmapAnnotation(
      Group = data$group_info,
      col = list(Group = colors),
      show_annotation_name = FALSE
    )
    
    Heatmap(
      mat,
      name = "expression",
      top_annotation = ha,
      cluster_columns = FALSE,
      cluster_rows = TRUE,
      show_column_names = show_col_names,
      show_row_names = show_row_names,
      column_title = title,
      use_raster = TRUE,
      raster_device = "CairoPNG",
      raster_quality = 3,
      column_split = data$group_info,
      row_names_max_width = unit(8, "cm"),
      row_names_gp = gpar(fontsize = 10)
    )
  } else {
    # Without annotations
    Heatmap(
      mat,
      name = "expression",
      cluster_rows = TRUE,
      cluster_columns = TRUE,
      show_column_names = show_col_names,
      show_row_names = show_row_names,
      column_title = title,
      use_raster = TRUE,
      raster_device = "CairoPNG",
      raster_quality = 3,
      row_names_max_width = unit(8, "cm"),
      row_names_gp = gpar(fontsize = 10)
    )
  }
  
  # Draw with proper sizing
  draw(ht, heatmap_width = unit(sz$width_cm, "cm"))
  
  return(ht)
}

# ==============================================================================
# GENE EXPRESSION PLOTTING
# ==============================================================================

#' Create gene violin/boxplot with consistent colors
#' 
#' @param gene Gene name
#' @param group_var Grouping variable (factor)
#' @param expr_data Expression values
#' @param plot_type "violin" or "boxplot"
#' @param show_points Show individual points?
#' @param log_scale Use log scale?
#' @param point_size Size of points
#' @param color_manager ColorManager instance
#' 
#' @return ggplot2 plot object
#' @export
create_gene_plot <- function(gene, group_var, expr_data, plot_type = "violin",
                            show_points = TRUE, log_scale = FALSE, 
                            point_size = 1, color_manager) {
  
  library(ggplot2)
  
  # Prepare data
  df <- data.frame(
    group = as.factor(group_var),
    expression = as.numeric(expr_data),
    stringsAsFactors = FALSE
  )
  
  # Remove NA values
  df <- na.omit(df)
  
  # Get consistent colors
  groups <- levels(df$group)
  colors <- color_manager$get_colors(groups)
  
  # Build plot
  p <- ggplot(df, aes(x = group, y = expression, fill = group)) +
    labs(
      x = "Group",
      y = "Expression Level",
      title = gene,
      fill = "Group"
    ) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "bottom"
    )
  
  # Add plot layer
  if (plot_type == "violin") {
    p <- p + geom_violin(scale = "width", trim = TRUE, alpha = 0.7)
  } else if (plot_type == "boxplot") {
    p <- p + geom_boxplot(alpha = 0.7, outlier.size = 0.5)
  }
  
  # Add points if requested
  if (show_points && nrow(df) < 10000) {
    # Use jitter for better visibility
    p <- p + geom_jitter(
      size = point_size,
      width = 0.2,
      alpha = 0.4,
      color = "black"
    )
  }
  
  # Apply log scale if requested
  if (log_scale) {
    p <- p + scale_y_continuous(trans = "log1p", labels = scales::scientific)
  }
  
  return(p)
}

#' Create multiple gene plots in grid
#' 
#' @param genes Character vector of gene names
#' @param group_var Grouping variable
#' @param expr_data Data frame with columns as genes, rows as cells
#' @param plot_type "violin" or "boxplot"
#' @param show_points Show points?
#' @param log_scale Log scale?
#' @param color_manager ColorManager instance
#' 
#' @return List of ggplot objects
#' @export
create_gene_plot_grid <- function(genes, group_var, expr_data, plot_type = "violin",
                                 show_points = TRUE, log_scale = FALSE,
                                 color_manager) {
  
  plots <- lapply(genes, function(gene) {
    if (!(gene %in% colnames(expr_data))) {
      return(ggplot() + theme_minimal() + 
             ggtitle(paste("Gene not found:", gene)))
    }
    
    create_gene_plot(
      gene = gene,
      group_var = group_var,
      expr_data = expr_data[[gene]],
      plot_type = plot_type,
      show_points = show_points,
      log_scale = log_scale,
      color_manager = color_manager
    )
  })
  
  names(plots) <- genes
  return(plots)
}

# ==============================================================================
# QC PLOT UTILITIES
# ==============================================================================

#' Create QC violin plot
#' 
#' @param data Data frame with $Metric and $Group columns
#' @param metric_name Name of metric
#' @param color_manager ColorManager instance
#' 
#' @return ggplot2 object
#' @export
create_qc_violin <- function(data, metric_name, color_manager) {
  
  library(ggplot2)
  
  groups <- levels(data$Group)
  colors <- color_manager$get_colors(groups)
  
  ggplot(data, aes(x = Group, y = Metric, fill = Group)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
    scale_fill_manual(values = colors) +
    labs(
      title = paste("QC Metric:", metric_name),
      x = "Group",
      y = metric_name,
      fill = "Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

#' Create QC scatter plot
#' 
#' @param data Data frame with $X, $Y, and optional $Group column
#' @param x_label X-axis label
#' @param y_label Y-axis label
#' @param color_manager ColorManager instance (optional)
#' 
#' @return ggplot2 object
#' @export
create_qc_scatter <- function(data, x_label, y_label, color_manager = NULL) {
  
  library(ggplot2)
  
  if (!is.null(data$Group) && !all(is.na(data$Group))) {
    groups <- unique(data$Group)
    colors <- if (!is.null(color_manager)) {
      color_manager$get_colors(groups)
    } else {
      setNames(seq_along(groups), groups)
    }
    
    plot <- ggplot(data, aes(x = X, y = Y, color = Group)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values = colors) +
      labs(
        title = paste(x_label, "vs", y_label),
        x = x_label,
        y = y_label,
        color = "Group"
      ) +
      theme_minimal() +
      theme(legend.position = "bottom")
  } else {
    plot <- ggplot(data, aes(x = X, y = Y)) +
      geom_point(alpha = 0.6, size = 2, color = "steelblue") +
      labs(
        title = paste(x_label, "vs", y_label),
        x = x_label,
        y = y_label
      ) +
      theme_minimal()
  }
  
  return(plot)
}

# ==============================================================================
# COLOR LEGEND UTILITIES
# ==============================================================================

#' Create color legend as standalone plot
#' 
#' @param categories Character vector of categories
#' @param colors Named character vector of colors
#' @param title Legend title
#' 
#' @return ggplot2 object
#' @export
create_color_legend <- function(categories, colors, title = "Categories") {
  
  library(ggplot2)
  
  legend_data <- data.frame(
    category = factor(categories, levels = categories),
    x = 1,
    y = seq_along(categories)
  )
  
  ggplot(legend_data, aes(x = x, y = y, fill = category)) +
    geom_tile(width = 0.5, height = 0.8) +
    scale_fill_manual(values = colors, guide = "none") +
    geom_text(aes(label = category), x = 1.5, hjust = 0) +
    xlim(0, 3) +
    ylim(0, length(categories) + 1) +
    labs(title = title) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

#' Export color mapping for JavaScript/front-end
#' 
#' @param categories Character vector of categories
#' @param colors Named character vector of colors
#' 
#' @return JSON string
#' @export
export_colors_as_json <- function(categories, colors) {
  
  color_map <- setNames(as.list(colors), names(colors))
  jsonlite::toJSON(color_map, pretty = TRUE)
}

# ==============================================================================
# TESTING
# ==============================================================================

#' Test visualization module
#' 
#' @export
test_visualization_module <- function() {
  
  cat("Testing visualization utilities...\n")
  
  # Test 1: Heatmap height calculation
  h <- calc_heatmap_height(100)
  cat("✓ Heatmap height (100 rows):", h$px, "px\n")
  
  # Test 2: Gene plot creation
  test_data <- data.frame(
    group = rep(c("A", "B", "C"), times = 33),
    expr = rnorm(99)
  )
  
  cm <- ColorManager$new()
  p <- create_gene_plot(
    gene = "TEST_GENE",
    group_var = test_data$group,
    expr_data = test_data$expr,
    color_manager = cm
  )
  
  cat("✓ Gene plot created:", class(p)[1], "\n")
  
  # Test 3: Color legend
  legend <- create_color_legend(
    categories = c("A", "B", "C"),
    colors = cm$get_colors(c("A", "B", "C")),
    title = "Test Categories"
  )
  
  cat("✓ Color legend created:", class(legend)[1], "\n")
  
  cat("\n✅ All visualization tests passed!\n")
}

# test_visualization_module()