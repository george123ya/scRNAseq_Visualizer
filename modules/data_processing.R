source("modules/utils.R")

detect_annotations <- function(df) {
  required <- c("x", "y", "magic_x", "magic_y")
  gene_patterns <- c("gene", "expr", "expression", "count", "value")
  all_cols <- names(df)
  
  # Detect multiple gene expression columns
  gene_cols <- unlist(lapply(gene_patterns, function(pattern) {
    grep(pattern, all_cols, ignore.case = TRUE, value = TRUE)
  }))
  gene_cols <- unique(gene_cols)

  # If none detected, pick a numeric column as fallback
  if (length(gene_cols) == 0) {
    numeric_cols <- names(Filter(is.numeric, df))[!names(Filter(is.numeric, df)) %in% required]
    if (length(numeric_cols) > 0) gene_cols <- numeric_cols[1]
  }
  
  # Annotation columns (categorical or low-unique numeric)
  annotation_cols <- setdiff(all_cols, c(required, gene_cols))

  # annotation_cols <- annotation_cols[sapply(df[annotation_cols], function(x)
  #   is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) <= 20))]

  # print("awo")

  list(
    gene_cols = gene_cols,
    annotation_cols = annotation_cols,
    required_cols = required
  )
}

process_dataframe <- function(df) {
  col_info <- detect_annotations(df)
  required <- col_info$required_cols

  # Check for required columns
  if (length(setdiff(required, names(df))) > 0)
    stop("Missing required columns: ", paste(setdiff(required, names(df)), collapse = ", "))

  processed_df <- data.frame(
    x = df$x, y = df$y,
    magic_x = df$magic_x, magic_y = df$magic_y,
    stringsAsFactors = FALSE
  )

  gene_expr_ranges <- list()
  
  # Scale each detected gene expression column
  if (length(col_info$gene_cols) > 0) {
    for (gene in col_info$gene_cols) {
      raw_gene <- as.numeric(df[[gene]])
      scaled_gene <- if (max(raw_gene, na.rm = TRUE) > min(raw_gene, na.rm = TRUE)) {
        (raw_gene - min(raw_gene, na.rm = TRUE)) / (max(raw_gene, na.rm = TRUE) - min(raw_gene, na.rm = TRUE))
      } else rep(0.5, length(raw_gene))
      
      processed_df[[paste0(gene, "_scaled")]] <- scaled_gene
      gene_expr_ranges[[gene]] <- range(raw_gene, na.rm = TRUE)
    }
  } #else {
  #  processed_df$gene_expr_scaled <- rep(0.5, nrow(df))
  #}

  # Process annotation columns
  annotation_info <- list()
  for (col_name in col_info$annotation_cols) {
    values <- as.factor(if (is.numeric(df[[col_name]])) as.factor(df[[col_name]]) else df[[col_name]] )
    color_map <- generate_colors(length(levels(values)))

    annotation_info[[col_name]] <- list(
      names = levels(values),
      colors = color_map,
      values = as.integer(values)
    )
    processed_df[[paste0(col_name, "_id")]] <- as.integer(values)
  }


  # Metacell colors
  metacell_colors <- if ("cluster" %in% names(annotation_info)) {
    unname(sapply(annotation_info$cluster$colors, darken_and_saturate_color))
  } else character(0)

  # print(metacell_colors)

  processed_df$dummy <- rep(1, nrow(df))

  # print(annotation_info)


  list(
    data = processed_df,
    annotations = annotation_info,
    gene_cols = col_info$gene_cols,
    n_annotations = length(annotation_info),
    cluster_colors = annotation_info$cluster$colors,
    metacell_colors = metacell_colors,
    gene_expr_ranges = gene_expr_ranges
  )
}

