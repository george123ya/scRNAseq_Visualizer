library(base64enc)

encode_data <- function(processed_result) {
  df <- processed_result$data
  
  # Detect annotation columns
  ann_cols <- names(df)[grepl("_id$", names(df)) & names(df) != "dummy"]
  
  # Detect all scaled gene expression columns dynamically
  gene_cols_scaled <- names(df)[grepl("_scaled$", names(df))]
  
  # Combine required, gene, annotation, and dummy columns
  all_cols <- c("x", "y", "magic_x", "magic_y", gene_cols_scaled, ann_cols, "dummy")

  # print(paste("Encoding columns:", paste(all_cols, collapse = ", ")))

  # Encode to base64
  mat <- as.matrix(df[all_cols])
  vec <- as.numeric(t(mat))
  raw_vec <- writeBin(vec, raw(), size = 4, endian = "little")
  base64encode(raw_vec)
}
