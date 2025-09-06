## ---------------------------------------------------------------
## Dependency checks and safe imports
## ---------------------------------------------------------------

# List of required R packages
required_pkgs <- c(
  "reticulate",
  "ComplexHeatmap",
  "RColorBrewer",
  "grid",
  "peakRAM"
)

# Check and load R packages
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("❌ Required R package '", pkg, "' is not installed. ",
         "Please install it with install.packages('", pkg, "')")
  } else {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }
}

# Python dependencies needed for fast Zarr/H5AD loading
py_required <- c("zarr", "fsspec", "pandas", "numpy", "requests", "aiohttp")

# Install Python deps if missing
for (pkg in py_required) {
  if (!reticulate::py_module_available(pkg)) {
    message("⚠️ Python module '", pkg, "' not found. Installing...")
    reticulate::py_install(pkg, pip = TRUE)
  }
}

# Safe imports with delay_load = TRUE
zarr   <- reticulate::import("zarr", delay_load = TRUE)
fsspec <- reticulate::import("fsspec", delay_load = TRUE)
pd     <- reticulate::import("pandas", delay_load = TRUE)
np     <- reticulate::import("numpy", delay_load = TRUE)


init_fast_zarr_h5ad <- function(file_path) {
  cat("Initializing fast loader:", file_path, "\n")
  
  # Detect file type
  is_zarr_file <- grepl("\\.zarr/?$", file_path, ignore.case = TRUE)
  is_h5ad_file <- grepl("\\.h5ad$", file_path, ignore.case = TRUE)
  
  if (is_zarr_file) {
    cat("📁 Detected Zarr format - using fast Zarr loader\n")
    return(load_zarr_fast(file_path))
  } else if (is_h5ad_file) {
    cat("📁 Detected H5AD format - using enhanced lazy loader\n")
    return(init_lazy_h5ad_enhanced(file_path))
  } else {
    # Try to auto-detect or default to lazy loader
    cat("❓ Unknown format - attempting lazy loader\n")
    return(init_lazy_h5ad_enhanced(file_path))
  }
}

load_zarr_fast <- function(file_path) {
  start_time <- Sys.time()
  
  # Import required Python modules
  zarr <- import("zarr")
  fsspec <- import("fsspec")
  pd <- import("pandas")
  np <- import("numpy")
  
  tryCatch({
    # Open Zarr store based on path type
    if (grepl("^https?://", file_path)) {
      cat("🌐 Loading from URL...\n")
      store <- fsspec$get_mapper(file_path)
      z <- zarr$open(store, mode = "r")
    } else {
      cat("💾 Loading from local path...\n")
      z <- zarr$open(file_path, mode = "r")
    }
    
    end_time <- Sys.time()
    cat("⏱ Zarr store opening took:", round(difftime(end_time, start_time, units = "secs"), 3), "seconds\n")
    
    return(process_zarr_data_fast(z, file_path))
    
  }, error = function(e) {
    cat("❌ Fast Zarr loading failed:", e$message, "\n")
    cat("🔄 Falling back to lazy AnnData loader...\n")
    return(init_lazy_h5ad_enhanced(file_path))
  })
}

process_zarr_data_fast <- function(z, file_path) {
  pd <- import("pandas")
  np <- import("numpy")
  
  # Explore Zarr structure first
  cat("🔍 Exploring Zarr structure...\n")
  zarr_keys <- py_to_r(reticulate::import_builtins()$list(z$keys()))
  cat("📋 Top-level groups:", paste(zarr_keys, collapse = ", "), "\n")
  
  # Handle different Zarr structures
  if ("obs" %in% zarr_keys && "var" %in% zarr_keys) {
    obs_group <- z[["obs"]]
    var_group <- z[["var"]]
    X_array <- z[["X"]]
  } else {
    cat("⚠️ Non-standard Zarr structure detected\n")
    stop("❌ Could not find standard obs/var/X structure in Zarr file")
  }
  
  # Get basic dimensions
  if (reticulate::py_has_attr(X_array, "shape")) {
    X_shape <- X_array$shape
    n_cells <- X_shape[[1]]
    n_genes <- X_shape[[2]]
    cat("📊 Dataset shape:", n_cells, "cells,", n_genes, "genes\n")
  } else {
    cat("⚠️ Could not determine dataset shape from X array\n")
    n_cells <- "unknown"
    n_genes <- "unknown"
  }
  
  # Fast UMAP coordinate loading
  umap_start <- Sys.time()
  umap_coords <- load_umap_fast(z)
  umap_end <- Sys.time()
  cat("⏱ Fast UMAP loading took:", round(difftime(umap_end, umap_start, units = "secs"), 3), "seconds\n")
  
  # Load gene information
  gene_start <- Sys.time()
  gene_info <- load_gene_info_fast(var_group)
  gene_end <- Sys.time()
  cat("⏱ Gene info loading took:", round(difftime(gene_end, gene_start, units = "secs"), 3), "seconds\n")
  
  # Get available obs columns
  obs_columns <- get_obs_columns_fast(obs_group)
  cat("🔑 Available obs columns:", paste(head(obs_columns, 10), collapse = ", "), 
      if(length(obs_columns) > 10) "..." else "", "\n")
  
  # Get available obsm and layers
  var_keys <- get_var_keys_fast(z, zarr_keys)
  obsm_keys <- get_obsm_keys_fast(z, zarr_keys)
  layer_keys <- get_layer_keys_fast(z, zarr_keys)
  
  cat("🗝️ Available obsm keys:", paste(obsm_keys, collapse = ", "), "\n")
  cat("📋 Available layers:", paste(layer_keys, collapse = ", "), "\n")
  
  # Load SEACells if available
  seacell_df <- load_seacells_fast(z, zarr_keys)
  
  # Load pathway data
  # if (!requireNamespace("msigdbr", quietly = TRUE)) {
  #   cat("⚠️ Installing msigdbr for pathway data...\n")
  #   install.packages("msigdbr")
  # }
  # library(msigdbr)
  
  # pathways_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
  # pathways <- split(pathways_df$gene_symbol, pathways_df$gs_name)
  # cat("🛤️ Loaded", length(pathways), "Reactome pathways\n")
  
  # # Filter pathways
  # valid_genes <- toupper(gene_info$display_names)
  # pathways <- lapply(pathways, function(genes) intersect(toupper(genes), valid_genes))
  # pathways <- pathways[sapply(pathways, length) > 0]
  # cat("🛤️ Retained", length(pathways), "pathways with valid genes\n")
  
  # Create shortened pathway names
  # pathway_names <- names(pathways)
  # pathway_display_names <- sub("^REACTOME_", "", pathway_names)  # Remove REACTOME_ prefix
  # pathway_display_names <- sapply(pathway_display_names, function(name) {
  #   if (nchar(name) > 30) paste0(substr(name, 1, 27), "...") else name
  # })
  # names(pathways) <- pathway_names  # Keep full names as keys
  # pathway_map <- data.frame(
  #   full_name = pathway_names,
  #   display_name = pathway_display_names,
  #   stringsAsFactors = FALSE
  # )
  # cat("🛤️ Sample display names:", paste(head(pathway_display_names, 5), collapse = ", "), "\n")
  
  # Create fast access functions
  get_obs_column_fast <- function(z, name) {
    to_r <- reticulate::py_to_r
    builtins <- reticulate::import_builtins()
    
    if (!inherits(z, "python.builtin.object")) stop("⚠️ z is not a Python object.")
    
    root_keys <- to_r(builtins$list(z$keys()))
    if (!("obs" %in% root_keys)) return(NULL)
    
    obs_keys <- to_r(builtins$list(z[["obs"]]$keys()))
    if (!(name %in% obs_keys)) return(NULL)
    
    col_obj <- z[["obs"]][[name]]
    
    tryCatch({
      zarr <- import("zarr")
      py$col_obj <- col_obj
      py$zarr <- zarr
      
      is_group <- py_eval("isinstance(col_obj, zarr.Group)")
      
      if (is_group) {
        codes <- to_r(col_obj[['codes']][])
        categories <- to_r(col_obj[['categories']][])
        return(categories[codes + 1])
      } else {
        return(to_r(col_obj[]))
      }
    }, error = function(e) {
      cat("⚠️ Error loading obs column", name, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  get_gene_expression_fast <- function(gene_name, layer = "X") {
    gene_idx <- match(gene_name, gene_info$display_names)
    cat("🔍 Looking for gene:", gene_name, "Found at index:", gene_idx, "\n")
    if (is.na(gene_idx)) {
      cat("❌ Gene", gene_name, "not found\n")
      return(NULL)
    }
    
    tryCatch({
      if (layer == "X") {
        peak <- peakRAM({
          expr_vec <- reticulate::py$get_gene_fast(z, as.integer(gene_idx - 1), 'X')
        })
        print(peak)
        return(py_to_r(expr_vec))
      } else {
        if (!(layer %in% layer_keys)) {
          cat("❌ Layer", layer, "not available\n")
          return(NULL)
        }
        expr_vec <- reticulate::py$get_gene_fast(z, as.integer(gene_idx - 1), layer)
        return(py_to_r(expr_vec))
      }
    }, error = function(e) {
      cat("⚠️ Error loading expression for", gene_name, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  get_genes_expression_batch <- function(gene_names, layer = "X") {
    gene_indices <- match(gene_names, gene_info$display_names)
    valid_mask <- !is.na(gene_indices)
    valid_indices <- as.integer(gene_indices[valid_mask] - 1)  # Ensure integer indices
    valid_genes <- gene_names[valid_mask]
    
    if (length(valid_indices) == 0) {
      cat("❌ No valid genes found in batch\n")
      return(NULL)
    }
    
    tryCatch({
      peak <- peakRAM({
        expr_mat <- reticulate::py$get_genes_batch(z, valid_indices, layer)
      })
      print(peak)
      expr_mat <- py_to_r(expr_mat)
      colnames(expr_mat) <- valid_genes
      cat("📏 Batch retrieved expression for", length(valid_genes), "genes\n")
      return(expr_mat)
    }, error = function(e) {
      cat("⚠️ Error loading batch expression:", e$message, "\n")
      cat("Run `reticulate::py_last_error()` for details.\n")
      return(NULL)
    })
  }
  
  get_obsm_fast <- function(key) {
    if (!(key %in% obsm_keys)) return(NULL)
    tryCatch({
      return(py_to_r(z$obsm[[key]][]))
    }, error = function(e) {
      cat("⚠️ Error loading obsm", key, ":", e$message, "\n")
      return(NULL)
    })
  }
  
  list(
    df = umap_coords,
    obs = NULL,
    var = gene_info$var_df,
    genes = gene_info$display_names,
    gene_info = gene_info,
    zarr_obj = z,
    seacell_df = seacell_df,
    obs_keys = obs_columns,
    available_layers = layer_keys,
    obsm_keys = obsm_keys,
    file_path = file_path,
    is_zarr = TRUE,
    get_obs_column = get_obs_column_fast,
    get_gene_expression = get_gene_expression_fast,
    get_genes_expression_batch = get_genes_expression_batch,
    get_obsm = get_obsm_fast,
    # pathways = pathways,
    # pathway_map = pathway_map,  # Add mapping of full to display names
    n_cells = n_cells,
    n_genes = n_genes
  )
}


# Helper function to process heatmap data
process_heatmap_data <- function(mat, score_type, cell_group_mode, cluster_by) {
    # Fill missing values
    mat[is.na(mat)] <- 0

    # Scale or log
    if (score_type == "z_score") {
        mat <- scale(mat)
        mat[is.na(mat)] <- 0
    }

    group_info <- NULL
    # Cell grouping logic
    if (cell_group_mode == "categorical") {
        req(cluster_by)
        group_var <- values$lazy_data$get_obs_column(values$lazy_data$zarr_obj, cluster_by)
        req(length(group_var) == nrow(mat))
        
        group_levels <- sort(unique(group_var))
        group_var <- factor(group_var, levels = group_levels)
        
        mat <- sapply(colnames(mat), function(g) tapply(mat[, g], group_var, mean, na.rm = TRUE))
        mat <- t(mat)
    } else if (cell_group_mode == "none") {
        req(cluster_by)
        group_var <- values$lazy_data$get_obs_column(values$lazy_data$zarr_obj, cluster_by)
        req(length(group_var) == nrow(mat))
        
        group_levels <- sort(unique(group_var))
        group_var <- factor(group_var, levels = group_levels)
        o <- order(group_var)
        
        mat <- t(mat[o, , drop = FALSE])
        group_info <- group_var[o]
    }

    list(mat = mat, group_info = group_info)
}

render_heatmap <- function(data, title) {
mat <- data$mat

if (!is.null(data$group_info)) {
    groups <- levels(data$group_info)
    library(RColorBrewer)

    # Base palette (Set1 has up to 9 distinct colors)
    max_colors <- 9
    base_pal <- brewer.pal(max_colors, "Set1")

    # Suppose groups is your vector of category names
    n <- length(groups)

    # If more than 9 groups, extend colors using colorRampPalette
    pal_colors <- if (n <= max_colors) {
    base_pal[1:n]
    } else {
    colorRampPalette(base_pal)(n)
    }

    pal <- setNames(pal_colors, groups)

    ha <- HeatmapAnnotation(
    Group = data$group_info,
    col = list(Group = pal),
    show_annotation_name = FALSE
    )

    Heatmap(
    mat,
    name = "expression",
    top_annotation = ha,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_column_names = FALSE,
    column_title = title,
    use_raster = TRUE,
    raster_device = "CairoPNG",
    raster_quality = 3,
    column_split = data$group_info,
    row_names_max_width = unit(8, "cm"),
    row_names_gp = gpar(fontsize = 10)
    )
} else {
    Heatmap(
    mat,
    name = "expression",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    column_title = title,
    use_raster = TRUE,
    raster_device = "CairoPNG",
    raster_quality = 3,
    row_names_max_width = unit(8, "cm"),
    row_names_gp = gpar(fontsize = 10)
    )
}
}

load_umap_fast <- function(z) {
  umap_keys <- c("X_umap_normal", "X_umap", "X_umap_magic")
  umap <- NULL
  umap_key_used <- NULL

  py_list <- reticulate::import_builtins()$list
  obsm_keys <- reticulate::py_to_r(py_list(z[["obsm"]]$keys()))

  for (k in umap_keys) {
    if (k %in% obsm_keys) {
      umap <- reticulate::py_to_r(z[["obsm"]][[k]][])
      umap_key_used <- k
      break
    }
  }
  
  if (is.null(umap)) stop("❌ No UMAP coordinates found.")
  
  df <- as.data.frame(umap)
  colnames(df) <- c("x", "y")
  # Try MAGIC UMAP if available
  if ("X_umap_magic" %in% obsm_keys && umap_key_used != "X_umap_magic") {
    umap_magic <- reticulate::py_to_r(z[["obsm"]][["X_umap_magic"]][] )
    df$magic_x <- umap_magic[, 1]
    df$magic_y <- umap_magic[, 2]
    cat("✨ Loaded MAGIC UMAP coordinates\n")
  } else if (umap_key_used == "X_umap_magic") {
    df$magic_x <- df$x
    df$magic_y <- df$y
    cat("✨ Using main UMAP as MAGIC coordinates\n")
  } else {
    df$magic_x <- NA
    df$magic_y <- NA
  }
  
  return(df)
}