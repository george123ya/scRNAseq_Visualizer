# ==============================================================================
# Enhanced scRNA-seq Interactive Visualizer
# ==============================================================================
# A Shiny application for visualizing single-cell RNA sequencing
# data with dynamic gene expression, multi-annotation support, and advanced
# analysis tools including MAGIC imputation, differential expression, and
# differential abundance testing.
#
# Features:
# - Interactive scatterplot with regl-scatterplot
# - Dynamic gene expression visualization
# - Multiple annotation layers (clusters, cell types, conditions, etc.)
# - MAGIC gene imputation support
# - SEACell metacell visualization
# - Differential gene expression analysis
# - Differential abundance testing with Milo
# - Gene set signature scoring
# ==============================================================================

# Download demo .h5ad file if not present

# if (!file.exists("processed_data/relaxed_epdsc_annotated_data.h5")) {
#   dir.create("processed_data", showWarnings = FALSE)
#   message("📦 .h5 file not found, downloading from Google Drive...")
  
#   # Use gdown with the FILE ID
#   system("gdown https://drive.google.com/uc?id=1DwykyMuvohpo1aFQsQerlGsmLl-IiHpy -O processed_data/relaxed_epdsc_annotated_data.h5")
# }

# Load required libraries
suppressPackageStartupMessages({
  library(shiny)
  # library(shinydashboard)
  # library(DT)
  library(tidyr)
  library(base64enc)
  library(viridisLite)
  library(colorspace)
  library(reticulate)
  library(dplyr)
})

# options(shiny.host = "0.0.0.0")
# options(shiny.port = 3838)

# Source global configurations and helper functions
source("global.R")

# Configure Python environment for scanpy/anndata

# use_condaenv("sc_rna_env_python2", required = TRUE)
# use_condaenv("shiny_app_env", conda = "/opt/conda/bin/conda", required = TRUE)
use_condaenv("shiny_app_env", required = TRUE)


# ==============================================================================
# PYTHON INTERFACE SETUP
# ==============================================================================

# Import Python libraries for h5ad file handling
# anndata <- import("anndata", delay_load = TRUE)
# scanpy <- import("scanpy", delay_load = TRUE)

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

#' Initialize lazy-loaded h5ad file with HDF5 support
#' 
#' @param file_path Path to the .h5ad file
#' @return List containing metadata, coordinates, and AnnData object reference
#' @description Efficiently loads metadata while keeping expression data accessible
init_lazy_h5ad_improved <- function(file_path) {
  
  print(file_path)

  cat("🚀 Initializing improved lazy-loaded .h5ad file:", file_path, "\n")
  
  tryCatch({
    print("hey")
    # Try experimental lazy loading first (for Zarr-backed files)
    tryCatch({
      print("huy")
      anndata_exp <- import("anndata.experimental", delay_load = TRUE)
      ad_lazy <- anndata_exp$read_lazy(file_path)
      cat("✅ Using experimental lazy loading (Zarr-backed)\n")
      return(process_lazy_anndata(ad_lazy, file_path, is_zarr = TRUE))
    }, error = function(e) {
      cat("⚠️  Experimental lazy loading failed:", e$message, "\n")
      cat("🔄 Trying standard approach with selective loading...\n")
    })
    
    # Standard AnnData approach with selective loading strategy
    anndata <- import("anndata", delay_load = TRUE)
    h5py <- import("h5py", delay_load = TRUE)
    
    # First, try to read only metadata using h5py directly
    tryCatch({
      # Open HDF5 file to inspect structure without loading matrices
      h5_file <- h5py$File(file_path, "r")
      
      cat("📋 Inspecting HDF5 structure...\n")
      
      # Check available keys
      main_keys <- py_to_r(list(h5_file$keys()))
      cat("🔑 Main keys:", paste(main_keys, collapse = ", "), "\n")
      
      # Try to get basic info without loading full matrices
      if ("obs" %in% main_keys) {
        n_obs <- h5_file[["obs"]]$shape[0]
        cat("📊 Found", n_obs, "observations\n")
      }
      
      if ("var" %in% main_keys) {
        n_vars <- h5_file[["var"]]$shape[0] 
        cat("🧬 Found", n_vars, "variables\n")
      }
      
      h5_file$close()
      
    }, error = function(e) {
      cat("⚠️  Direct HDF5 inspection failed:", e$message, "\n")
    })
    
    # Load with standard AnnData but implement our own lazy strategy
    cat("📖 Loading AnnData object with selective strategy...\n")
    ad_obj <- anndata$read_h5ad(file_path)
    
    return(process_lazy_anndata(ad_obj, file_path, is_zarr = FALSE))
    
  }, error = function(e) {
    cat("❌ All loading methods failed:", e$message, "\n")
    stop("Could not initialize lazy loading: ", e$message)
  })
}

#' Process AnnData object for lazy-style access
#' 
#' @param ad_obj AnnData object (lazy or regular)
#' @param file_path Original file path
#' @param is_zarr Whether this is a Zarr-backed lazy object
#' @return Processed lazy data structure
process_lazy_anndata <- function(ad_obj, file_path, is_zarr = NULL) {
  pd <- import("pandas")
  
  # Auto-detect Zarr from filename if not provided
  if (is.null(is_zarr)) {
    is_zarr <- grepl("\\.zarr$", file_path, ignore.case = TRUE)
  }
  cat("ℹ️  File detected as", ifelse(is_zarr, "Zarr-backed", "HDF5-backed"), "\n")
  
  # Get obs shape without loading full obs
  shape <- ad_obj$obs$shape
  n_cells <- shape[[1]]
  n_cols <- shape[[2]]
  cat("📋 Observation metadata shape:", n_cells, "cells,", n_cols, "columns\n")
  
  # Get obs column names lazily
  obs_cols <- py_to_r(ad_obj$obs$columns$to_list())
  cat("🔑 Available obs columns:", paste(obs_cols, collapse = ", "), "\n")
  
  # Helper to get single obs column safely as R vector
  get_obs_column <- function(name) {
    if (!(name %in% obs_cols)) return(NULL)
    col_py <- ad_obj$obs[[name]]
    # Convert lazy-backed col to pandas Series (materialize)
    col_pd <- pd$Series(col_py)
    return(py_to_r(col_pd))
  }
  
  # Load var info and genes fully
  var_df <- py_to_r(ad_obj$var)
  gene_names <- as.character(py_to_r(ad_obj$var_names$to_list()))
  cat("🧬 Found", length(gene_names), "genes\n")
  
  # Get obsm keys
  obsm_keys <- py_to_r(ad_obj$obsm_keys())
  cat("🗝️ Available obsm keys:", paste(obsm_keys, collapse = ", "), "\n")
  
  # Load UMAP coordinates
  umap_keys <- c("X_umap_normal", "X_umap", "X_umap_magic")
  umap <- NULL
  umap_key_used <- NULL
  for (k in umap_keys) {
    if (k %in% obsm_keys) {
      umap_py <- ad_obj$obsm[[k]]
      if (reticulate::py_has_attr(umap_py, "compute")) {
        umap_py <- umap_py$compute()
      }
      umap <- py_to_r(umap_py)
      umap_key_used <- k
      cat("📍 Loaded UMAP coordinates from:", k, "\n")
      break
    }
  }
  if (is.null(umap)) stop("❌ No UMAP coordinates found.")
  
  df <- as.data.frame(umap)
  colnames(df) <- c("x", "y")
  
  # Helper to get obs columns by possible names, returns factor vector or NULL
  get_obs_col <- function(possible_names) {
    for (nm in possible_names) {
      val <- get_obs_column(nm)
      if (!is.null(val)) {
        cat("📊 Found annotation:", nm, "\n")
        return(as.factor(val))
      }
    }
    return(NULL)
  }
  
  # Extract common annotations
  df$cluster <- get_obs_col(c("pheno_100", "leiden", "cluster_feature", "Cluster_Fine"))
  df$cell_type <- get_obs_col(c("celltype", "Annotation"))
  df$condition <- get_obs_col(c("Condition", "Sample"))
  df$time <- get_obs_col(c("Time"))
  
  # Extract MAGIC UMAP coordinates if available
  if ("X_umap_magic" %in% obsm_keys) {
    if (umap_key_used != "X_umap_magic") {
      umap_magic <- py_to_r(ad_obj$obsm[["X_umap_magic"]])
      df$magic_x <- umap_magic[, 1]
      df$magic_y <- umap_magic[, 2]
      cat("✨ Loaded MAGIC UMAP coordinates\n")
    } else {
      df$magic_x <- df$x
      df$magic_y <- df$y
      cat("✨ Using main UMAP as MAGIC coordinates\n")
    }
  } else {
    df$magic_x <- NA
    df$magic_y <- NA
    cat("⚠️ MAGIC UMAP coordinates not found\n")
  }
  
  # Expression layers
  layers_py <- ad_obj$layers$keys()
  available_layers <- reticulate::py_to_r(reticulate::import_builtins()$list(layers_py))
  cat("📋 Available layers:", paste(available_layers, collapse = ", "), "\n")
  
  # SEACells summary
  seacell_df <- NULL
  uns_keys <- py_to_r(ad_obj$uns_keys())
  if ("SEACells_summary" %in% uns_keys) {
    tryCatch({
      seacell_summary <- ad_obj$uns[["SEACells_summary"]]
      seacell_umap <- py_to_r(seacell_summary$obsm[["X_umap"]])
      seacell_df <- as.data.frame(seacell_umap)
      colnames(seacell_df) <- c("UMAP_1", "UMAP_2")
      seacell_df$SEACell <- rownames(seacell_df)
      seacell_df$cell_type <- py_to_r(seacell_summary$obs[["cluster_feature"]])
      cat("🔬 Loaded SEACells metacells:", nrow(seacell_df), "\n")
    }, error = function(e) {
      cat("⚠️ Could not load SEACells data:", e$message, "\n")
    })
  }
  
  memory_strategy <- ifelse(is_zarr, "true_lazy", "selective_loading")
  if (!is_zarr) {
    cat("💾 Using selective loading strategy for HDF5 file\n")
  }
  
  cat("✅ Successfully processed:", nrow(df), "cells with", memory_strategy, "strategy\n\n")
  
  list(
    df = df,
    obs = NULL,          # Don't keep full obs loaded (too large, use get_obs_column)
    var = var_df,
    genes = gene_names,
    ad_obj = ad_obj,
    seacell_df = seacell_df,
    obs_keys = obs_cols,
    available_layers = available_layers,
    obsm_keys = obsm_keys,
    file_path = file_path,
    is_zarr = is_zarr,
    memory_strategy = memory_strategy,
    get_obs_column = get_obs_column  # export this function for user to get columns lazily
  )
}



#' Extract gene expression data with smart caching
#' 
#' @param gene_names Character vector of gene names to extract
#' @param lazy_data List containing processed AnnData data
#' @param use_layer Character, which expression layer to use
#' @param use_cache Whether to use expression cache
#' @return List with encoded expression data and statistics
extract_gene_data_smart <- function(gene_names, lazy_data, use_layer = "auto", use_cache = TRUE) {
  if (length(gene_names) == 0) {
    return(list(
      genes = character(0), 
      data = character(0),
      magic_data = character(0), 
      ranges = list(), 
      magic_ranges = list()
    ))
  }

  tryCatch({
    # DEBUG: print classes before coercion
    cat("DEBUG: class of gene_names before coercion:", class(gene_names), "\n")
    cat("DEBUG: class of lazy_data$genes before coercion:", class(lazy_data$genes), "\n")

    # Force full coercion to character vector
    gene_names_char <- as.character(gene_names)
    genes_vec <- lazy_data$genes
    
    # If genes_vec is python object, convert to R character vector forcibly
    if (inherits(genes_vec, "python.builtin.object")) {
      genes_vec <- reticulate::py_to_r(genes_vec)
    }
    genes_vec_char <- as.character(genes_vec)

    # Now match
    gene_indices <- match(gene_names_char, genes_vec_char)
    valid_mask <- !is.na(gene_indices)

    if (!any(valid_mask)) {
      cat("⚠️  No valid genes found\n")
      return(list(
        genes = character(0), 
        data = character(0),
        magic_data = character(0), 
        ranges = list(), 
        magic_ranges = list()
      ))
    }

    valid_genes <- gene_names_char[valid_mask]
    valid_indices <- gene_indices[valid_mask]

    # Continue rest of function ...
    cat("🔍 Smart extracting", length(valid_genes), "genes:", paste(head(valid_genes, 5), collapse = ", "), 
        if(length(valid_genes) > 5) "..." else "", "\n")
    
    # Determine which layer to use
    if (use_layer == "auto") {
      layers_priority <- c("log1p_data", "raw", "lognorm_pseudocount.1")
      use_layer <- NULL
      
      for (layer in layers_priority) {
        if (layer %in% lazy_data$available_layers) {
          use_layer <- layer
          cat("📈 Auto-selected expression layer:", layer, "\n")
          break
        }
      }
      
      if (is.null(use_layer)) {
        use_layer <- "X"
        cat("📈 Using main expression matrix (X)\n")
      }
    }

    
    # Extract expression data efficiently
    if (lazy_data$is_zarr && lazy_data$memory_strategy == "true_lazy") {
      # True lazy loading for Zarr files
      cat("⚡ True lazy extraction from Zarr\n")
      expr_py <- NULL
      if (use_layer == "X") {
        expr_py <- lazy_data$ad_obj$X
      } else {
        expr_py <- lazy_data$ad_obj$layers[[use_layer]]
      }
      expr_subset_py <- expr_py[, valid_indices - 1]  # python 0-based

      # If lazy/dask array, compute first
      if (reticulate::py_has_attr(expr_subset_py, "compute")) {
        expr_subset_py <- expr_subset_py$compute()
      }

      expr_subset <- py_to_r(expr_subset_py)
    } else {
      # Selective extraction from loaded object
      cat("📊 Selective extraction from loaded object\n")
      if (use_layer == "X") {
        if (inherits(lazy_data$ad_obj$X, "dgRMatrix")) {
          expr_subset <- lazy_data$ad_obj$X[valid_indices, , drop = FALSE]
          expr_subset <- Matrix::t(expr_subset)
        } else if (inherits(lazy_data$ad_obj$X, "dgCMatrix")) {
          expr_subset <- lazy_data$ad_obj$X[, valid_indices, drop = FALSE]
        } else {
          expr_py <- lazy_data$ad_obj$X[, valid_indices - 1]
          if (reticulate::py_has_attr(expr_py, "compute")) {
            expr_py <- expr_py$compute()
          }
          expr_subset <- py_to_r(expr_py)
        }
      } else {
        expr_py <- lazy_data$ad_obj$layers[[use_layer]][, valid_indices - 1]
        if (reticulate::py_has_attr(expr_py, "compute")) {
          expr_py <- expr_py$compute()
        }
        expr_subset <- py_to_r(expr_py)
      }
    }


    # Convert to dense matrix if sparse
    if (inherits(expr_subset, c("dgRMatrix", "dgCMatrix"))) {
      expr_matrix <- as.matrix(expr_subset)
    } else {
      expr_matrix <- as.matrix(expr_subset)
    }
    
    # Ensure correct dimensions
    if (nrow(expr_matrix) != nrow(lazy_data$df)) {
      if (ncol(expr_matrix) == nrow(lazy_data$df)) {
        expr_matrix <- t(expr_matrix)
        cat("🔄 Transposed expression matrix to correct orientation\n")
      }
    }
    
    cat("📊 Final expression matrix:", nrow(expr_matrix), "cells ×", ncol(expr_matrix), "genes\n")
    
    # Calculate gene expression statistics
    gene_ranges <- calculate_gene_ranges(expr_matrix, valid_genes)
    
    # Extract MAGIC data if available
    magic_data <- ""
    magic_ranges <- list()
    
    if ("MAGIC_imputed_data" %in% lazy_data$available_layers) {
      tryCatch({
        cat("✨ Extracting MAGIC imputed data...\n")
        magic_subset <- py_to_r(lazy_data$ad_obj$layers[["MAGIC_imputed_data"]][, valid_indices - 1])
        magic_matrix <- as.matrix(magic_subset)
        
        # Ensure correct orientation for MAGIC data too
        if (nrow(magic_matrix) != nrow(lazy_data$df)) {
          if (ncol(magic_matrix) == nrow(lazy_data$df)) {
            magic_matrix <- t(magic_matrix)
          }
        }
        
        # Encode MAGIC data
        magic_data_vector <- as.vector(magic_matrix)
        binary_magic_data <- writeBin(as.numeric(magic_data_vector), raw())
        magic_data <- base64enc::base64encode(binary_magic_data)
        
        # Calculate MAGIC ranges
        magic_ranges <- calculate_gene_ranges(magic_matrix, valid_genes)
        
        cat("✨ MAGIC data extracted successfully\n")
        
      }, error = function(e) {
        cat("⚠️  Could not extract MAGIC data:", e$message, "\n")
        magic_data <- ""
        magic_ranges <- list()
      })
    }
    
    # Encode expression data for JavaScript
    gene_data_vector <- as.vector(expr_matrix)
    binary_data <- writeBin(as.numeric(gene_data_vector), raw())
    encoded_data <- base64enc::base64encode(binary_data)
    
    cat("💾 Encoded data size:", nchar(encoded_data), "characters\n")
    cat("💾 MAGIC data size:", nchar(magic_data), "characters\n")
    
    return(list(
      genes = valid_genes,
      data = encoded_data,
      magic_data = magic_data,
      ranges = gene_ranges,
      magic_ranges = magic_ranges,
      nrows = nrow(expr_matrix),
      ncols = ncol(expr_matrix),
      layer_used = use_layer,
      memory_strategy = lazy_data$memory_strategy
    ))
    
  }, error = function(e) {
    cat("❌ Error in smart gene extraction:", e$message, "\n")
    return(list(
      genes = character(0),
      data = character(0),
      magic_data = character(0),
      ranges = list(),
      magic_ranges = list()
    ))
  })
}


#' Helper function to calculate gene expression ranges
#' 
#' @param expr_matrix Expression matrix
#' @param gene_names Gene names
#' @return List of ranges for each gene
calculate_gene_ranges <- function(expr_matrix, gene_names) {
  gene_ranges <- list()
  
  if (ncol(expr_matrix) > 1) {
    col_mins <- apply(expr_matrix, 2, min, na.rm = TRUE)
    col_maxs <- apply(expr_matrix, 2, max, na.rm = TRUE)
    col_means <- colMeans(expr_matrix, na.rm = TRUE)
    
    for (i in seq_along(gene_names)) {
      gene_ranges[[i]] <- list(
        min = col_mins[i],
        max = col_maxs[i],
        mean = col_means[i]
      )
    }
  } else if (ncol(expr_matrix) == 1) {
    vals <- expr_matrix[, 1]
    gene_ranges[[1]] <- list(
      min = min(vals, na.rm = TRUE),
      max = max(vals, na.rm = TRUE),
      mean = mean(vals, na.rm = TRUE)
    )
  }
  
  names(gene_ranges) <- gene_names
  return(gene_ranges)
}

#' Get comprehensive dataset information
#' 
#' @param lazy_data List containing processed AnnData data
#' @return List with detailed dataset information
get_dataset_info_smart <- function(lazy_data) {
  # Safely get available layers
  available_layers <- tryCatch({
    if (is.null(lazy_data$available_layers)) {
      character(0)
    } else if (is.list(lazy_data$available_layers)) {
      as.character(unlist(lazy_data$available_layers))
    } else {
      as.character(lazy_data$available_layers)
    }
  }, error = function(e) {
    cat("⚠️ Error getting available layers:", e$message, "\n")
    character(0)
  })
  
  # Safely get obsm keys
  obsm_keys <- tryCatch({
    if (is.null(lazy_data$obsm_keys)) {
      character(0)
    } else if (is.list(lazy_data$obsm_keys)) {
      as.character(unlist(lazy_data$obsm_keys))
    } else {
      as.character(lazy_data$obsm_keys)
    }
  }, error = function(e) {
    cat("⚠️ Error getting obsm keys:", e$message, "\n")
    character(0)
  })
  
  # Safely check for MAGIC coordinates
  has_magic_coords <- tryCatch({
    !all(is.na(lazy_data$df$magic_x))
  }, error = function(e) {
    cat("⚠️ Error checking MAGIC coords:", e$message, "\n")
    FALSE
  })
  
  # Safely check for MAGIC expression
  has_magic_expr <- tryCatch({
    "MAGIC_imputed_data" %in% available_layers
  }, error = function(e) {
    cat("⚠️ Error checking MAGIC expression:", e$message, "\n")
    FALSE
  })
  
  # Safely check for SEACells
  has_seacells <- tryCatch({
    !is.null(lazy_data$seacell_df) && nrow(lazy_data$seacell_df) > 0
  }, error = function(e) {
    cat("⚠️ Error checking SEACells:", e$message, "\n")
    FALSE
  })

  get_local_size_mb <- function(path) {
    if (dir.exists(path)) {
      # Recursively sum all file sizes
      sum(file.info(list.files(path, recursive = TRUE, full.names = TRUE))$size, na.rm = TRUE) / (1024^2)
    } else if (file.exists(path)) {
      # Normal file
      file.size(path) / (1024^2)
    } else {
      0
    }
  }

  get_cloud_size_mb <- function(base_url, zarr_name) {
    total_size <- 0
    marker <- NULL
    
    repeat {
      url <- paste0(base_url, "?prefix=", zarr_name, if (!is.null(marker)) paste0("&marker=", URLencode(marker, reserved = TRUE)) else "")
      
      doc <- xml2::read_xml(url)
      ns <- xml2::xml_ns(doc)
      
      # Add sizes from <Size>
      sizes <- as.numeric(xml2::xml_text(xml2::xml_find_all(doc, ".//d1:Size", ns)))
      total_size <- total_size + sum(sizes, na.rm = TRUE)
      
      truncated <- xml2::xml_text(xml2::xml_find_first(doc, ".//d1:IsTruncated", ns)) == "true"
      if (!truncated) break
      
      marker <- xml2::xml_text(xml2::xml_find_first(doc, ".//d1:NextMarker", ns))
    }
    
    total_size / (1024^2)  # in MB
  }

  
  # Safely get file size
  get_dataset_size_mb <- function(path_or_url) {
    tryCatch({
      if (grepl("^https://storage\\.googleapis\\.com/", path_or_url)) {
        # Cloud Zarr
        base_url <- "https://storage.googleapis.com/scrna-seqbrowser/"
        zarr_name <- sub(".*/", "", path_or_url)
        get_cloud_size_mb(base_url, zarr_name)
      } else {
        # Local file or folder
        get_local_size_mb(path_or_url)
      }
    }, error = function(e) {
      cat("⚠️ Error getting file size:", e$message, "\n")
      0
    })
  }

  file_size_mb <- get_dataset_size_mb(lazy_data$file_path)
  file_size_mb <- round(file_size_mb, 2)

  
  # Safely get memory strategy
  memory_strategy <- tryCatch({
    if (!is.null(lazy_data$memory_strategy)) {
      as.character(lazy_data$memory_strategy)
    } else {
      "unknown"
    }
  }, error = function(e) {
    cat("⚠️ Error getting memory strategy:", e$message, "\n")
    "unknown"
  })
  
  # Safely get is_zarr flag
  is_zarr <- tryCatch({
    if (!is.null(lazy_data$is_zarr)) {
      as.logical(lazy_data$is_zarr)
    } else {
      FALSE
    }
  }, error = function(e) {
    cat("⚠️ Error getting zarr flag:", e$message, "\n")
    FALSE
  })
  
  info <- list(
    n_cells = nrow(lazy_data$df),
    n_genes = length(lazy_data$genes),
    available_layers = available_layers,
    obsm_keys = obsm_keys,
    has_magic_coords = has_magic_coords,
    has_magic_expr = has_magic_expr,
    has_seacells = has_seacells,
    file_size_mb = file_size_mb,
    memory_strategy = memory_strategy,
    is_zarr = is_zarr,
    file_type = if(is_zarr) "Zarr-backed" else "HDF5"
  )
  
  cat("✅ Dataset info generated successfully\n")
  return(info)
}


# Fixed server function with proper initialization

server <- function(input, output, session) {
  # Reactive values to store processed data
  values <- reactiveValues(
    processed_data = NULL,
    annotation_names = NULL,
    lazy_data = NULL,
    current_dataset = NULL,
    dataset_info = NULL
  )
  
  observeEvent(input$activateMAGIC, {
    session$sendCustomMessage("toggleMAGIC", input$activateMAGIC)
  })
  
  # Process data when lazy_data is available
  observeEvent(values$lazy_data, {
    req(values$lazy_data)
    
    tryCatch({
      session$sendCustomMessage("showSpinner", TRUE)
      output$status <- renderText("Processing metadata...")

      # Use the loaded lazy_data
      raw_data <- values$lazy_data$df

      # print(head(raw_data))
      
      # Process the dataframe
      processed_result <- process_dataframe(raw_data)
      values$processed_data <- processed_result
      
      # Create color by choices
      annotation_names <- names(processed_result$annotations)
      values$annotation_names <- annotation_names

      # print(head(processed_result$data))
      
      # Get dataset info using improved function with error handling
      values$dataset_info <- tryCatch({
        get_dataset_info_smart(values$lazy_data)
      }, error = function(e) {
        cat("⚠️ Error getting dataset info:", e$message, "\n")
        # Return minimal info structure
        list(
          n_cells = nrow(values$lazy_data$df),
          n_genes = length(values$lazy_data$genes),
          available_layers = character(0),
          obsm_keys = character(0),
          has_magic_coords = FALSE,
          has_magic_expr = FALSE,
          has_seacells = FALSE,
          file_size_mb = 0,
          memory_strategy = "unknown",
          is_zarr = FALSE,
          file_type = "Unknown"
        )
      })
      
      # Encode data for JavaScript
      encoded_data <- encode_data(processed_result)
      n_cols <- ncol(processed_result$data) 

      # Send to JavaScript
      session$sendCustomMessage("updateData", list(
        base64 = encoded_data,
        annotationData = processed_result$annotations,
        numCols = n_cols,
        clusters = processed_result$annotations$cluster$names,
        colors = processed_result$cluster_colors,
        metacellColors = processed_result$metacell_colors,
        geneExprRanges = processed_result$gene_expr_ranges
      ))
      
      session$sendCustomMessage("showSpinner", FALSE)
      
      strategy_text <- if(values$lazy_data$memory_strategy == "true_lazy") "truly lazy loaded" else "selectively loaded"
      output$status <- renderText(paste0("Rendered ", nrow(raw_data), " cells with ", 
                                        length(annotation_names), " annotations (", strategy_text, ")"))
      
    }, error = function(e) {
      session$sendCustomMessage("showSpinner", FALSE)
      showNotification(
        paste("Error processing data:", e$message), 
        type = "error", 
        duration = 10
      )
      output$status <- renderText(paste("Error processing data:", e$message))
      cat("❌ Error in data processing observer:", e$message, "\n")
    })
  })

  # UPDATED: Enhanced dataset loading with improved lazy loading
  observeEvent(input$dataset, {
    req(input$dataset)
    
    # Skip if placeholder values
    if (input$dataset %in% c("", "❌ Folder not found", "❌ No .h5ad files found")) {
      return()
    }
    
    # Check if file actually exists
    # if (!file.exists(input$dataset)) {
    #   showNotification(
    #     paste("File not found:", input$dataset), 
    #     type = "error", 
    #     duration = 5
    #   )
    #   return()
    # }
    
    # Show loading notification
    showNotification(
      paste("Initializing smart loading for:", basename(input$dataset)), 
      type = "message", 
      duration = 3
    )
    
    output$status <- renderText("Initializing smart loading...")
    session$sendCustomMessage("showSpinner", TRUE)
    
    tryCatch({
      # Use improved initialization
      values$current_dataset <- input$dataset
      print(input$dataset)
      values$lazy_data <- init_lazy_h5ad_improved(input$dataset)
      
      # Determine success message based on strategy
      strategy_msg <- if(values$lazy_data$is_zarr) {
        "✅ True lazy loading active (Zarr-backed)"
      } else {
        "✅ Smart selective loading active (HDF5)"
      }
      
      showNotification(strategy_msg, type = "message", duration = 4)
      
      cat("✅ Initialized smart loading for:", basename(input$dataset), "\n")
      cat("   Path:", input$dataset, "\n")
      cat("   Cells:", nrow(values$lazy_data$df), "\n")
      cat("   Genes:", length(values$lazy_data$genes), "\n")
      cat("   Strategy:", values$lazy_data$memory_strategy, "\n")
      cat("   File type:", if(values$lazy_data$is_zarr) "Zarr" else "HDF5", "\n")
      
    }, error = function(e) {
      session$sendCustomMessage("showSpinner", FALSE)
      showNotification(
        paste("❌ Error initializing loading:", e$message), 
        type = "error", 
        duration = 10
      )
      cat("❌ Error initializing smart loading:", e$message, "\n")
      
      # Clear any partial data
      values$lazy_data <- NULL
      values$processed_data <- NULL
      values$dataset_info <- NULL
      output$status <- renderText("Error loading dataset")
    })
  })

  library(xml2)

  get_top_level_zarr <- function(base_url) {
    zarr_files <- c()
    marker <- NULL
    
    repeat {
      url <- if (!is.null(marker)) {
        paste0(base_url, "?marker=", URLencode(marker, reserved = TRUE))
      } else {
        base_url
      }
      
      doc <- read_xml(url)
      ns <- xml_ns(doc)
      keys <- xml_text(xml_find_all(doc, ".//d1:Key", ns))
      
      # Keep only base zarr names
      zarr_base <- sub("^(.+?\\.zarr).*", "\\1", keys)
      zarr_base <- zarr_base[grepl("\\.zarr$", zarr_base)]
      
      zarr_files <- unique(c(zarr_files, zarr_base))
      
      truncated <- xml_text(xml_find_first(doc, ".//d1:IsTruncated", ns)) == "true"
      if (!truncated) break
      
      marker <- xml_text(xml_find_first(doc, ".//d1:NextMarker", ns))
    }
    
    zarr_files
  }

  output$datasetUI <- renderUI({
    folder_path <- "./processed_data"
    local_files <- character()
    
    # Local file scan
    if (dir.exists(folder_path)) {
      local_files <- list.files(folder_path, pattern = "\\.(zarr|h5ad|h5)$", 
                                full.names = TRUE)
    }
    
    # Cloud Zarr list
    cloud_files <- tryCatch(
      get_top_level_zarr("https://storage.googleapis.com/scrna-seqbrowser/"),
      error = function(e) character()
    )
    
    # Prepare local choices
    local_choices <- NULL
    if (length(local_files) > 0) {
      display_names <- gsub(paste0("^", folder_path, "/?"), "", local_files)
      display_names <- gsub("\\.(h5ad|h5)$", "", basename(display_names))
      local_choices <- setNames(local_files, display_names)
    }
    
    # Prepare cloud choices — friendly label, full URL value
    cloud_choices <- NULL
    if (length(cloud_files) > 0) {
      cloud_choices <- setNames(
        paste0("https://storage.googleapis.com/scrna-seqbrowser/", cloud_files),
        paste0("[Cloud] ", cloud_files)
      )
    }
    
    # Combine into grouped list
    choices <- list()
    if (!is.null(local_choices)) choices[["Local"]] <- local_choices
    if (!is.null(cloud_choices)) choices[["Cloud"]] <- cloud_choices
    
    if (length(choices) == 0) {
      return(
        div(
          selectInput("dataset", "Select Dataset:", 
                      choices = list("❌ No datasets found" = ""), 
                      selected = "")
        )
      )
    }

    print(choices)
    
    selectInput("dataset", "Select Dataset:", 
                choices = choices,
                selected = choices[[1]][1],
                width = "100%")
  })


  # UI: Gene search
  output$geneSearchUI <- renderUI({
    if (is.null(values$lazy_data)) {
      selectizeInput("geneSearch", "🔍 Search Gene:",
                    choices = list("Initialize dataset first..." = "loading"),
                    multiple = TRUE)
    } else {
      selectizeInput(
        "geneSearch",
        "🔍 Search Gene:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = paste0("Type gene names (", values$lazy_data$memory_strategy, ")..."),
          openOnFocus = FALSE,
          closeAfterSelect = TRUE,
          plugins = list("remove_button"),
          maxItems = 15
        )
      )
    }
  })

  # UPDATED: Gene search observer with smart extraction
  observeEvent(input$geneSearch, {
    req(values$lazy_data)
    
    if (!is.null(input$geneSearch) && length(input$geneSearch) > 0) {
      gene_list <- as.character(input$geneSearch)
      
      print(paste("🔍 Smart searching for genes:", paste(gene_list, collapse = ", ")))
      
      # Show loading status
      loading_msg <- paste("Loading expression data for", length(gene_list), "genes...")
      output$status <- renderText(loading_msg)
      session$sendCustomMessage("showSpinner", TRUE)
      
      # Use improved smart extraction
      gene_data <- extract_gene_data_smart(gene_list, values$lazy_data)
      
      session$sendCustomMessage("showSpinner", FALSE)
      
      print(paste("📊 Smart extracted data summary:"))
      print(paste("  - Genes found:", length(gene_data$genes)))
      print(paste("  - Normal data size:", nchar(gene_data$data)))
      print(paste("  - MAGIC data size:", nchar(gene_data$magic_data)))
      print(paste("  - Layer used:", gene_data$layer_used))
      print(paste("  - Strategy:", gene_data$memory_strategy))
      
      session$sendCustomMessage("geneSearchChange", list(
        genes = gene_list,
        expression_data = gene_data
      ))
      
      # Update status with strategy info
      strategy_info <- if(values$lazy_data$is_zarr) "(truly lazy)" else "(selective)"
      output$status <- renderText(paste0("Loaded expression for ", length(gene_data$genes), 
                                        " genes ", strategy_info))
      
    } else {
      session$sendCustomMessage("geneSearchChange", list(
        genes = character(0),
        expression_data = list(
          genes = character(0), 
          data = character(0), 
          magic_data = character(0),
          ranges = list(),
          magic_ranges = list()
        )
      ))
      
      if (!is.null(values$dataset_info)) {
        strategy_text <- paste0("(", values$dataset_info$memory_strategy, ")")
        output$status <- renderText(paste0("Ready - ", values$dataset_info$n_cells, " cells, ", 
                                          values$dataset_info$n_genes, " genes ", strategy_text))
      }
    }
  }, ignoreNULL = FALSE)

  # Update gene choices when data is loaded
  observeEvent(values$lazy_data, {
    req(values$lazy_data)
    
    updateSelectizeInput(
      session,
      "geneSearch",
      choices = values$lazy_data$genes,
      server = TRUE
    )
  })
  
  # Color by UI
  output$colorByUI <- renderUI({
    if (is.null(values$annotation_names)) {
      selectInput("colorBy", "Color cells by:", choices = list("Initialize dataset first..." = "loading"))
    } else {
      choices <- setNames(values$annotation_names, tools::toTitleCase(gsub("_", " ", values$annotation_names)))
      selectInput("colorBy", "Color cells by:", choices = choices, selected = choices[1])
    }
  })
  
  # Handle color by changes
  observeEvent(input$colorBy, {
    if (!is.null(input$colorBy) && input$colorBy != "loading") {
      session$sendCustomMessage("colorByChange", input$colorBy)
    }
  })

  # MAGIC toggle management
  observe({
    if (!is.null(input$colorBy)) {
        if (input$colorBy == "seacell") {
            print("Seacell selected, disabling MAGIC")
            disable("activateMAGIC")
        } else {
            enable("activateMAGIC")
        }
    }
  })

  # UPDATED: Enhanced data info output
  output$dataInfo <- renderUI({
    req(values$lazy_data, values$processed_data, values$dataset_info)
    
    info <- values$dataset_info
    result <- values$processed_data
    
    # Create strategy-specific messaging
    strategy_color <- if(info$memory_strategy == "true_lazy") "green" else "blue"
    strategy_icon <- if(info$memory_strategy == "true_lazy") "⚡" else "🧠"
    strategy_text <- if(info$memory_strategy == "true_lazy") "True Lazy Loading" else "Smart Selective Loading"
    
    info_text <- tagList(
      h5(paste("📋", info$file_type, "Dataset Info:")),
      p(paste("🔹 Total cells:", info$n_cells)),
      p(paste("🔹 Total genes:", info$n_genes)),
      p(paste("🔹 File size:", info$file_size_mb, "MB")),
      p(paste("🔹 Annotations:", length(result$annotations))),
      p(paste("🔹 MAGIC coordinates:", ifelse(info$has_magic_coords, "✓ Available", "❌ Not available"))),
      p(paste("🔹 MAGIC expression:", ifelse(info$has_magic_expr, "✓ Available", "❌ Not available"))),
      p(paste("🔹 SEACells:", ifelse(info$has_seacells, "✓ Available", "❌ Not available"))),
      hr(),
      p(strong(paste(strategy_icon, strategy_text, ":")), style = paste0("color: ", strategy_color, ";")),
      if(info$memory_strategy == "true_lazy") {
        div(
          p("• Expression data loaded on-demand", style = "font-size: 11px; color: green;"),
          p("• Minimal memory footprint", style = "font-size: 11px; color: green;"),
          p("• Zarr-backed lazy access", style = "font-size: 11px; color: green;")
        )
      } else {
        div(
          p("• Metadata loaded, expressions selective", style = "font-size: 11px; color: blue;"),
          p("• Smart matrix access", style = "font-size: 11px; color: blue;"),
          p("• HDF5 optimized loading", style = "font-size: 11px; color: blue;")
        )
      },
      hr(),
      if (length(info$available_layers) > 0) {
        div(
          strong("Available expression layers:"),
          br(),
          paste(info$available_layers, collapse = ", "),
          style = "font-size: 11px;"
        )
      }
    )
    
    div(style = "background: #f8f9fa; padding: 10px; border-radius: 5px; font-size: 12px;", info_text)
  })

  # UPDATED: Enhanced gene set signature with smart loading
  observeEvent(input$calc_gene_score, {
    req(values$lazy_data, input$gene_set_input)
    
    if (length(input$gene_set_input) == 0) {
      showNotification("Please select genes for signature calculation", type = "warning")
      return()
    }
    
    showNotification("Calculating gene set signature...", type = "message", duration = 2)
    output$status <- renderText("Calculating gene set signature...")
    session$sendCustomMessage("showSpinner", TRUE)
    
    tryCatch({
      # Extract expression data using smart method
      gene_data <- extract_gene_data_smart(input$gene_set_input, values$lazy_data)
      
      session$sendCustomMessage("showSpinner", FALSE)
      
      if (length(gene_data$genes) == 0) {
        showNotification("No valid genes found in dataset", type = "error")
        return()
      }
      
      # Decode the expression data
      binary_data <- base64enc::base64decode(gene_data$data)
      expr_vector <- readBin(binary_data, "numeric", n = gene_data$nrows * gene_data$ncols)
      expr_matrix <- matrix(expr_vector, nrow = gene_data$nrows, ncol = gene_data$ncols)
      
      # Calculate gene set signature (mean expression)
      if (ncol(expr_matrix) > 1) {
        signature_scores <- rowMeans(expr_matrix, na.rm = TRUE)
      } else {
        signature_scores <- expr_matrix[, 1]
      }
      
      # Encode signature scores
      signature_binary <- writeBin(as.numeric(signature_scores), raw())
      signature_encoded <- base64enc::base64encode(signature_binary)
      
      # Calculate signature statistics
      signature_range <- list(
        min = min(signature_scores, na.rm = TRUE),
        max = max(signature_scores, na.rm = TRUE),
        mean = mean(signature_scores, na.rm = TRUE)
      )
      
      # Send to JavaScript
      session$sendCustomMessage("updateGeneSignature", list(
        signature_data = signature_encoded,
        range = signature_range,
        genes_used = gene_data$genes,
        n_genes = length(gene_data$genes)
      ))
      
      strategy_info <- paste0("(", gene_data$memory_strategy, ")")
      output$status <- renderText(paste0("Gene signature calculated using ", length(gene_data$genes), 
                                        " genes ", strategy_info))
      showNotification(paste("Signature calculated using", length(gene_data$genes), "genes"), type = "message")
      
    }, error = function(e) {
      session$sendCustomMessage("showSpinner", FALSE)
      showNotification(paste("Error calculating signature:", e$message), type = "error", duration = 5)
      output$status <- renderText("Error calculating gene signature")
    })
  })

  # UPDATED: Enhanced differential expression with smart loading
  observeEvent(input$run_dge, {
    req(values$lazy_data, input$dge_group_by, input$dge_condition_1, input$dge_condition_2)
    
    if (input$dge_condition_1 == input$dge_condition_2) {
      showNotification("Please select different conditions for comparison", type = "warning")
      return()
    }
    
    showNotification("Running differential expression analysis...", type = "message", duration = 3)
    output$status <- renderText("Running differential expression analysis...")
    session$sendCustomMessage("showSpinner", TRUE)
    
    tryCatch({
      # Get cell indices for each condition
      group_col <- input$dge_group_by
      obs_data <- values$lazy_data$obs
      
      if (!group_col %in% names(obs_data)) {
        showNotification(paste("Column not found:", group_col), type = "error")
        session$sendCustomMessage("showSpinner", FALSE)
        return()
      }
      
      condition1_cells <- which(obs_data[[group_col]] == input$dge_condition_1)
      condition2_cells <- which(obs_data[[group_col]] == input$dge_condition_2)
      
      if (length(condition1_cells) == 0 || length(condition2_cells) == 0) {
        showNotification("No cells found for one or both conditions", type = "error")
        session$sendCustomMessage("showSpinner", FALSE)
        return()
      }
      
      cat("📊 DE analysis:", length(condition1_cells), "vs", length(condition2_cells), "cells\n")
      
      # Adaptive gene selection based on dataset size and strategy
      max_genes <- if(values$lazy_data$memory_strategy == "true_lazy") 500 else 200
      test_genes <- head(values$lazy_data$genes, max_genes)
      
      cat("🧬 Testing", length(test_genes), "genes with", values$lazy_data$memory_strategy, "strategy\n")
      
      # Extract expression data using smart method
      gene_data <- extract_gene_data_smart(test_genes, values$lazy_data)
      
      if (length(gene_data$genes) == 0) {
        showNotification("No gene expression data found", type = "error")
        session$sendCustomMessage("showSpinner", FALSE)
        return()
      }
      
      # Decode expression data
      binary_data <- base64enc::base64decode(gene_data$data)
      expr_vector <- readBin(binary_data, "numeric", n = gene_data$nrows * gene_data$ncols)
      expr_matrix <- matrix(expr_vector, nrow = gene_data$nrows, ncol = gene_data$ncols)
      colnames(expr_matrix) <- gene_data$genes
      
      # Perform differential expression analysis
      de_results <- perform_de_analysis(expr_matrix, condition1_cells, condition2_cells, gene_data$genes)
      
      # Display results
      output$de_results <- renderDT({
        datatable(de_results,
          options = list(pageLength = 10, scrollX = TRUE),
          selection = "single",
          caption = paste("DE Analysis:", input$dge_condition_1, "vs", input$dge_condition_2)
        ) %>%
        formatRound(c("mean_cond1", "mean_cond2", "log2fc"), 3) %>%
        formatSignif(c("pvalue", "padj"), 3)
      })
      
      session$sendCustomMessage("showSpinner", FALSE)
      strategy_info <- paste0("(", gene_data$memory_strategy, ")")
      output$status <- renderText(paste0("DE analysis complete: ", nrow(de_results), " genes tested ", strategy_info))
      showNotification("Differential expression analysis complete", type = "message")
      
    }, error = function(e) {
      session$sendCustomMessage("showSpinner", FALSE)
      showNotification(paste("Error in DE analysis:", e$message), type = "error", duration = 5)
      output$status <- renderText("Error in differential expression analysis")
    })
  })

  # Helper function for DE analysis
  perform_de_analysis <- function(expr_matrix, condition1_cells, condition2_cells, gene_names) {
    de_results <- data.frame(
      gene = gene_names,
      mean_cond1 = numeric(length(gene_names)),
      mean_cond2 = numeric(length(gene_names)),
      log2fc = numeric(length(gene_names)),
      pvalue = numeric(length(gene_names)),
      stringsAsFactors = FALSE
    )
    
    for (i in seq_along(gene_names)) {
      expr1 <- expr_matrix[condition1_cells, i]
      expr2 <- expr_matrix[condition2_cells, i]
      
      de_results$mean_cond1[i] <- mean(expr1, na.rm = TRUE)
      de_results$mean_cond2[i] <- mean(expr2, na.rm = TRUE)
      de_results$log2fc[i] <- log2((de_results$mean_cond2[i] + 0.001) / (de_results$mean_cond1[i] + 0.001))
      
      # Simple t-test
      if (var(expr1, na.rm = TRUE) > 0 || var(expr2, na.rm = TRUE) > 0) {
        tryCatch({
          test_result <- t.test(expr1, expr2)
          de_results$pvalue[i] <- test_result$p.value
        }, error = function(e) {
          de_results$pvalue[i] <- 1.0
        })
      } else {
        de_results$pvalue[i] <- 1.0
      }
    }
    
    # Adjust p-values and sort
    de_results$padj <- p.adjust(de_results$pvalue, method = "BH")
    de_results <- de_results[order(de_results$pvalue), ]
    
    return(de_results)
  }

  # Update gene set choices when data is loaded
  observeEvent(values$lazy_data, {
    req(values$lazy_data)
    
    updateSelectizeInput(
      session,
      "gene_set_input",
      choices = values$lazy_data$genes,
      server = TRUE
    )
  })

  # Update DGE condition choices
  observe({
    req(values$lazy_data, input$dge_group_by)
    
    if (input$dge_group_by %in% names(values$lazy_data$obs)) {
      unique_conditions <- unique(values$lazy_data$obs[[input$dge_group_by]])
      unique_conditions <- unique_conditions[!is.na(unique_conditions)]
      
      updateSelectInput(session, "dge_condition_1", 
                       choices = unique_conditions, 
                       selected = unique_conditions[1])
      updateSelectInput(session, "dge_condition_2", 
                       choices = unique_conditions, 
                       selected = if(length(unique_conditions) > 1) unique_conditions[2] else unique_conditions[1])
    }
  })

  # Selected points output
  output$selected_points <- renderPrint({
    if (is.null(input$selectedPoints) || length(input$selectedPoints) == 0) {
      "No cells selected. Use lasso tool to select cells."
    } else {
      paste0(
        "Selected cells: ", length(input$selectedPoints),
        "\nIndices (first 15): ", paste(head(input$selectedPoints, 15), collapse = ", "),
        if (length(input$selectedPoints) > 15) " ..."
      )
    }
  })

  # Helper to classify obs_keys
  categorical_vars <- reactive({
    req(values$lazy_data)
    keys <- values$lazy_data$obs_keys
    
    cats <- Filter(function(k) {
      col <- values$lazy_data$get_obs_column(k)
      uniq_vals <- unique(col)
      
      # Detect integer-like numeric
      is_integer_like <- is.numeric(col) && all(!is.na(col)) && all(col == floor(col))
      
      (is.factor(col) || is.character(col) || is_integer_like) &&
        length(uniq_vals) <= 20
    }, keys)
    
    sort(cats)
  })

  # Helper to get numeric obs columns
  numeric_obs_keys <- reactive({
    req(values$lazy_data)
    cols <- values$lazy_data$obs_keys
    
    nums <- Filter(function(k) {
      col <- values$lazy_data$get_obs_column(k)
      is_real_numeric <- is.numeric(col) && !all(col == floor(col))  # exclude integer-like
      length(unique(col)) > 1 && is_real_numeric
    }, cols)
    
    sort(nums)
  })

  # Single metric (for violin/box/hist)
  output$qc_metric_ui <- renderUI({
    selectInput("qc_metric", "QC Metric:",
      choices = numeric_obs_keys(),
      selected = numeric_obs_keys()[1]
    )
  })

  # Group by
  output$qc_group_by_ui <- renderUI({
    selectInput("qc_group_by", "Group by:",
      choices = categorical_vars(),
      selected = categorical_vars()[1]
    )
  })

  # Scatter X
  output$qc_x_metric_ui <- renderUI({
    selectInput("qc_x_metric", "X Metric:",
      choices = numeric_obs_keys(),
      selected = numeric_obs_keys()[1]
    )
  })

  # Scatter Y
  output$qc_y_metric_ui <- renderUI({
    selectInput("qc_y_metric", "Y Metric:",
      choices = numeric_obs_keys(),
      selected = numeric_obs_keys()[2]
    )
  })

  # Scatter color
  output$qc_color_by_ui <- renderUI({
    selectInput("qc_color_by", "Color by:",
      choices = c("None", categorical_vars()),
      selected = "None"
    )
  })
  
  # qc_main_plot

  # Reactive data for main plot — returns different data frames depending on plot type
  qc_main_data <- reactive({
    req(values$lazy_data, input$qc_plot_type)
    
    if (input$qc_plot_type == "scatter") {
      req(input$qc_x_metric, input$qc_y_metric)
      x <- as.numeric(values$lazy_data$get_obs_column(input$qc_x_metric))
      y <- as.numeric(values$lazy_data$get_obs_column(input$qc_y_metric))
      
      if (!is.null(input$qc_color_by) && input$qc_color_by != "None") {
        color <- as.factor(values$lazy_data$get_obs_column(input$qc_color_by))
        data.frame(X = x, Y = y, Color = color)
      } else {
        data.frame(X = x, Y = y)
      }
      
    } else {
      req(input$qc_metric, input$qc_group_by)
      metric <- as.numeric(values$lazy_data$get_obs_column(input$qc_metric))
      group <- as.factor(values$lazy_data$get_obs_column(input$qc_group_by))
      data.frame(Metric = metric, Group = group)
    }
  })

  # Single renderPlot that handles all plot types
  output$qc_main_plot <- renderPlot({
    dat <- qc_main_data()
    
    # Scatter plot branch
    if (input$qc_plot_type == "scatter") {
      # Downsample if too large
      if (nrow(dat) > 50000) {
        set.seed(123)
        dat <- dat[sample(nrow(dat), 50000), ]
      }
      
      p <- ggplot(dat, aes(x = X, y = Y))

      if ("Color" %in% colnames(dat)) {
        p <- p + aes(color = Color)
      }

      p <- p + geom_point(alpha = 0.5, size = input$qc_point_size) +
        labs(
          title = paste(input$qc_x_metric, "vs", input$qc_y_metric),
          x = input$qc_x_metric,
          y = input$qc_y_metric
        ) +
        theme_minimal()
      
      if (input$qc_log_scale) {
        p <- p + scale_x_log10() + scale_y_log10()
      }
      
    } else {
      # Downsample for jitter points only (violin/box)
      if (input$qc_show_points && nrow(dat) > 20000) {
        set.seed(123)
        dat_points <- dat[sample(nrow(dat), 20000), ]
      } else {
        dat_points <- dat
      }
      
      p <- ggplot(dat, aes(x = Group, y = Metric, fill = Group))
      
      if (input$qc_plot_type == "violin") {
        p <- p + geom_violin(trim = FALSE)
        if (input$qc_show_points) {
          p <- p + geom_jitter(data = dat_points, width = 0.2, alpha = 0.4, size = input$qc_point_size)
        }
      } else if (input$qc_plot_type == "box") {
        p <- p + geom_boxplot()
        if (input$qc_show_points) {
          p <- p + geom_jitter(data = dat_points, width = 0.2, alpha = 0.4, size = input$qc_point_size)
        }
      } else if (input$qc_plot_type == "histogram") {
        p <- ggplot(dat, aes(x = Metric, fill = Group)) +
          geom_histogram(position = "identity", alpha = 0.5, bins = 50)
      }
      
      if (input$qc_log_scale) {
        p <- p + scale_y_log10()
      }
      
      p <- p + theme_minimal() +
        labs(
          title = paste("QC Metric:", input$qc_metric),
          x = input$qc_group_by,
          y = input$qc_metric
        )
    }
    
    p
  })


  # Helper: get first existing obs column from a list of possible names
  get_first_available_column <- function(possible_names) {
    for (nm in possible_names) {
      col_data <- values$lazy_data$get_obs_column(nm)
      if (!is.null(col_data)) {
        return(col_data)
      }
    }
    return(NULL)
  }

  output$qc_stats_table <- renderTable({
    req(values$lazy_data)
    
    n_genes_col <- get_first_available_column(c("n_genes", "nGenes", "genes_count"))
    total_counts_col <- get_first_available_column(c("total_counts", "totalUMIs"))
    pct_mito_col <- get_first_available_column(c("pct_mito", "percent_mito", "mitochondrial_percent"))
    pct_ribo_col <- get_first_available_column(c("pct_ribo", "percent_ribo", "ribosomal_percent"))
    
    obs_list <- list()
    if (!is.null(n_genes_col)) obs_list[["n_genes"]] <- n_genes_col
    if (!is.null(total_counts_col)) obs_list[["total_counts"]] <- total_counts_col
    if (!is.null(pct_mito_col)) obs_list[["pct_mito"]] <- pct_mito_col
    if (!is.null(pct_ribo_col)) obs_list[["pct_ribo"]] <- pct_ribo_col
    
    obs_df <- as.data.frame(obs_list)
    
    safe_median <- function(col) median(col, na.rm = TRUE)
    
    stats_list <- list(
      "Total Cells" = if (length(obs_list) > 0) length(obs_list[[1]]) else NA,
      "Median Genes/Cell" = if ("n_genes" %in% names(obs_list)) safe_median(obs_list[["n_genes"]]) else NA,
      # "Median UMIs/Cell" = if ("total_counts" %in% names(obs_list)) safe_median(obs_list[["total_counts"]]) else NA,
      "Median % Mitochondrial Genes" = if ("pct_mito" %in% names(obs_list)) safe_median(obs_list[["pct_mito"]]) else NA,
      "Median % Ribosomal Genes" = if ("pct_ribo" %in% names(obs_list)) safe_median(obs_list[["pct_ribo"]]) else NA
    )
    
    stats_df <- data.frame(
      Metric = names(stats_list),
      Value = unlist(stats_list),
      stringsAsFactors = FALSE
    )
    
    stats_df <- stats_df[!is.na(stats_df$Value), ]
    stats_df
  })

  output$filter_preview <- renderText({
    req(values$lazy_data)
    
    n_genes_col <- get_first_available_column(c("n_genes", "nGenes", "genes_count"))
    total_counts_col <- get_first_available_column(c("total_counts", "totalUMIs"))
    pct_mito_col <- get_first_available_column(c("pct_mito", "percent_mito", "mitochondrial_percent"))
    
    obs_list <- list()
    if (!is.null(n_genes_col)) obs_list[["n_genes"]] <- n_genes_col
    if (!is.null(total_counts_col)) obs_list[["total_counts"]] <- total_counts_col
    if (!is.null(pct_mito_col)) obs_list[["pct_mito"]] <- pct_mito_col
    
    obs_df <- as.data.frame(obs_list)
    
    filters <- c()
    passing <- rep(TRUE, nrow(obs_df))
    
    if ("n_genes" %in% names(obs_df)) {
      filters <- c(filters, "n_genes > 200")
      passing <- passing & (obs_df$n_genes > 200)
    }
    
    if ("pct_mito" %in% names(obs_df)) {
      filters <- c(filters, "pct_mito < 5")
      passing <- passing & (obs_df$pct_mito < 5)
    }
    
    if ("total_counts" %in% names(obs_df)) {
      filters <- c(filters, "total_counts > 1000")
      passing <- passing & (obs_df$total_counts > 1000)
    }
    
    passing_count <- sum(passing)
    total_count <- nrow(obs_df)
    
    paste0(passing_count, " / ", total_count, " cells (",
          round(passing_count / total_count * 100, 2), "%) pass filters: ",
          paste(filters, collapse = ", "))
  })

  output$qc_detailed_table <- DT::renderDataTable({
    req(values$lazy_data)
    
    # Get all obs column names
    all_cols <- values$lazy_data$obs_keys
    
    # Load each column lazily
    obs_list <- lapply(all_cols, function(colname) {
      values$lazy_data$get_obs_column(colname)
    })
    names(obs_list) <- all_cols
    
    # Convert to data.frame
    obs_df <- as.data.frame(obs_list)
    
    # Optionally, you can format or subset obs_df here before showing
    
    DT::datatable(
      obs_df,
      options = list(pageLength = 10, scrollX = TRUE, lengthChange = TRUE),
      filter = 'top',
      rownames = FALSE
    )
  })

}

# shinyApp(ui, server)
# shinyApp(ui = ui, server = server)