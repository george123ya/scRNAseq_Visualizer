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
  library(ps)
  library(promises)
  library(future)
  library(ggplot2)
  library(later)
  library(xml2)
  library(peakRAM)
})

# Source global configurations and helper functions
source("global.R")
source("modules/utils.R")

# Configure Python environment for scanpy/anndata

use_condaenv("sc_rna_env_python2", required = TRUE)

reticulate::py_run_string("
import zarr
import numpy as np

def get_gene_fast(z_csc, gene_index):

  n_genes = z_csc['X']['indptr'].shape[0] - 1
  n_cells = z_csc['obs'][list(z_csc['obs'].keys())[0]].shape[0]
    
  # Get start and end positions for this gene's data in CSC format
  start_idx = z_csc['X']['indptr'][gene_index]
  end_idx = z_csc['X']['indptr'][gene_index + 1]

  # If gene has no expression values
  if start_idx == end_idx:
      return np.zeros(n_cells, dtype=np.float32)

  # Get the cell indices and expression values for this gene
  cell_indices = z_csc['X']['indices'][start_idx:end_idx]
  values = z_csc['X']['data'][start_idx:end_idx]

  # Create full expression vector
  expression = np.zeros(n_cells, dtype=np.float32)
  expression[cell_indices] = values
    
  return expression
")


# use_condaenv("shiny_app_env", conda = "/opt/conda/bin/conda", required = TRUE)
# use_condaenv("shiny_app_env", required = TRUE)

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

#' Initialize fast loader for h5ad/zarr files
#' 
#' @param file_path Path to the .zarr or .h5ad file (local or URL)
#' @return List containing metadata, coordinates, and fast access functions
#' @description Efficiently loads data using the best method based on file type
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

#' Fast Zarr file loader
#' 
#' @param file_path Path to .zarr file
#' @return Processed data structure
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

#' Process Zarr data with direct access
#' 
#' @param z Zarr group object
#' @param file_path Original file path
#' @return Processed data structure with fast access functions
process_zarr_data_fast <- function(z, file_path) {
  pd <- import("pandas")
  np <- import("numpy")
  
  # Explore Zarr structure first
  cat("🔍 Exploring Zarr structure...\n")
  zarr_keys <- py_to_r(reticulate::import_builtins()$list(z$keys()))
  cat("📋 Top-level groups:", paste(zarr_keys, collapse = ", "), "\n")
  
  # Handle different Zarr structures
  if ("obs" %in% zarr_keys && "var" %in% zarr_keys) {
    # Standard AnnData Zarr structure
    obs_group <- z[["obs"]]
    var_group <- z[["var"]]
    X_array <- z[["X"]]
  } else {
    # Look for alternative structures
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
  # print(obs_columns)
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
      # If it has 'codes' and 'categories', treat as categorical
      zarr <- import("zarr")
      # Assign to Python environment
      py$col_obj <- col_obj
      py$zarr <- zarr
      
      # Check if it's a group (categorical)
      is_group <- py_eval("isinstance(col_obj, zarr.Group)")

      if (is_group) {
        codes <- to_r(col_obj[['codes']][])
        categories <- to_r(col_obj[['categories']][])
        return(categories[codes + 1])  # R is 1-indexed
      } else {
        # Numeric / array column
        return(to_r(col_obj[]))
      }
    }, error = function(e) {
      cat("⚠️ Could not load obs column", name, ":", e$message, "\n")
      return(NULL)
    })
  }


  get_gene_expression_fast <- function(gene_name, layer = "X") {
    gene_idx <- match(gene_name, gene_info$display_names)
    if (is.na(gene_idx)) {
      cat("❌ Gene", gene_name, "not found\n")
      return(NULL)
    }
    
    tryCatch({
      if (layer == "X") {
        # Ultra-fast lazy extraction via Dask-style block processing
        # expr_vec <- reticulate::py$get_gene_fast(z, as.integer(gene_idx - 1))
        peak <- peakRAM({
          expr_vec <- reticulate::py$get_gene_fast(z, as.integer(gene_idx - 1))
        })
        print(peak)
        return(py_to_r(expr_vec))
      } else {
        if (!(layer %in% layer_keys)) {
          cat("❌ Layer", layer, "not available\n")
          return(NULL)
        }
        expr_vec <- py_to_r(z$layers[[layer]][, as.integer(gene_idx - 1)])
        return(expr_vec)
      }
      
    }, error = function(e) {
      cat("⚠️ Error loading expression for", gene_name, ":", e$message, "\n")
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
    obs = NULL,  # Don't preload obs data
    var = gene_info$var_df,
    genes = gene_info$display_names,
    gene_info = gene_info,
    zarr_obj = z,  # Keep Zarr object for direct access
    seacell_df = seacell_df,
    obs_keys = obs_columns,
    available_layers = layer_keys,
    obsm_keys = obsm_keys,
    file_path = file_path,
    is_zarr = TRUE,
    # Fast access functions
    get_obs_column = get_obs_column_fast,
    get_gene_expression = get_gene_expression_fast,
    get_obsm = get_obsm_fast
  )
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


#' Fast gene information loading
load_gene_info_fast <- function(var_group) {
  tryCatch({
    # Try different ways to get var index/gene names
    var_attrs <- py_to_r(reticulate::import_builtins()$list(var_group$array_keys()))
    cat("🔍 Available var attributes:", paste(var_attrs, collapse = ", "), "\n")
    
    # Try multiple possible index/name keys
    index_keys <- c("_index", "index", "gene_ids", "var_names")
    symbol_keys <- c("gene_symbols", "gene_symbol", "Symbol", "symbol", 'gene')

    # symbol_keys <- c("gene_symbols", "gene_name", "symbol", "Symbol", 
    #                       "gene_symbol", "Gene_Symbol", "SYMBOL", "name", 
    #                       "Name", "feature_name")
    
    var_index <- NULL
    for (key in index_keys) {
      if (key %in% var_attrs) {
        tryCatch({
          var_index <- py_to_r(var_group[[key]][])
          cat("✅ Found var index using key:", key, "\n")
          break
        }, error = function(e) {
          next
        })
      }
    }
    
    # If no index found, create sequential names
    if (is.null(var_index)) {
      # Try to get var shape from any available array
      var_shape <- NULL
      for (attr in var_attrs) {
        tryCatch({
          var_shape <- var_group[[attr]]$shape[[1]]
          break
        }, error = function(e) {
          next
        })
      }
      
      if (is.null(var_shape)) {
        stop("Could not determine number of genes")
      }
      
      var_index <- paste0("GENE_", seq_len(var_shape))
      cat("⚠️ No gene index found, created sequential names\n")
    }
    
    # Try to get gene symbols
    display_names <- var_index  # default to using index
    for (key in symbol_keys) {
      if (key %in% var_attrs) {
        tryCatch({
          gene_symbols <- py_to_r(var_group[[key]][])
          # Use symbols where available, fall back to index
          display_names <- ifelse(is.na(gene_symbols) | gene_symbols == "" | gene_symbols == "nan", 
                                var_index, gene_symbols)
          cat("✅ Found gene symbols using key:", key, "\n")
          break
        }, error = function(e) {
          next
        })
      }
    }
    
    var_df <- data.frame(
      gene_id = var_index,
      gene_symbol = display_names,
      stringsAsFactors = FALSE
    )
    
    # Add other var columns if they exist
    for (attr in var_attrs) {
      if (!(attr %in% c(index_keys, symbol_keys))) {
        tryCatch({
          var_df[[attr]] <- py_to_r(var_group[[attr]][])
        }, error = function(e) {
          cat("⚠️ Could not load var attribute:", attr, "-", e$message, "\n")
        })
      }
    }
    
    return(list(
      display_names = display_names,
      var_df = var_df,
      gene_ids = var_index
    ))
    
  }, error = function(e) {
    cat("⚠️ Gene info loading failed, using fallback:", e$message, "\n")
    # Create minimal fallback
    return(list(
      display_names = paste0("GENE_", 1:1000),  # placeholder
      var_df = data.frame(gene_id = paste0("GENE_", 1:1000), stringsAsFactors = FALSE),
      gene_ids = paste0("GENE_", 1:1000)
    ))
  })
}

#' Get obs columns quickly
get_obs_columns_fast <- function(obs_group) {
  tryCatch({
    # materialize the keys generator
    py_list <- reticulate::import_builtins()$list
    obs_keys <- py_to_r(py_list(obs_group$keys()))
    
    cat("🔍 Available obs columns:", paste(head(obs_keys, 10), collapse = ", "), 
        if(length(obs_keys) > 10) "..." else "", "\n")
    return(obs_keys)
  }, error = function(e) {
    cat("⚠️ Could not get obs columns:", e$message, "\n")
    return(character(0))
  })
}


#' Get obsm keys quickly  
get_obsm_keys_fast <- function(z, zarr_keys) {
  tryCatch({
    if ("obsm" %in% zarr_keys) {
      obsm_group <- z[["obsm"]]
      return(py_to_r(reticulate::import_builtins()$list(obsm_group$array_keys())))
    } else {
      return(character(0))
    }
  }, error = function(e) {
    cat("⚠️ Could not get obsm keys:", e$message, "\n")
    return(character(0))
  })
}

#' Get var keys quickly
get_var_keys_fast <- function(z, zarr_keys) {
  tryCatch({
    if ("var" %in% zarr_keys) {
      var_group <- z[["var"]]
      return(py_to_r(reticulate::import_builtins()$list(var_group$keys())))
    } else {
      return(character(0))
    }
  }, error = function(e) {
    cat("⚠️ Could not get var keys:", e$message, "\n")
    return(character(0))
  })
}

#' Get layer keys quickly
get_layer_keys_fast <- function(z, zarr_keys) {
  tryCatch({
    if ("layers" %in% zarr_keys) {
      layers_group <- z[["layers"]]
      return(py_to_r(reticulate::import_builtins()$list(layers_group$group_keys())))
    } else {
      return(character(0))
    }
  }, error = function(e) {
    cat("⚠️ Could not get layer keys:", e$message, "\n")
    return(character(0))
  })
}

#' Fast SEACells loading
load_seacells_fast <- function(z, zarr_keys) {
  tryCatch({
    if ("uns" %in% zarr_keys) {
      uns_group <- z[["uns"]]
      uns_subkeys <- py_to_r(reticulate::import_builtins()$list(uns_group$group_keys()))
      
      if ("SEACells_summary" %in% uns_subkeys) {
        seacell_summary <- uns_group[["SEACells_summary"]]
        seacell_umap <- py_to_r(seacell_summary[["obsm"]][["X_umap"]][])
        
        seacell_df <- as.data.frame(seacell_umap)
        colnames(seacell_df) <- c("UMAP_1", "UMAP_2")
        seacell_df$SEACell <- py_to_r(seacell_summary[["obs"]][["_index"]][])
        seacell_df$cell_type <- py_to_r(seacell_summary[["obs"]][["cluster_feature"]][])
        
        cat("🔬 Loaded SEACells metacells:", nrow(seacell_df), "\n")
        return(seacell_df)
      }
    }
    return(NULL)
  }, error = function(e) {
    cat("⚠️ Could not load SEACells data:", e$message, "\n")
    return(NULL)
  })
}

# Usage example:
# data <- init_fast_zarr_h5ad("path/to/file.zarr")
# data <- init_fast_zarr_h5ad("path/to/file.h5ad")  
# data <- init_fast_zarr_h5ad("https://example.com/file.zarr")
# cell_type <- data$get_obs_column("cell_type")
# gene_expr <- data$get_gene_expression("CD34")
# umap_magic <- data$get_obsm("X_umap_magic")

#' Fallback function if init_lazy_h5ad_enhanced is not available
#' This is a simplified version for H5AD files
init_lazy_h5ad_enhanced_fallback <- function(file_path) {
  cat("Using fallback H5AD loader for:", file_path, "\n")
  
  tryCatch({
    ad <- import("anndata", convert = FALSE)
    adata <- ad$read_h5ad(file_path)
    
    # Basic processing similar to your original function
    # Extract UMAP coordinates
    umap <- py_to_r(adata$obsm[["X_umap"]])
    df <- as.data.frame(umap)
    colnames(df) <- c("x", "y")
    df$magic_x <- NA
    df$magic_y <- NA
    
    # Extract gene info
    var_df <- py_to_r(adata$var)
    genes <- rownames(var_df)
    
    # Get obs column names
    obs_cols <- py_to_r(adata$obs$columns$tolist())
    
    get_obs_column <- function(name) {
      if (!(name %in% obs_cols)) return(NULL)
      return(py_to_r(adata$obs[[name]]))
    }
    
    return(list(
      df = df,
      obs = NULL,
      var = var_df,
      genes = genes,
      gene_info = list(display_names = genes, var_df = var_df),
      ad_obj = adata,
      seacell_df = NULL,
      obs_keys = obs_cols,
      available_layers = c("X"),
      obsm_keys = py_to_r(adata$obsm_keys()),
      file_path = file_path,
      is_zarr = FALSE,
      get_obs_column = get_obs_column
    ))
    
  }, error = function(e) {
    stop("❌ Could not load H5AD file: ", e$message)
  })
}

# Check if the enhanced function exists, if not define fallback
if (!exists("init_lazy_h5ad_enhanced")) {
  init_lazy_h5ad_enhanced <- init_lazy_h5ad_enhanced_fallback
}

#' Initialize lazy-loaded h5ad file with HDF5 support
#' 
#' @param file_path Path to the .h5ad file
#' @return List containing metadata, coordinates, and AnnData object reference
#' @description Efficiently loads metadata while keeping expression data accessible
init_lazy_h5ad_enhanced <- function(file_path) {
  cat("Initializing lazy-loaded file:", file_path, "\n")
  
  start_time <- Sys.time()  # start timing
  
  tryCatch({
    anndata_exp <- import("anndata.experimental", delay_load = TRUE)
    ad_lazy <- anndata_exp$read_lazy(file_path)
    
    end_time <- Sys.time()  # end timing
    cat("⏱ Experimental lazy loading took:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
    
    return(process_lazy_anndata_enhanced(ad_lazy, file_path, is_zarr = TRUE))
  }, error = function(e) {
    cat("⚠️  Experimental lazy loading failed:", e$message, "\n")
  })
}

#' Process AnnData object for lazy-style access
#' 
#' @param ad_obj AnnData object (lazy or regular)
#' @param file_path Original file path
#' @param is_zarr Whether this is a Zarr-backed lazy object
#' @return Processed lazy data structure
process_lazy_anndata_enhanced <- function(ad_obj, file_path, is_zarr = NULL) {
  pd <- import("pandas")
  
  # Auto-detect Zarr from filename if not provided
  if (is.null(is_zarr)) {
    is_zarr <- grepl("\\.zarr$", file_path, ignore.case = TRUE)
  }
  
  # Get obs shape without loading full obs
  shape <- ad_obj$obs$shape
  n_cells <- shape[[1]]
  n_cols <- shape[[2]]
  # cat("📋 Observation metadata shape:", n_cells, "cells,", n_cols, "columns\n")
  
  # Get obs column names lazily
  obs_cols <- py_to_r(ad_obj$obs$columns$to_list())
  # cat("🔑 Available obs columns:", paste(obs_cols, collapse = ", "), "\n")
  
  # Helper to get single obs column safely as R vector
  get_obs_column <- function(name) {
    if (!(name %in% obs_cols)) return(NULL)
    col_py <- ad_obj$obs[[name]]
    # Convert lazy-backed col to pandas Series (materialize)
    col_pd <- pd$Series(col_py)
    return(py_to_r(col_pd))
  }
  
  # Load var info and process gene names
  var_df <- py_to_r(ad_obj$var)
  gene_info <- process_gene_names(ad_obj)
  
  # cat("🧬 Found", length(gene_info$display_names), "genes\n")
  
  # Get obsm keys
  obsm_keys <- py_to_r(ad_obj$obsm_keys())
  # cat("🗝️ Available obsm keys:", paste(obsm_keys, collapse = ", "), "\n")
  
  # Load UMAP coordinates
  umap_keys <- c("X_umap_normal", "X_umap", "X_umap_magic")
  umap <- NULL
  umap_key_used <- NULL

  # Time the UMAP coordinate loading
  umap_start <- Sys.time()
  for (k in umap_keys) {
    if (k %in% obsm_keys) {
      umap_py <- ad_obj$obsm[[k]]
      if (reticulate::py_has_attr(umap_py, "compute")) {
        cat("  🔄 Computing UMAP coordinates from lazy array...\n")
        umap_py <- umap_py$compute()
      }
      umap <- py_to_r(umap_py)
      umap_key_used <- k
      break
    }
  }
  umap_end <- Sys.time()
  cat("  ⏱ UMAP coordinate loading took:", round(difftime(umap_end, umap_start, units = "secs"), 3), "seconds\n")

  if (is.null(umap)) stop("❌ No UMAP coordinates found.")
  
  df <- as.data.frame(umap)
  colnames(df) <- c("x", "y")
  
  # Extract MAGIC UMAP coordinates if available
  if ("X_umap_magic" %in% obsm_keys) {
    if (umap_key_used != "X_umap_magic") {
      umap_magic <- py_to_r(ad_obj$obsm[["X_umap_magic"]])
      df$magic_x <- umap_magic[, 1]
      df$magic_y <- umap_magic[, 2]
      # cat("✨ Loaded MAGIC UMAP coordinates\n")
    } else {
      df$magic_x <- df$x
      df$magic_y <- df$y
      # cat("✨ Using main UMAP as MAGIC coordinates\n")
    }
  } else {
    df$magic_x <- NA
    df$magic_y <- NA
    # cat("⚠️ MAGIC UMAP coordinates not found\n")
  }
  
  # Expression layers
  layers_py <- ad_obj$layers$keys()
  available_layers <- reticulate::py_to_r(reticulate::import_builtins()$list(layers_py))
  # cat("📋 Available layers:", paste(available_layers, collapse = ", "), "\n")
  
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
  
  list(
    df = df,
    obs = NULL,          # Don't keep full obs loaded (too large, use get_obs_column)
    var = var_df,
    genes = gene_info$display_names,         # What user sees in UI
    gene_info = gene_info,                   # Complete gene mapping info
    ad_obj = ad_obj,
    seacell_df = seacell_df,
    obs_keys = obs_cols,
    available_layers = available_layers,
    obsm_keys = obsm_keys,
    file_path = file_path,
    is_zarr = is_zarr,
    get_obs_column = get_obs_column  # export this function for user to get columns lazily
  )
}

process_gene_names <- function(ad_obj) {
  pd <- reticulate::import("pandas")
  
  # Get var column names
  var_cols <- py_to_r(ad_obj$var$columns$to_list())

  # print(var_cols)
  
  # Helper to safely fetch a var column
  get_var_column <- function(name) {
    if (!(name %in% var_cols)) {
      return(NULL)
    }
    col_py <- ad_obj$var[[name]]
    col_pd <- pd$Series(col_py)  # materialize lazy-backed column
    return(as.character(py_to_r(col_pd)))
  }
  
  # Materialize var_names
  var_names <- as.character(py_to_r(ad_obj$var_names$to_list()))
  
  # List of possible gene symbol columns
  possible_gene_cols <- c("gene_symbols", "gene_name", "symbol", "Symbol", 
                          "gene_symbol", "Gene_Symbol", "SYMBOL", "name", 
                          "Name", "feature_name")
  
  # Find the first matching column
  gene_symbols <- NULL
  gene_symbols_col <- NULL
  for (col in possible_gene_cols) {
    col_data <- get_var_column(col)
    if (!is.null(col_data)) {
      gene_symbols <- col_data
      gene_symbols_col <- col
      # cat("📋 Found gene symbols in column:", col, "\n")
      break
    }
  }
  
  # Determine if var_names contains Ensembl IDs
  has_ensembl_ids <- any(grepl("^ENS[A-Z]*[0-9]+", var_names, ignore.case = TRUE))
  
  # Decide display names and mapping
  if (has_ensembl_ids && !is.null(gene_symbols)) {
    display_names <- gene_symbols
    missing <- is.na(display_names) | display_names == "" | display_names == "NA"
    display_names[missing] <- var_names[missing]
    
    name_to_index <- setNames(var_names, display_names)
    
    # Handle duplicated symbols
    dup <- duplicated(display_names) | duplicated(display_names, fromLast = TRUE)
    if (any(dup)) {
      display_names[dup] <- paste0(gene_symbols[dup], " (", var_names[dup], ")")
      name_to_index <- setNames(var_names, display_names)
      cat("⚠️ Found", sum(dup), "duplicated gene symbols, added ENS IDs as suffix\n")
    }
    
    # cat("✅ Using gene symbols for display, ENS IDs for indexing\n")
    
  } else if (!has_ensembl_ids && !is.null(gene_symbols)) {
    display_names <- gene_symbols
    missing <- is.na(display_names) | display_names == "" | display_names == "NA"
    display_names[missing] <- var_names[missing]
    
    name_to_index <- setNames(var_names, display_names)
    # cat("✅ Using cleaned gene symbols for display\n")
    
  } else {
    display_names <- var_names
    name_to_index <- setNames(var_names, display_names)
    # cat("✅ Using var_names directly for display and indexing\n")
  }
  
  # cat("📊 Total genes:", length(display_names), "\n")
  # cat("📊 Unique display names:", length(unique(display_names)), "\n")
  
  return(list(
    display_names = display_names,
    var_names = var_names,
    name_to_index = name_to_index,
    has_ensembl_ids = has_ensembl_ids,
    gene_symbols_col = gene_symbols_col
  ))
}



#' Extract gene expression data with smart caching
#' 
#' @param gene_names Character vector of gene names to extract
#' @param lazy_data List containing processed AnnData data
#' @param use_layer Character, which expression layer to use
#' @param use_cache Whether to use expression cache
#' @return List with encoded expression data and statistics
extract_gene_data_enhanced <- function(gene_names, lazy_data, use_layer = "X", use_cache = TRUE) {
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
    # Map display names to var_names (columns in Zarr)
    var_names_for_extraction <- character(0)
    valid_display_names <- character(0)
    for (gene_name in gene_names) {
      idx <- which(lazy_data$gene_info$display_names == gene_name)
      if (length(idx) > 0) {
        var_names_for_extraction <- c(var_names_for_extraction, lazy_data$gene_info$gene_ids[idx])
        valid_display_names <- c(valid_display_names, gene_name)
      } else {
        cat("⚠️ Gene not found:", gene_name, "\n")
      }
    }

    if (length(var_names_for_extraction) == 0) {
      cat("⚠️ No valid genes found\n")
      return(list(
        genes = character(0),
        data = character(0),
        magic_data = character(0),
        ranges = list(),
        magic_ranges = list()
      ))
    }

    cat("🔍 Extracting", length(var_names_for_extraction), "genes:", 
        paste(head(valid_display_names, 5), collapse = ", "), 
        if(length(valid_display_names) > 5) "..." else "", "\n")

    # Extract expression data via helper (handles Zarr CSR)
    expr_list <- lapply(valid_display_names, function(g) {
      lazy_data$get_gene_expression(gene_name = g, layer = use_layer)
    })
    expr_matrix <- do.call(cbind, expr_list)

    # Ensure correct orientation
    if (nrow(expr_matrix) != nrow(lazy_data$df) && ncol(expr_matrix) == nrow(lazy_data$df)) {
      expr_matrix <- t(expr_matrix)
      cat("🔄 Transposed expression matrix to correct orientation\n")
    }

    cat("📊 Final expression matrix:", nrow(expr_matrix), "cells ×", ncol(expr_matrix), "genes\n")
    gene_ranges <- calculate_gene_ranges(expr_matrix, valid_display_names)

    # MAGIC layer
    magic_data <- ""
    magic_ranges <- list()
    if ("MAGIC_imputed_data" %in% lazy_data$available_layers) {
      tryCatch({
        cat("✨ Extracting MAGIC imputed data...\n")
        magic_matrix <- lazy_data$get_gene_expression(var_names = var_names_for_extraction, layer = "MAGIC_imputed_data")
        if (nrow(magic_matrix) != nrow(lazy_data$df) && ncol(magic_matrix) == nrow(lazy_data$df)) {
          magic_matrix <- t(magic_matrix)
        }
        magic_data_vector <- as.vector(magic_matrix)
        magic_data <- base64enc::base64encode(writeBin(as.numeric(magic_data_vector), raw()))
        magic_ranges <- calculate_gene_ranges(magic_matrix, valid_display_names)
        cat("✨ MAGIC data extracted successfully\n")
      }, error = function(e) {
        cat("⚠️ Could not extract MAGIC data:", e$message, "\n")
      })
    }

    # Encode main expression matrix
    gene_data_vector <- as.vector(expr_matrix)
    encoded_data <- base64enc::base64encode(writeBin(as.numeric(gene_data_vector), raw()))
    cat("💾 Encoded data size:", nchar(encoded_data), "characters\n")
    cat("💾 MAGIC data size:", nchar(magic_data), "characters\n")

    return(list(
      genes = valid_display_names,
      data = encoded_data,
      magic_data = magic_data,
      ranges = gene_ranges,
      magic_ranges = magic_ranges,
      nrows = nrow(expr_matrix),
      ncols = ncol(expr_matrix),
      layer_used = use_layer,
      var_names_used = var_names_for_extraction
    ))

  }, error = function(e) {
    cat("❌ Error in enhanced gene extraction:", e$message, "\n")
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
    is_zarr = is_zarr,
    file_type = if(is_zarr) "Zarr-backed" else "HDF5"
  )
  
  # cat("✅ Dataset info generated successfully\n")
  return(info)
}

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

# Helper: get first existing obs column from a list of possible names
get_first_available_column <- function(lazy_data, possible_names) {
  for (nm in possible_names) {
    col_data <- lazy_data$get_obs_column(lazy_data$z, nm)
    if (!is.null(col_data)) {
      return(col_data)
    }
  }
  return(NULL)
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

  # --- Dataset Selection ---

  # UPDATED: Enhanced dataset loading with improved lazy loading
  observeEvent(input$dataset, {
    req(input$dataset)
    
    # Skip if placeholder values
    if (input$dataset %in% c("", "❌ Folder not found", "❌ No .h5ad files found")) {
      return()
    }

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
      values$lazy_data <- init_fast_zarr_h5ad(input$dataset)
      
      # Determine success message based on strategy
      strategy_msg <- if(values$lazy_data$is_zarr) {
        "✅ True lazy loading active (Zarr-backed)"
      } else {
        "✅ Smart selective loading active (HDF5)"
      }
      
      showNotification(strategy_msg, type = "message", duration = 4)
      
      # cat("✅ Initialized smart loading for:", basename(input$dataset), "\n")
      # cat("   Path:", input$dataset, "\n")
      # cat("   Cells:", nrow(values$lazy_data$df), "\n")
      # cat("   Genes:", length(values$lazy_data$genes), "\n")
      # cat("   Strategy:", values$lazy_data$memory_strategy, "\n")
      # cat("   File type:", if(values$lazy_data$is_zarr) "Zarr" else "HDF5", "\n")
      
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
    
    selectInput("dataset", "Select Dataset:", 
                choices = choices,
                selected = choices[[1]][1],
                width = "100%")
  })

  # --- Lazy Data Loading ---

  # Process data when lazy_data is available
  observeEvent(values$lazy_data, {
    req(values$lazy_data)
    
    tryCatch({
      session$sendCustomMessage("showSpinner", TRUE)
      output$status <- renderText("Processing metadata...")
      
      t_start <- Sys.time()
      cat("⏱ Starting lazy_data processing at", format(t_start), "\n")
      
      # Use the loaded lazy_data
      raw_data <- values$lazy_data$df
      
      # --- Process dataframe ---
      t1 <- Sys.time()
      processed_result <- process_dataframe(raw_data)
      values$processed_data <- processed_result
      t2 <- Sys.time()
      cat("  ⏱ process_dataframe took", round(difftime(t2, t1, units = "secs"), 3), "seconds\n")
      
      # --- Create color by choices ---
      annotation_names <- names(processed_result$annotations)
      values$annotation_names <- annotation_names
      
      # --- Get dataset info ---
      t3 <- Sys.time()
      values$dataset_info <- tryCatch({
        NULL
        # get_dataset_info_smart(values$lazy_data)
      }, error = function(e) {
        cat("⚠️ Error getting dataset info:", e$message, "\n")
        list(
          n_cells = nrow(values$lazy_data$df),
          n_genes = length(values$lazy_data$genes),
          available_layers = character(0),
          obsm_keys = character(0),
          has_magic_coords = FALSE,
          has_magic_expr = FALSE,
          has_seacells = FALSE,
          file_size_mb = 0,
          is_zarr = FALSE,
          file_type = "Unknown"
        )
      })
      t4 <- Sys.time()
      cat("  ⏱ get_dataset_info_smart took", round(difftime(t4, t3, units = "secs"), 3), "seconds\n")
      
      # --- Encode data for JavaScript ---
      t5 <- Sys.time()
      encoded_data <- encode_data(processed_result)
      n_cols <- ncol(processed_result$data)
      t6 <- Sys.time()
      cat("  ⏱ encode_data took", round(difftime(t6, t5, units = "secs"), 3), "seconds\n")
      
      # --- Send to JavaScript ---
      # Send to JavaScript with precise timing
      t7 <- Sys.time()
      
      # Generate unique timing ID
      timing_id <- paste0("update_", format(Sys.time(), "%Y%m%d_%H%M%S_"), 
                        sprintf("%03d", as.integer(runif(1) * 1000)))
      
      # Convert R time to JavaScript milliseconds (Unix timestamp * 1000)
      send_time_js <- as.numeric(Sys.time()) * 1000
      
      cat("📤 Sending updateData at:", format(Sys.time()), "\n")
      cat("📤 Timing ID:", timing_id, "\n")
      cat("📤 Send time (JS format):", send_time_js, "\n")
      
      session$sendCustomMessage("updateData", list(
        base64 = encoded_data,
        annotationData = processed_result$annotations,
        numCols = n_cols,
        clusters = processed_result$annotations$cluster$names,
        colors = processed_result$cluster_colors,
        metacellColors = processed_result$metacell_colors,
        geneExprRanges = processed_result$gene_expr_ranges,
        sendTime = send_time_js,      # JavaScript-compatible timestamp
        timingId = timing_id          # Unique ID for tracking
      ))
      
      t8 <- Sys.time()
      cat("⏱ sendCustomMessage preparation took:", 
          round(difftime(t8, t7, units = "secs"), 3), "seconds\n")
      
      # Store timing info for callback verification
      timing_info <- list(
        r_send_time = Sys.time(),
        js_send_time = send_time_js,
        timing_id = timing_id
      )
      
      assign(paste0("timing_", timing_id), timing_info, envir = .GlobalEnv)
      cat("  ⏱ sendCustomMessage(updateData) took", round(difftime(t8, t7, units = "secs"), 3), "seconds\n")
      
      session$sendCustomMessage("showSpinner", FALSE)

      output$status <- renderText(paste0("Rendered ", nrow(raw_data), " cells with ", 
                                        length(annotation_names), " annotations."))
      
      t_end <- Sys.time()
      cat("⏱ Total observer time:", round(difftime(t_end, t_start, units = "secs"), 3), "seconds\n")
      
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

  # Enhanced callback with timing ID verification
  observeEvent(input$updateData_done, {
    callback_receive_time <- Sys.time()
    
    # Handle both old and new formats
    if (is.numeric(input$updateData_done)) {
      # Old simple format - just log the timestamp
      cat("⏱ Simple callback received:", input$updateData_done, "\n")
      return()
    }
    
    if (is.list(input$updateData_done)) {
      timing_data <- input$updateData_done
      timing_id <- timing_data$timingId
      
      if (!is.null(timing_id)) {
        # Retrieve stored timing info
        timing_var <- paste0("timing_", timing_id)
        
        if (exists(timing_var, envir = .GlobalEnv)) {
          stored_timing <- get(timing_var, envir = .GlobalEnv)
          r_send_time <- stored_timing$r_send_time
          total_r_to_r_time <- difftime(callback_receive_time, r_send_time, units = "secs")
          
          # Determine data source
          is_cloud <- grepl("^https://", values$current_dataset %||% "")
          data_source <- if(is_cloud) "CLOUD" else "LOCAL"
          
          cat("\n🎯 Timing Breakdown (", data_source, ") for", timing_id, ":\n")
          cat("  📤 R send time:", format(r_send_time), "\n")
          cat("  📥 R callback time:", format(callback_receive_time), "\n")
          
          if (!is.null(timing_data$transferTime)) {
            cat("  ⏱ R→JS transfer:", round(timing_data$transferTime, 2), "ms\n")
          }
          
          if (!is.null(timing_data$jsProcessingTime)) {
            cat("  ⏱ JavaScript TOTAL:", round(timing_data$jsProcessingTime, 2), "ms\n")
          }
          
          if (!is.null(timing_data$redrawTime)) {
            cat("  ⏱ redrawAllPlots:", round(timing_data$redrawTime, 2), "ms\n")
          }
          
          if (!is.null(timing_data$totalEndToEndTime)) {
            cat("  ⏱ End-to-end (JS):", round(timing_data$totalEndToEndTime, 2), "ms\n")
          }
          
          cat("  ⏱ R-to-R roundtrip:", round(total_r_to_r_time * 1000, 2), "ms\n")
          
          if (!is.null(timing_data$dataStats)) {
            stats <- timing_data$dataStats
            cat("  📊 Data:", stats$pointCount, "points,", stats$numCols, "cols\n")
          }
          
          if (!is.null(timing_data$error)) {
            cat("  ❌ Error:", timing_data$error, "\n")
          }
          
          cat("  ================================\n\n")
          
          # Clean up
          rm(list = timing_var, envir = .GlobalEnv)
          
        } else {
          cat("⚠️ Could not find stored timing for ID:", timing_id, "\n")
        }
      } else {
        cat("⚠️ Timing data missing timingId\n")
        print(str(timing_data))
      }
    } else {
      cat("⚠️ Invalid updateData_done format:", class(input$updateData_done), "\n")
      print(input$updateData_done)
    }
  })

  
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

  # --- Gene Search ---
  
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
          placeholder = paste0("Type gene names ..."),
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
      gene_data <- extract_gene_data_enhanced(gene_list, values$lazy_data)
      
      session$sendCustomMessage("showSpinner", FALSE)
      
      # print(paste("📊 Smart extracted data summary:"))
      # print(paste("  - Genes found:", length(gene_data$genes)))
      # print(paste("  - Normal data size:", nchar(gene_data$data)))
      # print(paste("  - MAGIC data size:", nchar(gene_data$magic_data)))
      # print(paste("  - Layer used:", gene_data$layer_used))
      # print(paste("  - Strategy:", gene_data$memory_strategy))
      
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
        output$status <- renderText(paste0("Ready - ", values$dataset_info$n_cells, " cells, ", 
                                          values$dataset_info$n_genes, " genes"))
      }
    }
  }, ignoreNULL = FALSE)

  
  # --- Annotation Color ---

  # Color by UI
  output$colorByUI <- renderUI({
    if (is.null(categorical_vars())) {
      selectInput("colorBy", "Color cells by:", choices = list("Initialize dataset first..." = "loading"))
    } else {
      # choices <- setNames(values$annotation_names, tools::toTitleCase(gsub("_", " ", values$annotation_names)))
      choices <- c(categorical_vars(), numeric_obs_keys())
      selectInput("colorBy", "Color cells by:", choices = choices, selected = choices[1])
    }
  })
  
  # Handle color by changes
  observeEvent(input$colorBy, {
    if (!is.null(input$colorBy) && input$colorBy != "loading") {
      
      # Get annotation column
      annotation_data <- values$lazy_data$get_obs_column(values$lazy_data$z, input$colorBy)
      
      ann_colors <- generate_colors(length(unique(annotation_data)))
      var_type <- ifelse(is.numeric(annotation_data), "gene", "categorical")
      
      if (var_type != 'gene') {
        numeric_annotation_data <- as.numeric(factor(annotation_data, exclude = NULL))
      } else {
        numeric_annotation_data <- as.numeric(annotation_data)
      }
      
      # Convert to binary and encode as base64
      binary_data <- writeBin(as.numeric(numeric_annotation_data), 
                            raw(), 
                            size = 4, 
                            endian = "little")
      annotation_raw <- base64enc::base64encode(binary_data)
      
      session$sendCustomMessage("colorByChange", list(
        colorBy = input$colorBy,
        annotation_raw = annotation_raw,  # Send binary data
        annotation_length = length(numeric_annotation_data),  # Include length for validation
        names = sort(unique(annotation_data)),
        colors = ann_colors,
        var_type = var_type
      ))
    }
  })


  # --- MAGIC Visual ---

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

  observeEvent(input$activateMAGIC, {
    session$sendCustomMessage("toggleMAGIC", input$activateMAGIC)
  })


  # --- Data Information Display ---

  # UPDATED: Enhanced data info output
  output$dataInfo <- renderUI({
    req(values$lazy_data, values$processed_data, values$dataset_info)
    
    info <- values$dataset_info
    result <- values$processed_data
    
    # Create strategy-specific messaging
    strategy_color <- "green"
    strategy_icon <- "⚡"
    strategy_text <- "True Lazy Loading"
    
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

  # --- Selection Point Output ---

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

  # ---- Precompute column statistics once ----
  obs_stats <- reactive({
    req(values$lazy_data)
    keys <- values$lazy_data$obs_keys
    keys <- keys[!keys %in% c("_index", "_categories")]
    # Load each obs column once
    cols <- lapply(keys, function(k) {
      values$lazy_data$get_obs_column(values$lazy_data$z, k)
    })

    names(cols) <- keys
    
    # Compute stats for each column
    lapply(cols, function(col) {
      nuniq <- length(unique(col))
      
      # Detect integer-like numeric (robust)
      is_numeric <- is.numeric(col)
      is_integer_like <- is_numeric && all(!is.na(col)) &&
        all(abs(col - round(col)) < .Machine$double.eps^0.5)
      is_real_numeric <- is_numeric && !is_integer_like
      
      is_categorical <- (is.factor(col) || is.character(col) || is_integer_like) &&
                        nuniq <= 20
      
      list(
        categorical = is_categorical,
        numeric     = is_real_numeric && nuniq > 1
      )
    })
  })

  # ---- Categorical obs keys ----
  categorical_vars <- reactive({
    stats <- obs_stats()
    sort(names(Filter(function(s) s$categorical, stats)))
  })

  # ---- Numeric obs keys ----
  numeric_obs_keys <- reactive({
    stats <- obs_stats()
    sort(names(Filter(function(s) s$numeric, stats)))
  })


  # --- QC Tab ---

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
      x <- as.numeric(values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_x_metric))
      y <- as.numeric(values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_y_metric))
      
      if (!is.null(input$qc_color_by) && input$qc_color_by != "None") {
        color <- as.factor(values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_color_by))
        data.frame(X = x, Y = y, Color = color)
      } else {
        data.frame(X = x, Y = y)
      }
      
    } else {
      req(input$qc_metric, input$qc_group_by)
      metric <- as.numeric(values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_metric))
      group <- as.factor(values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_group_by))
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


  output$qc_stats_table <- renderTable({
    req(values$lazy_data)

    n_genes_col <- get_first_available_column(values$lazy_data, c("n_genes", "nGenes", "genes_count"))
    total_counts_col <- get_first_available_column(values$lazy_data, c("total_counts", "totalUMIs"))
    pct_mito_col <- get_first_available_column(values$lazy_data, c("pct_mito", "percent_mito", "mitochondrial_percent"))
    pct_ribo_col <- get_first_available_column(values$lazy_data, c("pct_ribo", "percent_ribo", "ribosomal_percent"))

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

    n_genes_col <- get_first_available_column(values$lazy_data, c("n_genes", "nGenes", "genes_count"))
    total_counts_col <- get_first_available_column(values$lazy_data, c("total_counts", "totalUMIs"))
    pct_mito_col <- get_first_available_column(values$lazy_data, c("pct_mito", "percent_mito", "mitochondrial_percent"))

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
      values$lazy_data$get_obs_column(values$lazy_data$z, colname)
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