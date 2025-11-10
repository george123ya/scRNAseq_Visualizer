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

# Load required libraries
suppressPackageStartupMessages({
  library(shiny)
  # library(shinydashboard)
  # library(DT)
  library(shinyjs)
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
  library(reshape2)
  library(rlang)
  library(promises)
  library(redux)
  library(future)
  library(future.redis)
  library(arrow)
})

# future::plan(multisession)
# Source global configurations and helper functions
source("global.R")
source("modules/utils.R")



# Configure Python environment for scanpy/anndata

# use_condaenv("sc_rna_env_python2", required = TRUE)
source_python("./scripts/rank_genes_zarr.py")  # Load ONCE at startup
use_condaenv("shiny_app_env", conda = "/opt/conda/bin/conda", required = TRUE)
# use_condaenv("shiny_app_env", required = TRUE)

# Enable direct S3 reads (no fsspec, no Python)
Sys.setenv(
  ARROW_S3_ALLOW_UNSAFE_URLS = "true",
  AWS_EC2_METADATA_DISABLED = "true"  # speeds up S3 client init
)

# Optional: If you use public buckets (no credentials)
s3_bucket_opts <- S3FileSystem$create(
  anonymous = TRUE
)

reticulate::py_run_string("
import zarr
import numpy as np
from scipy.sparse import csc_matrix, hstack

def get_gene_fast(z_csc, gene_index, layer='X'):
    # Get number of cells
    n_cells = z_csc['obs'][list(z_csc['obs'].array_keys())[0]].shape[0]
    
    # Select correct CSC structure
    if layer == 'X':
        csc = z_csc['X']
    else:
        csc = z_csc['layers'][layer]

    # Get start and end positions for this gene in the selected layer
    start_idx = csc['indptr'][gene_index]
    end_idx = csc['indptr'][gene_index + 1]

    # If gene has no expression values
    if start_idx == end_idx:
        return np.zeros(n_cells, dtype=np.float32)

    # Get the cell indices and expression values for this gene
    cell_indices = csc['indices'][start_idx:end_idx]
    values = csc['data'][start_idx:end_idx]

    # Create full expression vector
    expression = np.zeros(n_cells, dtype=np.float32)
    expression[cell_indices] = values
    
    return expression

# Optimized get_genes_batch for non-contiguous indices
def get_genes_batch(z_csc, gene_indices, layer='X', batch_size=1000):
    n_cells = z_csc['obs'][next(iter(z_csc['obs'].array_keys()))].shape[0]
    csc = z_csc['X'] if layer == 'X' else z_csc['layers'][layer]
    
    # Cache indptr to reduce cloud fetches
    indptr = csc['indptr'][:]
    indices = csc['indices']
    data = csc['data']
    
    # Sort gene indices
    sorted_indices = sorted(enumerate(gene_indices), key=lambda x: x[1])
    sorted_genes = [x[1] for x in sorted_indices]
    original_order = [x[0] for x in sorted_indices]
    
    if not sorted_genes:
        return csc_matrix((n_cells, 0), dtype=np.float32)
    
    # Initialize lists for non-zero entries
    all_row_indices = []
    all_col_indices = []
    all_data_values = []
    
    # Process genes in batches
    for batch_start in range(0, len(sorted_genes), batch_size):
        batch_end = min(batch_start + batch_size, len(sorted_genes))
        batch_genes = sorted_genes[batch_start:batch_end]
        
        if not batch_genes:
            continue
        
        # Get the range of data we need to fetch
        min_gene = batch_genes[0]
        max_gene = batch_genes[-1]
        
        # Fetch indptr for the full range
        range_indptr = indptr[min_gene:max_gene+2]
        start_idx = range_indptr[0]
        end_idx = range_indptr[-1]
        
        if start_idx < end_idx:
            # Fetch the contiguous block of data
            batch_indices = indices[start_idx:end_idx]
            batch_data = data[start_idx:end_idx]
        else:
            batch_indices = np.array([], dtype=np.int32)
            batch_data = np.array([], dtype=np.float32)
        
        # Process each gene in this batch
        for local_i, gene_idx in enumerate(batch_genes):
            # Calculate relative position within the fetched range
            rel_pos = gene_idx - min_gene
            s_rel = range_indptr[rel_pos] - start_idx
            e_rel = range_indptr[rel_pos + 1] - start_idx
            
            if s_rel < e_rel:
                # Global column index (position in the sorted list)
                global_col_idx = batch_start + local_i
                
                # Extract data for this specific gene
                gene_row_indices = batch_indices[s_rel:e_rel]
                gene_data = batch_data[s_rel:e_rel]
                
                # Add to global lists
                all_row_indices.append(gene_row_indices)
                all_col_indices.append(np.full(len(gene_row_indices), global_col_idx, dtype=np.int32))
                all_data_values.append(gene_data)
    
    # Combine all non-zero entries
    if not all_row_indices:
        return csc_matrix((n_cells, len(gene_indices)), dtype=np.float32)
    
    all_row_indices = np.concatenate(all_row_indices)
    all_col_indices = np.concatenate(all_col_indices)
    all_data_values = np.concatenate(all_data_values)
    
    # Create final sparse matrix
    result_mat = csc_matrix(
        (all_data_values, (all_row_indices, all_col_indices)),
        shape=(n_cells, len(gene_indices)),
        dtype=np.float32
    )
    
    # Reorder columns to match original gene_indices order
    inverse_order = np.argsort(original_order)
    result_mat = result_mat[:, inverse_order]
    
    return result_mat.toarray()
")

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
  
  # Create fast access functions
  get_obs_column_fast <- function(z, name, cells = NULL) {
    to_r <- reticulate::py_to_r
    builtins <- reticulate::import_builtins()
    
    if (!inherits(z, "python.builtin.object")) stop("⚠️ z is not a Python object.")
    
    root_keys <- to_r(builtins$list(z$keys()))
    if (!("obs" %in% root_keys)) return(NULL)
    
    obs_keys <- to_r(builtins$list(z[["obs"]]$keys()))
    if (!(name %in% obs_keys)) return(NULL)
    
    col_obj <- z[["obs"]][[name]]
    
    tryCatch({
      zarr <- reticulate::import("zarr")
      py$col_obj <- col_obj
      py$zarr <- zarr
      
      is_group <- reticulate::py_eval("isinstance(col_obj, zarr.Group)")
      full_col <- NULL
      
      if (is_group) {
        full_codes <- to_r(col_obj[['codes']][])
        full_categories <- to_r(col_obj[['categories']][])
        
        # infer whether categories are numeric or character
        if (is.numeric(full_categories) && length(unique(full_categories)) < 30) {
          # treat as numeric factor
          full_col <- vapply(
            full_codes,
            function(code) if (code < 0) NA_real_ else full_categories[code + 1],
            FUN.VALUE = numeric(1)
          )
        } else {
          # treat as character factor
          full_col <- vapply(
            full_codes,
            function(code) if (code < 0) NA_character_ else as.character(full_categories[code + 1]),
            FUN.VALUE = character(1)
          )
        }
      } else {
        full_col <- to_r(col_obj[])
      }
      
      # handle cell subsetting
      if (is.null(cells) || identical(cells, "all")) {
        return(full_col)
      } else if (is.numeric(cells) && all(cells >= 1) && all(cells <= length(full_col))) {
        return(full_col[cells])
      } else {
        cat("⚠️ Invalid cells argument of length", length(cells), "\n")
        valid_cells <- intersect(cells, seq_along(full_col))
        return(full_col[valid_cells])
      }
    }, error = function(e) {
      cat("⚠️ Error loading obs column", name, ":", e$message, "\n")
      return(NULL)
    })
  }

  
  get_gene_expression_fast <- function(gene_name, layer = "X", cells = NULL) {
      gene_idx <- match(gene_name, gene_info$display_names)
      cat("🔍 Looking for gene:", gene_name, "Found at index:", gene_idx, "\n")
      if (is.na(gene_idx)) {
        cat("❌ Gene", gene_name, "not found\n")
        return(NULL)
      }
      
      tryCatch({
        gene_idx_py <- as.integer(gene_idx - 1)
        full_expr_vec_py <- NULL
        if (layer == "X") {
          peak <- peakRAM({
            full_expr_vec_py <- reticulate::py$get_gene_fast(z, gene_idx_py, 'X')
          })
          print(peak)
          full_expr_vec <- py_to_r(full_expr_vec_py)
        } else {
          if (!(layer %in% layer_keys)) {
            cat("❌ Layer", layer, "not available\n")
            return(NULL)
          }
          full_expr_vec_py <- reticulate::py$get_gene_fast(z, gene_idx_py, layer)
          full_expr_vec <- py_to_r(full_expr_vec_py)
        }
        
        # Handle cell subsetting
        if (is.null(cells) || identical(cells, "all")) {
          return(full_expr_vec)
        } else if (is.numeric(cells) && all(cells >= 1) && all(cells <= length(full_expr_vec))) {
          return(full_expr_vec[cells])
        } else {
          cat("2")
          cat("⚠️ Invalid cells argument: must be NULL, 'all', or a valid list of 1-based indices\n")
          return(full_expr_vec)  # Fallback to full
        }
      }, error = function(e) {
        cat("⚠️ Error loading expression for", gene_name, ":", e$message, "\n")
        return(NULL)
      })
  }
  
  get_genes_expression_batch <- function(gene_names, layer = "X", cells = NULL) {
      
      gene_indices <- match(gene_names, toupper(gene_info$display_names))
      valid_mask <- !is.na(gene_indices)
      valid_indices <- as.integer(gene_indices[valid_mask] - 1)  # Ensure integer indices
      valid_genes <- gene_names[valid_mask]
      
      if (length(valid_indices) == 0) {
        cat("❌ No valid genes found in batch\n")
        return(NULL)
      }
      
      tryCatch({
        peak <- peakRAM({
          full_expr_mat_py <- reticulate::py$get_genes_batch(z, valid_indices, layer)
        })
        print(peak)
        full_expr_mat <- py_to_r(full_expr_mat_py)
        colnames(full_expr_mat) <- valid_genes
        cat("📏 Batch retrieved expression for", length(valid_genes), "genes\n")
        
        # Handle cell subsetting (rows)
        if (is.null(cells) || identical(cells, "all")) {
          return(full_expr_mat)
        } else if (is.numeric(cells) && all(cells >= 1) && all(cells <= nrow(full_expr_mat))) {
          return(full_expr_mat[cells, ])
        } else {
          cat("3")
          cat("⚠️ Invalid cells argument: must be NULL, 'all', or a valid list of 1-based indices\n")
          return(full_expr_mat)  # Fallback to full
        }
      }, error = function(e) {
        cat("⚠️ Error loading batch expression:", e$message, "\n")
        cat("Run `reticulate::py_last_error()` for details.\n")
        return(NULL)
      })
  }
  
  get_obsm_fast <- function(key) {
    if (!(key %in% obsm_keys)) return(NULL)
    tryCatch({
      return(py_to_r(z[['obsm']][[key]][]))
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
    n_cells = n_cells,
    n_genes = n_genes
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
      return(py_to_r(reticulate::import_builtins()$list(layers_group$keys())))
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
  cat("📋 Observation metadata shape:", n_cells, "cells,", n_cols, "columns\n")
  
  # Get obs column names lazily
  obs_cols <- py_to_r(ad_obj$obs$columns$to_list())
  cat("🔑 Available obs columns:", paste(obs_cols, collapse = ", "), "\n")
  
  # Helper to get single obs column safely as R vector
  get_obs_column <- function(z, name) {
    if (!(name %in% obs_cols)) return(NULL)
    col_py <- ad_obj$obs[[name]]
    col_pd <- pd$Series(col_py)
    return(py_to_r(col_pd))
  }
  
  # Load var info and process gene names
  var_df <- py_to_r(ad_obj$var)
  gene_info <- process_gene_names(ad_obj)
  cat("🧬 Found", length(gene_info$display_names), "genes\n")
  
  # Get obsm keys
  obsm_keys <- py_to_r(ad_obj$obsm_keys())
  cat("🗝️ Available obsm keys:", paste(obsm_keys, collapse = ", "), "\n")
  
  # Load UMAP coordinates
  umap_keys <- c("X_umap_normal", "X_umap", "X_umap_magic")
  umap <- NULL
  umap_key_used <- NULL
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
  
  list(
    df = df,
    obs = NULL,
    var = var_df,
    genes = gene_info$display_names,
    gene_info = gene_info,
    ad_obj = ad_obj,
    seacell_df = seacell_df,
    obs_keys = obs_cols,
    available_layers = available_layers,
    obsm_keys = obsm_keys,
    file_path = file_path,
    is_zarr = is_zarr,
    get_obs_column = get_obs_column,
    # pathways = pathways  # Add pathways to the returned list
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
extract_gene_data_enhanced <- function(gene_names, lazy_data, use_layer = "X", use_cache = TRUE,
                                       only_vector = FALSE) {
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
        expr_list <- lapply(var_names_for_extraction, function(g) {
          lazy_data$get_gene_expression(gene_name = g, layer = "MAGIC_imputed_data")
        })
        magic_matrix <- do.call(cbind, expr_list)
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

    if (only_vector) {
      return(gene_data_vector)
    }

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

  # Safely get file or dataset size in megabytes for local paths or various cloud URLs
  get_dataset_size_mb <- function(path_or_url) {
    tryCatch({
      if (grepl("^https?://", path_or_url)) {
        if (grepl("storage\\.googleapis\\.com", path_or_url)) {
          base_url <- sub("/[^/]*$", "/", path_or_url)
          resource_name <- sub(".*/", "", path_or_url)
          get_cloud_size_mb(base_url, resource_name)
        } else if (grepl("\\.s3\\.[a-z0-9-]+\\.amazonaws\\.com", path_or_url)) {
          base_url <- sub("/[^/]*$", "/", path_or_url)
          resource_name <- sub(".*/", "", path_or_url)
          get_cloud_size_mb(base_url, resource_name)
        } else if (grepl("\\.blob\\.core\\.windows\\.net", path_or_url)) {
          base_url <- sub("/[^/]*$", "/", path_or_url)
          resource_name <- sub(".*/", "", path_or_url)
          get_cloud_size_mb(base_url, resource_name)
        } else {
          get_generic_url_size_mb(path_or_url)
        }
      } else {
        get_local_size_mb(path_or_url)
      }
    }, error = function(e) {
      cat("⚠️ Error getting size:", e$message, "\n")
      0
    })
  }

  # Placeholder for generic URL size fetching
  get_generic_url_size_mb <- function(url) {
    require(httr)
    response <- HEAD(url)
    if (status_code(response) == 200) {
      size_bytes <- as.numeric(headers(response)$`content-length`)
      if (is.na(size_bytes)) {
        stop("Content-Length header not available")
      }
      size_mb <- size_bytes / (1024 * 1024)  # Convert bytes to MB
      return(size_mb)
    } else {
      stop("Failed to retrieve size from URL")
    }
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
  tryCatch({
    zarr_files <- c()
    marker <- NULL
    
    # Normalize base_url to ensure it ends with "/"
    if (!grepl("/$", base_url)) base_url <- paste0(base_url, "/")
    
    # Determine provider and set query parameters
    if (grepl("storage\\.googleapis\\.com", base_url)) {
      # Google Cloud Storage (restore original logic)
      repeat {
        url <- if (!is.null(marker)) {
          paste0(base_url, "?marker=", URLencode(marker, reserved = TRUE))
        } else {
          base_url
        }
        
        doc <- read_xml(url)  # Use original read_xml for GCS
        ns <- xml_ns(doc)
        keys <- xml_text(xml_find_all(doc, ".//d1:Key", ns))
        
        # Keep only base Zarr names
        zarr_base <- sub("^(.+?\\.zarr).*", "\\1", keys)
        zarr_base <- zarr_base[grepl("\\.zarr$", zarr_base)]
        
        zarr_files <- unique(c(zarr_files, zarr_base))
        
        truncated <- xml_text(xml_find_first(doc, ".//d1:IsTruncated", ns)) == "true"
        if (!truncated) break
        
        marker <- xml_text(xml_find_first(doc, ".//d1:NextMarker", ns))
        if (is.na(marker) || marker == "") break
      }
      
    } else if (grepl("\\.s3\\.[a-z0-9-]+\\.amazonaws\\.com", base_url)) {
      # AWS S3
      repeat {
        query <- list(
          "list-type" = "2",
          delimiter = "/"
        )
        if (!is.null(marker)) {
          query$`continuation-token` <- marker
        }
        
        response <- GET(base_url, query = query)
        if (status_code(response) != 200) {
          stop(sprintf("Failed to list S3 bucket: %s", status_code(response)))
        }
        
        doc <- read_xml(content(response, as = "text"))
        ns <- xml_ns(doc)
        keys <- xml_text(xml_find_all(doc, ".//d1:Key", ns))
        
        # Keep only base Zarr names
        zarr_base <- sub("^(.+?\\.zarr).*", "\\1", keys)
        zarr_base <- zarr_base[grepl("\\.zarr$", zarr_base)]
        
        zarr_files <- unique(c(zarr_files, zarr_base))
        
        truncated <- xml_text(xml_find_first(doc, ".//d1:IsTruncated", ns)) == "true"
        if (!truncated) break
        
        marker <- xml_text(xml_find_first(doc, ".//d1:NextContinuationToken", ns))
        if (is.na(marker) || marker == "") break
      }
      
    } else if (grepl("\\.blob\\.core\\.windows\\.net", base_url)) {
      # Azure Blob Storage
      repeat {
        query <- list(
          restype = "container",
          comp = "list",
          delimiter = "/"
        )
        if (!is.null(marker)) {
          query$marker <- marker
        }
        
        response <- GET(base_url, query = query)
        if (status_code(response) != 200) {
          stop(sprintf("Failed to list Azure container: %s", status_code(response)))
        }
        
        doc <- read_xml(content(response, as = "text"))
        ns <- xml_ns(doc)
        keys <- xml_text(xml_find_all(doc, ".//d1:Name", ns))
        
        # Keep only base Zarr names
        zarr_base <- sub("^(.+?\\.zarr).*", "\\1", keys)
        zarr_base <- zarr_base[grepl("\\.zarr$", zarr_base)]
        
        zarr_files <- unique(c(zarr_files, zarr_base))
        
        marker <- xml_text(xml_find_first(doc, ".//d1:NextMarker", ns))
        truncated <- !is.na(marker) && marker != ""
        if (!truncated) break
      }
      
    } else {
      stop("Unsupported cloud provider")
    }
    
    return(zarr_files)
    
  }, error = function(e) {
    cat("⚠️ Error listing Zarr datasets:", e$message, "\n")
    character()
  })
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

show_prefix_modal <- function() {
  modalDialog(
    title = "Version prefix (optional)",
    textInput("version_prefix_input",
              label = "Prefix (e.g. my_analysis)",
              placeholder = "Leave blank for auto-date"),
    footer = tagList(
      modalButton("Cancel"),
      actionButton("prefix_ok", "OK", class = "btn-primary")
    ),
    easyClose = FALSE,
    fade = TRUE
  )
}


# Fixed server function with proper initialization

server <- function(input, output, session) {

  volcano_selected_indices <- reactiveVal(integer(0))
  volcano_updating_from_plot <- reactiveVal(FALSE)
  volcano_updating_from_table <- reactiveVal(FALSE)
  volcano_table_data <- reactiveVal(NULL)

  # Helper function for preparing the table data (outside the server function)
  prepare_volcano_table <- function(df, selected_indices, mode, max_rows = 1000) {
    cat("\n=== prepare_volcano_table DEBUG ===\n")
    cat("Mode:", mode, "\n")
    cat("Selected indices:", paste(head(selected_indices, 10), collapse = ", "), "\n")
    
    # Add original row index for mapping selection back to full data
    df$original_index <- seq_len(nrow(df))
    
    # Add significance score for sorting
    if (!"significance_score" %in% colnames(df)) {
      df$significance_score <- df$neg_log10_padj * abs(df$logfoldchanges)
    }
    
    # 1. Filter and sort
    if (mode == "selected") {
      if (length(selected_indices) > 0) {
        display_df <- df[selected_indices, ]
      } else {
        display_df <- df[FALSE, ]
      }
      display_df <- display_df[order(-display_df$neg_log10_padj), ]
      
    } else { # mode == "all"
      df_sorted <- df[order(-df$significance_score), ]
      is_selected <- df_sorted$original_index %in% selected_indices
      selected_genes <- df_sorted[is_selected, ]
      unselected_genes <- df_sorted[!is_selected, ]
      
      n_selected <- nrow(selected_genes)
      n_top_unselected <- max_rows - n_selected
      
      if (n_top_unselected > 0) {
        top_unselected <- head(unselected_genes, n_top_unselected)
        display_df <- rbind(selected_genes, top_unselected)
      } else {
        display_df <- selected_genes
      }
      
      display_df$is_selected <- display_df$original_index %in% selected_indices
      display_df <- display_df[order(-display_df$is_selected, -display_df$neg_log10_padj), ]
    }
    
    cat("Display_df rows:", nrow(display_df), "\n")
    cat("First gene:", display_df$gene[1], "\n")
    cat("First log2FC:", display_df$logfoldchanges[1], "\n")
    cat("First selected status:", display_df$original_index[1] %in% selected_indices, "\n")
    
    # 2. Build display data frame step by step
    gene_col <- as.character(display_df$gene)
    sel_col <- display_df$original_index %in% selected_indices
    log2fc_col <- round(display_df$logfoldchanges, 3)
    pval_col <- formatC(display_df$pvals_adj, format = "e", digits = 2)
    neglog10_col <- round(display_df$neg_log10_padj, 2)
    status_col <- as.character(display_df$significant)
    
    cat("Column lengths - Gene:", length(gene_col), "Sel:", length(sel_col), 
        "log2FC:", length(log2fc_col), "\n")
    cat("First 3 Sel values:", sel_col[1:3], "\n")
    cat("First 3 log2FC values:", log2fc_col[1:3], "\n")
    
    # Create data frame explicitly
    display_cols <- data.frame(
      Gene = gene_col,
      Sel = sel_col,
      log2FC = log2fc_col,
      pvalue = pval_col,
      neglog10p = neglog10_col,
      Status = status_col,
      stringsAsFactors = FALSE
    )
    
    cat("Data frame created, dimensions:", nrow(display_cols), "x", ncol(display_cols), "\n")
    cat("Column names:", paste(names(display_cols), collapse = ", "), "\n")
    cat("First row:\n")
    print(display_cols[1, ])
    
    # Rename columns for display
    names(display_cols) <- c("Gene", "Sel", "log2FC", "p-value", "-log10(p)", "Status")
    
    cat("After rename, column names:", paste(names(display_cols), collapse = ", "), "\n")
    cat("First row after rename:\n")
    print(display_cols[1, ])
    cat("=== END DEBUG ===\n\n")
    
    # Return display data
    return(list(
      display = display_cols,
      original_indices = display_df$original_index,
      selected_rows = which(sel_col == TRUE)
    ))
  }

  # Reactive values to store processed data
  values <- reactiveValues(
    processed_data = NULL,
    annotation_names = NULL,
    lazy_data = NULL,
    current_dataset = NULL,
    dataset_info = NULL
  )

  results <- reactiveVal(NULL)
  status <- reactiveVal("No analysis run yet.")

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

  # UI: Dataset selection
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

    # Load cloud configuration
    load_cloud_config <- function(config_file = "cloud_config.json") {
      if (file.exists(config_file)) {
        config <- fromJSON(config_file)
        return(config)
      } else {
        warning("Config file not found, using default empty config")
        return(list(folders = character(), urls = character()))
      }
    }

    config <- load_cloud_config("cloud_config.json")
    
    # Cloud folders and Zarr datasets
    cloud_choices <- NULL
    if (length(config$folders) > 0) {
      for (folder_url in config$folders) {
        cloud_files <- get_top_level_zarr(folder_url)
        if (length(cloud_files) > 0) {
          provider <- if (grepl("storage\\.googleapis\\.com", folder_url)) "GCS"
                    else if (grepl("\\.s3\\.", folder_url)) "S3"
                    else if (grepl("\\.blob\\.core\\.windows\\.net", folder_url)) "Azure"
                    else "Cloud"
          cloud_choices <- c(cloud_choices, setNames(
            paste0(folder_url, cloud_files),
            paste0("[", provider, "] ", cloud_files)
          ))
        }
      }
    }
    
    # Predefined cloud URLs from config
    if (length(config$urls) > 0) {
      cloud_choices <- c(cloud_choices, setNames(
        config$urls,
        paste0("[Predefined] ", basename(config$urls))
      ))
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

  # Store fetched versions here
  versions_list <- reactiveVal(list(Original = "original"))

  # R function to fetch versions from the proxy server (No changes needed here, 
  # as it handles HTTP errors via stop_for_status, which throws to the catch block)
  fetch_versions_async <- function(zarr_url) {
      LIST_URL <- "http://localhost:8080/list_versions"
      
      response <- httr::GET(
          url = LIST_URL,
          query = list(zarr_url = zarr_url)
      )
      
      # This throws an error if status is 4xx/5xx, immediately jumping to the catch block
      httr::stop_for_status(response, task = paste("fetch versions for", zarr_url))
      
      return(httr::content(response, as = "parsed"))
  }
  
  # R function to robustly convert S3 HTTPS URLs to s3:// format
  normalize_url <- function(url) {
      # 1. Check if it's already s3://
      if (grepl("^s3://", url)) {
          return(url)
      }
      
      # 2. Handle GCS URLs (often work as is)
      if (grepl("^https://storage.googleapis.com/", url)) {
          return(url)
      }
      
      # 3. Handle AWS HTTPS URLs (The main fix for the observed error)
      if (grepl("\\.amazonaws\\.com/", url)) {
          # Pattern 1: Virtual hosted style (bucket.s3.region.amazonaws.com)
          # Capture Group 1: Bucket Name (scrnaseq-browser)
          # Capture Group 2: The path/key (/strict_epdsc_annotated_data_csc_full.zarr)
          
          # We need to extract the part before ".s3." (the bucket) and the path after ".com/"
          match <- regexec("https?://([^.]+)\\.s3\\.[^/]+\\.amazonaws\\.com/(.*)", url)
          
          if (match[[1]][1] != -1) {
              # Extract the captured groups (requires substring and attr(match[[1]], "match.length"))
              matches <- regmatches(url, match)
              
              if (length(matches[[1]]) >= 3) {
                  bucket_name <- matches[[1]][2]
                  object_key <- matches[[1]][3]
                  return(paste0("s3://", bucket_name, "/", object_key))
              }
          }
          
          # Fallback for path-style or other complex formats: simply replace the domain base
          # This is a bit risky, but better than the error you currently have.
          s3_url <- sub("https?://[^/]*s3[^/]*\\.amazonaws\\.com/", "s3://", url)
          return(s3_url)
      }

      # 4. Default return for other types (e.g., local files, which should be filtered elsewhere)
      return(url) 
  }

  version_refresh_trigger <- reactiveVal(0)

  observe({
      req(input$dataset)
      version_refresh_trigger()

      default_choices <- list(Original = "original")
      
      is_zarr <- grepl("\\.zarr$", input$dataset)
      is_remote <- grepl("^s3://|^https?://", input$dataset) 

      if (!is_zarr || !is_remote) {
          if (is_zarr) {
              shiny::showNotification(
                  "Local Zarr file selected. Versions are only tracked for cloud datasets.", 
                  duration = 5, 
                  type = "default"
              )
          }
          versions_list(default_choices)
          return()
      }
      
      zarr_url <- input$dataset
      zarr_url_normalized <- normalize_url(zarr_url)
      
      print(paste("Fetching versions for Zarr URL:", zarr_url_normalized))

      # FIX: Specify exactly what to export and load
      future({
          # Define the function inline to avoid environment capture
          LIST_URL <- "http://localhost:8080/list_versions"
          
          response <- httr::GET(
              url = LIST_URL,
              query = list(zarr_url = zarr_url_normalized)
          )
          
          httr::stop_for_status(response, task = paste("fetch versions for", zarr_url_normalized))
          
          return(httr::content(response, as = "parsed"))
          
      }, globals = list(zarr_url_normalized = zarr_url_normalized),  # Only export the URL
        packages = c("httr")  # Explicitly specify required packages
      ) %>%
      promises::then(function(result) {
          if (length(result$versions) > 0) {
              version_choices <- sapply(result$versions, function(v) {
                  display_name <- paste0(
                      v$version_prefix, 
                      " (", v$metadata$analysis_date, 
                      ") - ", v$metadata$user_id
                  )
                  return(display_name)
              }, USE.NAMES = FALSE)
              
              version_values <- sapply(result$versions, function(v) v$version_prefix)
              choices <- c(default_choices, setNames(version_values, version_choices))
              versions_list(choices)
              
          } else {
              versions_list(default_choices)
          }
      }) %>% 
      promises::catch(function(error) {
          error_msg <- sub(".*:\\s*", "", error$message) 
          
          shiny::showNotification(
              paste("Failed to load versions (", error_msg, "). Defaulting to 'Original'."),
              duration = 8,
              type = "warning"
          )
          
          versions_list(default_choices)
      })
  })

  output$versionUI <- renderUI({
      req(input$dataset)
      if (!grepl("\\.zarr$", input$dataset)) {
          return(NULL)
      }
      
      choices <- versions_list()
      
      selectInput(
          "version_prefix",
          "Select Version:", 
          choices = choices,
          selected = choices[1],
          width = "100%"
      )
  })

  observeEvent(input$version_prefix, {
    req(values$current_dataset)
    req(values$lazy_data)
    req(input$version_prefix)

    print(values$current_dataset)

    # Only proceed if the dataset is a Zarr file and not local
    if (
      !grepl("\\.zarr$", values$current_dataset) ||
      grepl("^(/|\\.|[a-zA-Z]:\\\\)", values$current_dataset)
    ) {
      return()
    }

    version_prefix <- input$version_prefix

    # 🆕 HELPER FUNCTIONS TO LOAD VERSION DATA
    get_version_object <- function(base_url, version_prefix, array_name) {
      fsspec <- import("fsspec")
      zarr <- import("zarr")
      np <- import("numpy", convert = FALSE)

      full_array_url <- paste0(
        base_url, "/versions/", 
        version_prefix, "/masks/", 
        array_name
      )

      array_store <- fsspec$get_mapper(full_array_url)
      z_array_py <- zarr$open_array(store = array_store, mode = "r")
      r_data <- z_array_py[]
      r_data_converted <- py_to_r(r_data)
    }

    get_version_metadata <- function(base_url, version_prefix) {
      fsspec <- import("fsspec")
      json_py <- import("json")

      zattrs_url <- paste0(
        base_url, "/versions/", 
        version_prefix, "/.zattrs"
      )

      file_obj <- fsspec$open(zattrs_url, mode = "rt")
      f <- file_obj$`__enter__`() 
      metadata_content_py <- f$read()
      file_obj$`__exit__`(NULL, NULL, NULL)

      metadata_py <- json_py$loads(metadata_content_py)
      version_metadata_r <- py_to_r(metadata_py)
    }

    if (version_prefix == "original") {
      cat("🔄 Restoring ORIGINAL version (full dataset)\n")
      
      # 🆕 CLEAR ALL FILTERS
      filtered_indices(NULL)
      
      # 🆕 GET CURRENT ANNOTATION TO RESTORE ALL VISIBLE
      current_anno <- input$colorBy
      if (!is.null(current_anno) && current_anno != "loading") {
        annotation_data <- values$lazy_data$get_obs_column(values$lazy_data$z, current_anno)
        all_categories <- seq(0, length(unique(annotation_data)) - 1)
        
        # ⭐ TRIGGER COLORBY CHANGE TO FORCE REDRAW
        annotation_numeric <- as.numeric(factor(annotation_data, exclude = NULL))
        
        # Get colors
        ann_colors <- generate_colors(length(unique(annotation_data)))
        
        # Convert to binary and encode
        binary_data <- writeBin(as.numeric(annotation_numeric), raw(), size = 4, endian = "little")
        annotation_raw <- base64enc::base64encode(binary_data)
        
        # Send to JS to redraw
        session$sendCustomMessage("colorByChange", list(
          colorBy = current_anno,
          annotation_raw = annotation_raw,
          annotation_length = length(annotation_numeric),
          names = sort(unique(annotation_data)),
          colors = ann_colors,
          var_type = "categorical",
          visibleCategories = all_categories  # ⭐ Include this
        ))
      }
      
      # Also clear numeric filter UI
      session$sendCustomMessage("clearAllFilters", list(
        timestamp = as.numeric(Sys.time()) * 1000
      ))
      
      cat("✅ Full dataset restored with", length(all_categories), "categories\n")
      return()
    }

    # 🆕 LOAD VERSION DATA (for non-original versions)
    cell_mask_version <- get_version_object(
      base_url = values$current_dataset,
      version_prefix = version_prefix,
      array_name = "cell_mask"
    )
    
    visible_categories_version <- get_version_object(
      base_url = values$current_dataset,
      version_prefix = version_prefix,
      array_name = "visible_categories"
    )
    
    selected_points_version <- get_version_object(
      base_url = values$current_dataset,
      version_prefix = version_prefix,
      array_name = "selected_points"
    )
    
    annotation_metadata <- get_version_metadata(
      base_url = values$current_dataset,
      version_prefix = version_prefix
    )

    # 🆕 DETERMINE FILTERING STRATEGY
    n_selected <- sum(selected_points_version)
    n_masked <- sum(cell_mask_version)
    n_total <- length(cell_mask_version)
    
    cat("🔄 Loading version:", version_prefix, "\n")
    cat("   Total cells:", n_total, "\n")
    cat("   Masked cells:", n_masked, "\n")
    cat("   Selected points:", n_selected, "\n")
    
    filter_to_selection <- FALSE
    
    if (n_selected > 0 && n_selected < n_masked) {
      # User made a lasso selection - show ONLY these cells
      filter_to_selection <- TRUE
      cat("   Strategy: FILTER to", n_selected, "selected cells\n")
      
      # Update filtered_indices to match selection
      selected_indices <- which(selected_points_version)
      filtered_indices(selected_indices)
      
    } else if (n_masked < n_total) {
      # Cells were filtered (category/numeric) but no lasso
      filter_to_selection <- FALSE
      cat("   Strategy: HIGHLIGHT", n_masked, "filtered cells\n")
      
      # Update filtered_indices to match mask
      masked_indices <- which(cell_mask_version)
      filtered_indices(masked_indices)
    } else {
      # No filtering
      cat("   Strategy: Show ALL cells\n")
      filtered_indices(NULL)
    }
    
    # 🆕 RESTORE NUMERIC FILTERS if they exist
    if (!is.null(annotation_metadata$numeric_filters) && 
        length(annotation_metadata$numeric_filters) > 0) {
      cat("   Restoring", length(annotation_metadata$numeric_filters), "numeric filters\n")
      
      # Send to JS to rebuild filter UI
      session$sendCustomMessage("restoreNumericFilters", list(
        filters = annotation_metadata$numeric_filters
      ))
    } else {
      # Clear numeric filters if version has none
      session$sendCustomMessage("clearAllFilters", list(
        timestamp = as.numeric(Sys.time()) * 1000
      ))
    }

    # 🆕 SEND RESTORE MESSAGE WITH FILTERING STRATEGY
    session$sendCustomMessage("restoreState", list(
      visibleCategories = visible_categories_version,
      selectedPoints = if (filter_to_selection) selected_points_version else NULL,
      annotationName = annotation_metadata$umap_colorby,
      filterToSelection = filter_to_selection,
      cellMask = cell_mask_version,
      isOriginal = FALSE
    ))
  })

  # --- Lazy Data Loading ---

  # Process data when lazy_data is available
  observeEvent(values$lazy_data, {
    req(values$lazy_data)
    
    tryCatch({
      session$sendCustomMessage("showSpinner", TRUE)
      # output$status <- renderText("Processing metadata...")
      
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

      # Set in plot storage that canvas exists
      plot_storage$canvases <- list(main = TRUE)
      
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

      # output$status <- renderText(paste0("Rendered ", nrow(raw_data), " cells with ", 
      #                                   length(annotation_names), " annotations."))
      
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

  # UI: Gene search for visualization tab
  output$geneSearchUI <- renderUI({
    message("Rendering geneSearchUI") # Debug
    if (is.null(values$lazy_data)) {
      selectizeInput("geneSearch", "🔍 Search Gene:",
                    choices = list("Initialize dataset first..." = "loading"),
                    multiple = TRUE
      )
    } else {
      selectizeInput("geneSearch", "🔍 Search Gene:",
                    choices = NULL, # Rely on server-side updates
                    multiple = TRUE,
                    options = list(
                      placeholder = "Type gene names...",
                      openOnFocus = FALSE,
                      closeAfterSelect = TRUE,
                      plugins = list("remove_button"),
                      maxItems = 15,
                      onInitialize = I("function() { console.log('geneSearch initialized'); }") # Debug
                    )
      )
    }
  })

  # UPDATED: Gene search observer with smart extraction
  observeEvent(list(input$geneSearch, matrix_for_viz()), {
    req(values$lazy_data)
    
    # Send layer change to JavaScript
    session$sendCustomMessage("updateLayer", list(
      layer = matrix_for_viz()
    ))
    
    if (!is.null(input$geneSearch) && length(input$geneSearch) > 0) {
      gene_list <- as.character(input$geneSearch)
      
      print(paste("🔍 Smart searching for genes:", paste(gene_list, collapse = ", ")))
      
      # Show loading status
      loading_msg <- paste("Loading expression data for", length(gene_list), "genes...")
      output$status <- renderText(loading_msg)
      session$sendCustomMessage("showSpinner", TRUE)
      
      # Extract gene data for the selected layer
      gene_data <- extract_gene_data_enhanced(gene_list, values$lazy_data, use_layer = matrix_for_viz())
      print(paste("Expr extracted from layer:", matrix_for_viz()))
      
      session$sendCustomMessage("showSpinner", FALSE)
      session$sendCustomMessage("geneSearchChange", list(
        genes = gene_list,
        expression_data = gene_data
      ))
      
      # Keep 'main', clear the rest
      plot_storage$canvases <- plot_storage$canvases["main"]
      
      # Update with current gene_list
      for (gene in gene_list) {
        plot_storage$canvases[[gene]] <- TRUE
      }
      
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

  # Track last dataset change time (in ms, to match JS Date.now())
  values$last_dataset_change <- reactiveVal(0)

  # Update timestamp when dataset changes
  observeEvent(input$dataset, {
    values$last_dataset_change(as.numeric(Sys.time()) * 1000)
  })

  observeEvent(input$dataset, {
    req(input$dataset)
    
    cat("\n🔄 DATASET CHANGED - Clearing all filters\n")
    
    # ⭐ Reset numeric filters
    filtered_indices(NULL)
    
    # Reset visible categories timestamp
    values$last_dataset_change(as.numeric(Sys.time()) * 1000)
    
    # Send clear message to JavaScript
    session$sendCustomMessage("clearAllFilters", list(timestamp = Sys.time()))
  })

  visible_categories <- reactive({
    vis_cats <- input$visibleCategories
    
    # Return NULL if no value or timestamp is before last dataset change
    if (is.null(vis_cats) || !is.numeric(vis_cats$timestamp) || vis_cats$timestamp < values$last_dataset_change()) {
      return(NULL)
    }
    
    vis_cats
  })

  # Reactive to store selected points
  selected_points <- reactive({
    input$selectedPoints
  })

  # Define a reactive expression for the subset indices
  subset_indices_reactive <- reactive({
    req(input$dataset, values$lazy_data)
    
    n_total <- nrow(values$lazy_data$df)
    
    cat("\n🔄 RECALCULATING SUBSET PIPELINE\n")
    cat("═════════════════════════════════════════════════════════\n")
    
    # ========================================================================
    # STAGE 0: MASTER DATASET
    # ========================================================================
    
    cat("Stage 0 - MASTER: ", format(n_total, big.mark=","), " cells\n")
    current_indices <- seq_len(n_total)
    
    # ========================================================================
    # STAGE 1: CATEGORY FILTERING
    # ========================================================================
    
    cat("Stage 1 - CATEGORY: ")
    
    vis_cats <- visible_categories()
    if (!is.null(vis_cats) && !is.null(vis_cats$annotationName)) {
      annotation_name <- vis_cats$annotationName
      visible_indices <- unlist(vis_cats$visibleIndices)
      
      # Get annotation data for ALL cells
      point_data <- values$lazy_data$get_obs_column(values$lazy_data$z, annotation_name)
      category_names <- sort(unique(point_data))
      visible_category_names <- category_names[visible_indices + 1]
      
      # Filter to cells in visible categories
      cat_mask <- point_data %in% visible_category_names
      cat_mask[is.na(cat_mask)] <- FALSE
      
      current_indices <- which(cat_mask)
      cat("Active - ", format(length(current_indices), big.mark=","), " cells\n")
    } else {
      cat("Inactive - ", format(length(current_indices), big.mark=","), " cells\n")
    }
    
    # ========================================================================
    # STAGE 2: NUMERIC FILTERING
    # ========================================================================
    
    cat("Stage 2 - NUMERIC: ")
    
    filtered <- filtered_indices()
    if (!is.null(filtered) && length(filtered) > 0) {
      # Keep only cells that pass numeric filters AND are in current_indices
      current_indices <- intersect(current_indices, filtered)
      cat("Active - ", format(length(current_indices), big.mark=","), " cells\n")
    } else {
      cat("Inactive - ", format(length(current_indices), big.mark=","), " cells\n")
    }
    
    # ========================================================================
    # STAGE 3: LASSO SELECTION (FINAL)
    # ========================================================================
    
    cat("Stage 3 - LASSO: ")
    
    sel_points <- selected_points()
    if (!is.null(sel_points) && !is.null(sel_points$selectedIndices) && 
        length(sel_points$selectedIndices) > 0) {
      selected_indices <- as.numeric(unlist(sel_points$selectedIndices)) + 1
      
      # Keep only cells that are selected AND passed all previous filters
      current_indices <- intersect(current_indices, selected_indices)
      cat("Active - ", format(length(current_indices), big.mark=","), " cells (FINAL)\n")
    } else {
      cat("Inactive - ", format(length(current_indices), big.mark=","), " cells (FINAL)\n")
    }
    
    cat("═════════════════════════════════════════════════════════\n")
    cat("✅ FINAL SUBSET:", format(length(current_indices), big.mark=","), "cells\n\n")
    
    list(
      combined = current_indices  # Use this for all visualizations
    )
  })


  filtered_indices <- reactiveVal(NULL)

  observeEvent(input$applyFilters, {
    req(values$lazy_data)
    
    raw_filters <- input$applyFilters
    
    cat("\n🔍 ========== FILTER RECEIVER ==========\n")
    cat("Raw class:", class(raw_filters), "\n")
    
    # Handle empty array (Clear button sends empty JSON [])
    if (is.null(raw_filters) || raw_filters == "" || raw_filters == "[]") {
      cat("✅ CLEAR BUTTON PRESSED - Resetting all numeric filters\n")
      cat("🔍 ==========================================\n\n")
      
      # ⭐ CRITICAL: Set to NULL to clear all numeric filtering
      filtered_indices(NULL)
      
      # Also update the plots immediately
      session$sendCustomMessage("updateFilteredIndices", list(
        filteredIndices = NULL,
        count = nrow(values$lazy_data$df),
        active = FALSE
      ))
      
      return()
    }
    
    # ⭐ PARSE JSON
    filters <- tryCatch({
      parsed <- jsonlite::fromJSON(raw_filters, simplifyDataFrame = FALSE)
      cat("✅ Parsed", length(parsed), "filter(s)\n")
      parsed
    }, error = function(e) {
      cat("❌ JSON parse error:", e$message, "\n")
      list()
    })
    
    if (length(filters) == 0) {
      cat("❌ No valid filters\n")
      cat("🔍 ==========================================\n\n")
      filtered_indices(NULL)
      return()
    }
    
    # Get current subset (from previous filters)
    current_subset <- seq_len(nrow(values$lazy_data$df))
    n_total <- nrow(values$lazy_data$df)
    
    # Start with all cells
    passing_mask <- rep(FALSE, n_total)
    passing_mask[current_subset] <- TRUE
    
    cat("📊 Starting with", sum(passing_mask), "cells\n\n")
    
    # Apply each filter with detailed logging
    for (i in seq_along(filters)) {
      f <- filters[[i]]
      
      cat("━━ Filter", i, "━━\n")
      
      # Extract parameters
      if (is.list(f)) {
        col_name <- f$col
        min_thresh <- as.numeric(f$minThresh)
        max_thresh <- as.numeric(f$maxThresh)
      } else {
        col_name <- f["col"]
        min_thresh <- as.numeric(f["minThresh"])
        max_thresh <- as.numeric(f["maxThresh"])
      }
      
      col_name <- as.character(col_name)
      
      cat("  Column: '", col_name, "'\n", sep="")
      cat("  Range: [", min_thresh, ", ", max_thresh, "]\n", sep="")
      
      # Validate
      if (is.null(col_name) || is.na(col_name) || col_name == "" || 
          is.na(min_thresh) || is.na(max_thresh)) {
        cat("  ❌ Invalid - skipping\n\n")
        next
      }
      
      # GET FULL COLUMN DATA
      col_data <- tryCatch({
        values$lazy_data$get_obs_column(values$lazy_data$z, col_name)
      }, error = function(e) {
        cat("  ❌ Error:", e$message, "\n")
        NULL
      })
      
      if (is.null(col_data) || length(col_data) != n_total) {
        cat("  ⚠️ Column invalid - skipping\n\n")
        next
      }
      
      # Apply threshold
      col_numeric <- as.numeric(col_data)
      filter_mask <- (col_numeric >= min_thresh) & (col_numeric <= max_thresh)
      filter_mask[is.na(filter_mask)] <- FALSE
      
      # Apply AND logic
      passing_mask <- passing_mask & filter_mask
      
      remaining <- sum(passing_mask)
      cat("  ✅ After filter:", remaining, "cells\n\n")
    }
    
    # Extract final indices and store
    passing_indices <- which(passing_mask)
    
    # ⭐ CRITICAL: Update filtered_indices() reactive
    filtered_indices(passing_indices)
    
    cat("═════════════════════════════════════════════\n")
    cat("✅ RESULT:", length(passing_indices), "cells pass all", length(filters), "filters\n")
    cat("🔍 ==========================================\n\n")
  })


  # Optionally store the combined in reactiveValues for broader access
  observe({
    values$global_subset_indices <- subset_indices_reactive()$combined
  })

  output$selectedInfo <- renderUI({
    subset_info <- subset_indices_reactive()
    vis_cats <- visible_categories()
    sel_points <- selected_points()
    filtered <- filtered_indices()
    
    n_master <- nrow(values$lazy_data$df)
    
    # ========================================================================
    # CALCULATE COUNTS IN PIPELINE ORDER
    # ========================================================================
    
    # Stage 0: Master dataset
    n_stage_0 <- n_master
    
    # Stage 1: After category filter
    n_stage_1 <- n_master
    cat_active <- FALSE
    cat_details <- ""
    
    if (!is.null(vis_cats) && !is.null(vis_cats$annotationName)) {
      annotation_name <- vis_cats$annotationName
      visible_indices <- unlist(vis_cats$visibleIndices)
      
      point_data <- values$lazy_data$get_obs_column(values$lazy_data$z, annotation_name)
      category_names <- sort(unique(point_data))
      visible_category_names <- category_names[visible_indices + 1]
      
      cat_mask <- point_data %in% visible_category_names
      n_stage_1 <- sum(cat_mask, na.rm = TRUE)
      cat_active <- TRUE
      
      # Build category detail string
      cat_counts <- table(point_data[cat_mask])
      cat_details_lines <- sapply(names(cat_counts), function(cat) {
        sprintf("  %s: %s", cat, format(as.numeric(cat_counts[cat]), big.mark = ","))
      })
      cat_details <- paste(cat_details_lines, collapse = "\n")
    }
    
    # Stage 2: After numeric filters
    n_stage_2 <- n_stage_1
    numeric_active <- FALSE
    
    if (!is.null(filtered) && length(filtered) > 0) {
      n_stage_2 <- length(filtered)
      numeric_active <- TRUE
    }
    
    # Stage 3: After lasso selection (FINAL)
    n_stage_3 <- n_stage_2
    lasso_active <- FALSE
    
    if (!is.null(sel_points) && !is.null(sel_points$selectedIndices) && 
        length(sel_points$selectedIndices) > 0) {
      n_stage_3 <- length(unlist(sel_points$selectedIndices))
      lasso_active <- TRUE
    }
    
    # ========================================================================
    # BUILD DISPLAY IN CORRECT ORDER
    # ========================================================================
    
    HTML(
      paste(
        # Main container
        "<div style='background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%); padding: 20px; border-radius: 10px; margin-bottom: 20px;'>",
        
        # Title
        "<div style='font-size: 16px; font-weight: bold; color: #2c3e50; margin-bottom: 15px;'>",
        "📊 CELL SCREENING PIPELINE",
        "</div>",
        
        # Stage 0: Master dataset
        "<div style='background: white; padding: 12px; border-radius: 6px; margin-bottom: 10px; border-left: 4px solid #3498db;'>",
        "<div style='font-size: 12px; color: #7f8c8d;'>🎯 MASTER DATASET</div>",
        "<div style='font-size: 24px; font-weight: bold; color: #2c3e50;'>",
        format(n_stage_0, big.mark = ","),
        " <span style='font-size: 14px; color: #7f8c8d;'>cells</span></div>",
        "</div>",
        
        # Stage 1: Category filter
        "<div style='background: ", if(cat_active) "#e8f5e9" else "#f5f5f5", "; padding: 12px; border-radius: 6px; margin-bottom: 10px; border-left: 4px solid ", if(cat_active) "#4caf50" else "#ccc", ";'>",
        "<div style='font-size: 12px; color: #7f8c8d; margin-bottom: 4px;'>",
        if(cat_active) "1️⃣ CATEGORIES (ACTIVE)" else "1️⃣ Categories",
        "</div>",
        "<div style='font-size: 20px; font-weight: bold; color: #2c3e50;'>",
        format(n_stage_1, big.mark = ","),
        " <span style='font-size: 14px; color: #7f8c8d;'>cells</span>",
        " <span style='font-size: 12px; color: ", if(cat_active) "#4caf50" else "#999", ";'>",
        if(cat_active) paste0("(", round(n_stage_1/n_stage_0*100, 1), "%)") else "—",
        "</span>",
        "</div>",
        if(nzchar(cat_details)) {
          paste0(
            "<div style='font-size: 11px; color: #666; margin-top: 8px; padding-top: 8px; border-top: 1px solid rgba(0,0,0,0.1);'>",
            gsub("\n", "<br/>", cat_details),
            "</div>"
          )
        } else {
          ""
        },
        "</div>",
        
        # Stage 2: Numeric filters
        "<div style='background: ", if(numeric_active) "#e3f2fd" else "#f5f5f5", "; padding: 12px; border-radius: 6px; margin-bottom: 10px; border-left: 4px solid ", if(numeric_active) "#2196f3" else "#ccc", ";'>",
        "<div style='font-size: 12px; color: #7f8c8d; margin-bottom: 4px;'>",
        if(numeric_active) "2️⃣ NUMERIC FILTERS (ACTIVE)" else "2️⃣ Numeric Filters",
        "</div>",
        "<div style='font-size: 20px; font-weight: bold; color: #2c3e50;'>",
        format(n_stage_2, big.mark = ","),
        " <span style='font-size: 14px; color: #7f8c8d;'>cells</span>",
        " <span style='font-size: 12px; color: ", if(numeric_active) "#2196f3" else "#999", ";'>",
        if(numeric_active) paste0("(", round(n_stage_2/n_stage_1*100, 1), "%)") else "—",
        "</span>",
        "</div>",
        "</div>",
        
        # Stage 3: Lasso selection (FINAL)
        "<div style='background: ", if(lasso_active) "#fff3e0" else "#f5f5f5", "; padding: 12px; border-radius: 6px; margin-bottom: 10px; border-left: 4px solid ", if(lasso_active) "#ff9800" else "#ccc", ";'>",
        "<div style='font-size: 12px; color: #7f8c8d; margin-bottom: 4px;'>",
        if(lasso_active) "3️⃣ LASSO SELECTION (ACTIVE - FINAL)" else "3️⃣ Lasso Selection",
        "</div>",
        "<div style='font-size: 20px; font-weight: bold; color: #2c3e50;'>",
        format(n_stage_3, big.mark = ","),
        " <span style='font-size: 14px; color: #7f8c8d;'>cells</span>",
        " <span style='font-size: 12px; color: ", if(lasso_active) "#ff9800" else "#999", ";'>",
        if(lasso_active) paste0("(", round(n_stage_3/n_stage_2*100, 1), "%)") else "—",
        "</span>",
        "</div>",
        "</div>",
        
        # Final summary bar
        "<div style='background: linear-gradient(90deg, #667eea 0%, #764ba2 100%); padding: 15px; border-radius: 6px; color: white; text-align: center;'>",
        "<div style='font-size: 12px; opacity: 0.9; margin-bottom: 5px;'>FINAL RESULT FOR VISUALIZATION</div>",
        "<div style='font-size: 28px; font-weight: bold;'>",
        format(n_stage_3, big.mark = ","),
        " <span style='font-size: 16px; opacity: 0.9;'>/ ", format(n_stage_0, big.mark = ","), "</span>",
        "</div>",
        "<div style='font-size: 12px; opacity: 0.9; margin-top: 5px;'>",
        round(n_stage_3/n_stage_0*100, 1), "% of dataset",
        "</div>",
        "</div>",
        
        "</div>"
      )
    )
  })

  observe({
    subset_indices_reactive()  # Trigger the reactive
    cat("🔄 Visible categories updated - recalculating subset\n")
  })

  # When lasso selection changes, plots update
  observe({
    subset_indices_reactive()  # Trigger the reactive
    cat("🔄 Lasso selection updated - recalculating subset\n")
  })

  # When numeric filters change, plots update
  observe({
    filtered_indices()  # Dependency
    subset_indices_reactive()  # Trigger the reactive
    cat("🔄 Numeric filters updated - recalculating subset\n")
  })


  # Selected points output
  output$selected_points <- renderPrint({
    if (is.null(input$selectedPoints) || length(input$selectedPoints) == 0) {
      "No cells selected. Use lasso tool to select cells."
    } else {
      paste0(
        "Selected cells: ", length(input$selectedPoints$selectedIndices),
        "\nIndices (first 15): ", paste(head(input$selectedPoints$selectedIndices, 15), collapse = ", "),
        if (length(input$selectedPoints$selectedIndices) > 15) " ..."
      )
    }
  })

  observeEvent(filtered_indices(), {
    filtered <- filtered_indices()
    
    if (is.null(filtered)) {
      cat("📤 NO NUMERIC FILTERS - showing all cells\n")
      session$sendCustomMessage("updateFilteredIndices", list(
        filteredIndices = NULL,
        count = nrow(values$lazy_data$df),
        active = FALSE
      ))
    } else {
      cat("📤 NUMERIC FILTERS ACTIVE -", length(filtered), "cells pass\n")
      session$sendCustomMessage("updateFilteredIndices", list(
        filteredIndices = as.integer(filtered - 1),  # 0-based for JS
        count = length(filtered),
        active = TRUE
      ))
    }
  }, priority = 100)  # High priority to update plots immediately

  observeEvent(input$clearAllFilters, {
    cat("🔄 Clearing all filters\n")
    filtered_indices(NULL)
    
    session$sendCustomMessage("updateFilteredIndices", list(
      filteredIndices = NULL,
      count = nrow(values$lazy_data$df),
      active = FALSE
    ))
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
      
      # Skip columns with only one unique value
      if (nuniq <= 1) {
        return(list(
          categorical = FALSE,
          numeric = FALSE
        ))
      }
      
      # Detect different data types
      is_factor_or_char <- is.factor(col) || is.character(col)
      is_numeric <- is.numeric(col)
      
      # Check if numeric values are integer-like (natural numbers)
      is_integer_like <- FALSE
      if (is_numeric && all(!is.na(col))) {
        is_integer_like <- all(abs(col - round(col)) < .Machine$double.eps^0.5)
      }
      
      # Categorical: factors, characters, OR integer-like numbers with ≤20 categories
      is_categorical <- (is_factor_or_char || is_integer_like) && nuniq <= 20
      
      # Numeric: only continuous (non-integer-like) numeric data
      is_truly_numeric <- is_numeric && !is_integer_like && nuniq > 1
      
      list(
        categorical = is_categorical,
        numeric = is_truly_numeric
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

  # --- Choose Matrix ---

  output$chooseMatrixUI <- renderUI({
    req(numeric_obs_keys())
    selectInput("chooseMatrix", "Choose Matrix:",
      choices = c('X', values$lazy_data$available_layers),
      selected = 'X'
    )
  })

  matrix_for_viz <- reactive({
    input$chooseMatrix
  })

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
      
      # Use the refined classification from obs_stats
      stats <- obs_stats()
      
      # Determine var_type using your refined logic
      if (input$colorBy %in% names(stats)) {
        var_type <- ifelse(stats[[input$colorBy]]$numeric, "gene", "categorical")
      } else {
        # Fallback for columns not in obs_stats (shouldn't happen normally)
        var_type <- "categorical"
      }
      
      ann_colors <- generate_colors(length(unique(annotation_data)))
      
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


  # Initialize reactive value for manual gene sets
  manual_genesets <- reactiveVal(list())

  # Add manual gene set
  observeEvent(input$add_manual_geneset, {
    req(input$manual_geneset_name, input$manual_geneset_genes)
    
    name <- trimws(input$manual_geneset_name)
    genes <- toupper(trimws(unlist(strsplit(input$manual_geneset_genes, ",|\\s+"))))
    genes <- genes[genes != ""]
    
    if (name != "" && length(genes) > 0) {
      current_sets <- manual_genesets()
      current_sets[[name]] <- genes
      manual_genesets(current_sets)
      
      # Clear inputs
      updateTextInput(session, "manual_geneset_name", value = "")
      updateTextInput(session, "manual_geneset_genes", value = "")
      
      # Update selectize choices and select all gene sets
      updateSelectizeInput(session, "selected_manual_genesets", 
                          choices = names(current_sets),
                          selected = names(current_sets)) # Select all gene sets
    }
  })

  # Display manual gene sets
  output$manual_genesets_display <- renderUI({
    sets <- manual_genesets()
    if (length(sets) == 0) {
      return(p("No custom gene sets added yet.", style = "color: #666;"))
    }
    
    set_displays <- lapply(names(sets), function(name) {
      genes_text <- paste(sets[[name]], collapse = ", ")
      if (nchar(genes_text) > 50) {
        genes_text <- paste0(substr(genes_text, 1, 47), "...")
      }
      
      div(
        style = "margin-bottom: 5px; padding: 5px; background-color: #f8f9fa; border-left: 3px solid #007bff;",
        strong(name), ": ", genes_text,
        actionButton(paste0("remove_", name), "×", 
                    class = "btn btn-danger btn-xs pull-right",
                    style = "padding: 1px 6px; font-size: 10px;",
                    onclick = paste0("Shiny.setInputValue('remove_manual_geneset', '", name, "');"))
      )
    })
    
    do.call(tagList, set_displays)
  })

  # Remove manual gene set
  observeEvent(input$remove_manual_geneset, {
    current_sets <- manual_genesets()
    current_sets[[input$remove_manual_geneset]] <- NULL
    manual_genesets(current_sets)
    
    # Update selectize choices and select all remaining gene sets
    updateSelectizeInput(session, "selected_manual_genesets", 
                        choices = names(current_sets),
                        selected = names(current_sets)) # Select all remaining gene sets
  })

  # Check if manual gene sets exist
  output$has_manual_genesets <- reactive({
    return(length(manual_genesets()) > 0)
  })
  outputOptions(output, "has_manual_genesets", suspendWhenHidden = FALSE)

  # Ensure all gene sets are selected whenever manual_genesets changes
  observe({
    current_sets <- manual_genesets()
    updateSelectizeInput(session, "selected_manual_genesets",
                        choices = names(current_sets),
                        selected = names(current_sets)) # Select all gene sets
  })
  # UI select for clustering (with alphabetical ordering)
  output$heatmap_cluster_by_ui_panel1 <- renderUI({
    req(values$lazy_data$obs_keys)
    choices <- sort(values$lazy_data$obs_keys)
    selectizeInput(
      "heatmap_cluster_by_panel1",
      "Cluster or color by:",
      choices = choices,
      options = list(
        placeholder = "Select a grouping variable..."
        # onInitialize = I('function() { this.setValue(""); }')
      )
    )
  })

  output$heatmap_cluster_by_ui_panel2 <- renderUI({
    req(values$lazy_data$obs_keys)
    choices <- sort(values$lazy_data$obs_keys)
    selectizeInput(
      "heatmap_cluster_by_panel2",
      "Cluster or color by:",
      choices = choices,
      options = list(
        placeholder = "Select a grouping variable..."
        # onInitialize = I('function() { this.setValue(""); }')
      )
    )
  })

  # GMT file loading indicator (Panel 2)
  output$gmt_loaded_panel2 <- reactive({
    return(!is.null(input$gmt_file_panel2))
  })
  outputOptions(output, "gmt_loaded_panel2", suspendWhenHidden = FALSE)

  # Parse .gmt file (Panel 2)
  gmt_data_panel2 <- reactive({
    req(input$gmt_file_panel2)
    gmt_file <- input$gmt_file_panel2$datapath
    gmt <- readLines(gmt_file)
    gene_sets <- lapply(gmt, function(line) {
      parts <- unlist(strsplit(line, "\t"))
      list(name = parts[1], genes = toupper(parts[3:length(parts)]))
    })
    names(gene_sets) <- sapply(gene_sets, function(x) x$name)
    gene_sets
  })

  output$umap_gmt_select_panel2 <- renderUI({
    req(gmt_data_panel2())
    gmt_sets <- gmt_data_panel2()
    set_names <- names(gmt_sets)
    
    display_names <- sapply(set_names, function(x) {
      if (nchar(x) > 60) {
        paste0(substr(x, 1, 57), "...")
      } else {
        x
      }
    })
    names(display_names) <- set_names
    
    selectizeInput(
      "umap_gmt_sets_panel2",
      "Select gene sets:",
      choices = display_names,
      multiple = TRUE,
      options = list(
        placeholder = "Select one or more gene sets...",
        maxOptions = 500,
        onInitialize = I('function() { this.setValue(""); }')
      )
    )
  })
  
  # --- Button-controlled pathway scoring and UMAP rendering ---

  # Disable button while computing (optional but recommended)
  observeEvent(input$update_umap_panel2, {
    shinyjs::disable("update_umap_panel2")
    on.exit(shinyjs::enable("update_umap_panel2"))
  })

  # 1️⃣ Pathway scores: only compute when button pressed
  pathway_scores <- eventReactive(input$update_umap_panel2, {
    req(input$umap_gmt_sets_panel2, gmt_data_panel2())
    
    selected_sets <- input$umap_gmt_sets_panel2
    all_gmt <- gmt_data_panel2()
    gene_groups <- all_gmt[selected_sets]
    req(length(gene_groups) > 0)

    # --- Safe cell count ---
    n_cells_safe <- function() {
      tryCatch({
        if (!is.null(values$lazy_data$zarr_obj)) {
          obs_keys <- py_to_r(reticulate::import_builtins()$list(values$lazy_data$zarr_obj[["obs"]]$keys()))
          if (length(obs_keys) > 0) {
            first_obs <- values$lazy_data$zarr_obj[["obs"]][[obs_keys[1]]]
            return(as.integer(py_to_r(first_obs$shape[1])))
          }
        }
        if (is.numeric(values$lazy_data$n_cells)) {
          return(as.integer(values$lazy_data$n_cells))
        }
        return(1000L)
      }, error = function(e) {
        cat("⚠️ Error getting n_cells, using fallback: 1000\n")
        return(1000L)
      })
    }
    
    n_cells_default <- n_cells_safe()

    subset_indices <- isolate({
      if (!is.null(subset_indices_reactive())) subset_indices_reactive()$combined else NULL
    })
    
    # Collect all unique genes from selected gene sets
    all_geneset_genes <- unique(unlist(gene_groups))
    valid_genes <- intersect(toupper(all_geneset_genes), toupper(values$lazy_data$genes))
    
    scores_matrix <- NULL
    if (length(valid_genes) > 0) {
      # Get expression for all genes at once (all cells for UMAP)
      all_expr_mat <- isolate(values$lazy_data$get_genes_expression_batch(
        valid_genes, 
        layer = matrix_for_viz(), 
        cells = subset_indices
      ))
      
      print(paste("Expression matrix loaded with", nrow(all_expr_mat), "cells and", ncol(all_expr_mat), "genes"))
      
      if (!is.null(all_expr_mat) && nrow(all_expr_mat) > 0) {
        gene_lookup <- setNames(1:ncol(all_expr_mat), toupper(colnames(all_expr_mat)))
        n_cells_actual <- as.integer(nrow(all_expr_mat))
        
        print(paste("Calculating scores for", length(gene_groups), "gene sets across", n_cells_actual, "cells"))
        
        mat <- sapply(seq_along(gene_groups), function(i) {
          if (is.list(gene_groups[[i]]) && "genes" %in% names(gene_groups[[i]])) {
            geneset_genes <- intersect(toupper(gene_groups[[i]]$genes), names(gene_lookup))
          } else {
            geneset_genes <- intersect(toupper(gene_groups[[i]]), names(gene_lookup))
          }
          
          if (length(geneset_genes) > 0) {
            gene_indices <- gene_lookup[geneset_genes]
            gene_indices <- gene_indices[!is.na(gene_indices)]
            
            if (length(gene_indices) > 0) {
              expr_subset <- all_expr_mat[, gene_indices, drop = FALSE]
              if (input$umap_score_method_panel2 == "mean") {
                rowMeans(expr_subset, na.rm = TRUE)
              } else {
                rowSums(expr_subset, na.rm = TRUE)
              }
            } else {
              rep(0, n_cells_actual)
            }
          } else {
            rep(0, n_cells_actual)
          }
        })
        scores_matrix <- round(mat, 2)
        colnames(scores_matrix) <- names(gene_groups)
      } else {
        print("⚠️ Fallback: Using zero matrix")
        scores_matrix <- matrix(0, nrow = n_cells_default, ncol = length(gene_groups))
        colnames(scores_matrix) <- names(gene_groups)
      }
    } else {
      print("⚠️ No valid genes found in selected gene sets, using zero matrix")
      scores_matrix <- matrix(0, nrow = n_cells_default, ncol = length(gene_groups))
      colnames(scores_matrix) <- names(gene_groups)
    }

    cat("Pathway scores updated for", length(selected_sets), "sets\n")
    scores_matrix
  })


  # 2️⃣ UMAP data: also only recompute when button is pressed
  umap_data <- eventReactive(input$update_umap_panel2, {
    req(pathway_scores())

    subset_indices <- isolate({
      if (!is.null(subset_indices_reactive())) {
        subset_indices_reactive()$combined
      } else {
        seq_len(nrow(values$lazy_data$df))
      }
    })

    umap_df <- isolate(values$lazy_data$df[subset_indices, , drop = FALSE])
    req(!is.null(umap_df) && nrow(umap_df) >= 2)

    umap_coords <- umap_df[, 1:2]
    colnames(umap_coords) <- c("x", "y")

    n_cells_umap <- nrow(umap_coords)
    n_cells_scores <- nrow(pathway_scores())

    print(paste("UMAP cells:", n_cells_umap, "Scores cells:", n_cells_scores))

    if (n_cells_umap != n_cells_scores) {
      if (n_cells_umap > n_cells_scores) {
        umap_coords <- umap_coords[1:n_cells_scores, ]
      } else {
        extra_scores <- matrix(0, nrow = n_cells_umap - n_cells_scores, ncol = ncol(pathway_scores()))
        padded_scores <- rbind(pathway_scores(), extra_scores)
      }
    }

    # Combine UMAP coords with all pathway scores
    plot_data <- cbind(umap_coords, as.data.frame(pathway_scores()))
    
    # Get list of pathways (up to 4 max)
    pathways <- colnames(pathway_scores())
    pathways_to_plot <- head(pathways, 4)
    
    list(data = plot_data, pathways = pathways_to_plot)
  })


  # 3️⃣ Track sync state for multi-plot
  multi_sync_enabled <- reactiveVal(FALSE)
  
  observeEvent(input$toggle_multi_sync, {
    req(umap_data())
    ud <- umap_data()
    
    multi_sync_enabled(!multi_sync_enabled())
    
    # Generate plot IDs based on current pathways
    plot_ids <- paste0("umapPlot_", seq_along(ud$pathways))
    
    cat("🔗 Sync toggled:", multi_sync_enabled(), "| Plot IDs:", paste(plot_ids, collapse=", "), "\n")
    
    session$sendCustomMessage("my_scatterplot_sync", list(
      enabled = multi_sync_enabled(),
      plotIds = plot_ids
    ))
  })


  # 4️⃣ Render UI container for multi-plots
  output$umapPlotContainer <- renderUI({
    req(umap_data())
    ud <- umap_data()
    n_pathways <- length(ud$pathways)
    
    if (n_pathways == 0) {
      return(div("No pathways to display", style = "color: #999; padding: 20px;"))
    }
    
    # Create grid layout
    grid_style <- if (n_pathways == 1) "grid-template-columns: 1fr;" else
                  if (n_pathways == 2) "grid-template-columns: 1fr 1fr;" else
                  "grid-template-columns: 1fr 1fr;"
    
    plot_outputs <- lapply(seq_along(ud$pathways), function(i) {
      pathway <- ud$pathways[i]
      plot_id <- paste0("umapPlot_", i)
      
      div(
        class = "plot-container",
        style = "border: 1px solid #ddd; border-radius: 4px; overflow: hidden; background: white;",
        my_scatterplotOutput(
          outputId = plot_id,
          width = "100%",
          height = "400px"
        )
      )
    })
    
    # Add sync button if multiple plots
    sync_button_ui <- NULL
    if (n_pathways > 1) {
      sync_button_ui <- actionButton(
        "toggle_multi_sync",
        if (multi_sync_enabled()) "🔗 Disable Sync" else "Enable Sync",
        style = paste(
          "margin-bottom: 10px;",
          "padding: 10px 20px;",
          "background:", if (multi_sync_enabled()) "#28a745" else "#6c757d", ";",
          "color: white;",
          "border: none;",
          "border-radius: 4px;",
          "cursor: pointer;",
          "font-weight: bold;"
        )
      )
    }
    
    div(
      sync_button_ui,
      div(
        style = paste("display: grid; gap: 15px;", grid_style),
        plot_outputs
      )
    )
  })


  # 5️⃣ Render individual pathway plots (up to 4)
  observe({
    req(umap_data())
    ud <- umap_data()
    
    cat("🔄 Rendering", length(ud$pathways), "plots\n")
    
    # Render up to 4 individual plots
    for (i in 1:min(4, length(ud$pathways))) {
      local({
        idx <- i
        pathway <- ud$pathways[idx]
        plot_id <- paste0("umapPlot_", idx)
        
        cat("Rendering plot", idx, "with pathway:", pathway, "\n")
        
        output[[plot_id]] <- renderMy_scatterplot({
          req(umap_data())
          ud_inner <- umap_data()
          
          plot_data <- data.frame(
            x = ud_inner$data$x,
            y = ud_inner$data$y,
            expression = ud_inner$data[[pathway]]
          )
          
          size_calc <- log(nrow(plot_data)) + 3
          
          my_scatterplot(
            data = plot_data,
            x = "x",
            y = "y",
            colorBy = "expression",
            size = size_calc,
            continuous_palette = "inferno",
            xlab = "UMAP 1",
            ylab = "UMAP 2",
            showAxes = FALSE,
            showTooltip = TRUE,
            opacity = 0.8,
            backgroundColor = "white",
            legend_title = pathway,
            enableDownload = TRUE,
            plotId = plot_id,
            syncPlots = paste0("umapPlot_", seq_along(ud_inner$pathways))
          )
        })
      })
    }
  })

  # UI select for GMT gene sets (Panel 2) - with "All" option
  output$heatmap_gmt_select_panel2 <- renderUI({
    req(gmt_data_panel2())
    gmt_sets <- gmt_data_panel2()
    set_names <- names(gmt_sets)
    
    display_names <- sapply(set_names, function(x) {
      if (nchar(x) > 60) {
        paste0(substr(x, 1, 57), "...")
      } else {
        x
      }
    })
    names(display_names) <- set_names
    
    # Add "All" option at the beginning
    all_choices <- c("All gene sets" = "ALL", display_names)
    
    selectizeInput(
      "heatmap_gmt_sets_panel2",
      "Select gene sets:",
      choices = all_choices,
      multiple = TRUE,
      options = list(
        placeholder = "Select one or more gene sets (or 'All gene sets')...",
        maxOptions = 500,
        onInitialize = I('function() { this.setValue(""); }'),
        onItemAdd = I('function(value, item) {
          if (value === "ALL") {
            // If "All" is selected, remove all other selections
            this.setValue(["ALL"]);
          } else if (this.getValue().includes("ALL")) {
            // If other item is selected and "All" is already selected, remove "All"
            var current = this.getValue();
            var filtered = current.filter(function(v) { return v !== "ALL"; });
            filtered.push(value);
            this.setValue(filtered);
          }
        }')
      )
    )
  })

  # Heatmap data for Panel 1 (Individual Genes) - with separate update triggers and ignoreNULL
  heatmap_data_panel1 <- eventReactive({
    input$update_heatmap_panel1
    input$update_both_heatmaps
  }, {
    withProgress(message = "Updating Panel 1 heatmap...", value = 0, {
      req(input$heatmap_genes_panel1)
      incProgress(0.1)

      genes <- toupper(trimws(unlist(strsplit(input$heatmap_genes_panel1, ",|\\s+"))))
      genes <- genes[genes %in% toupper(values$lazy_data$genes)]
      req(length(genes) > 0, "No valid genes found in dataset for Panel 1")

      if (!is.null(subset_indices_reactive())) {
        subset_indices <- subset_indices_reactive()$combined
      } else {
        subset_indices <- NULL
      }

      expr_mat <- values$lazy_data$get_genes_expression_batch(genes, layer = matrix_for_viz(), cells = subset_indices)

      if (is.null(expr_mat)) {
        mat <- matrix(NA, nrow = values$lazy_data$n_cells, ncol = length(genes))
        colnames(mat) <- genes
      } else {
        mat <- expr_mat
      }
      
      incProgress(0.5)
      process_heatmap_data(mat, input$heatmap_score_panel1, input$heatmap_cell_group_mode_panel1, input$heatmap_cluster_by_panel1, subset_indices)
    })
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  # Heatmap data for Panel 2 (Gene Sets/Pathways) - with separate update triggers and ignoreNULL
  heatmap_data_panel2 <- eventReactive({
    input$update_heatmap_panel2
    input$update_both_heatmaps
  }, {
    withProgress(message = "Updating Panel 2 heatmap...", value = 0, {
      req(input$heatmap_gene_group_mode_panel2)
      incProgress(0.1)

      if (!is.null(subset_indices_reactive())) {
        subset_indices <- subset_indices_reactive()$combined
      } else {
        subset_indices <- NULL
      }

      gene_groups <- NULL
      gene_group_names <- NULL
      
      # Get actual number of cells from the data structure
      n_cells_safe <- function() {
        tryCatch({
          # Try to get from zarr object first
          if (!is.null(values$lazy_data$zarr_obj)) {
            obs_keys <- py_to_r(reticulate::import_builtins()$list(values$lazy_data$zarr_obj[["obs"]]$keys()))
            if (length(obs_keys) > 0) {
              first_obs <- values$lazy_data$zarr_obj[["obs"]][[obs_keys[1]]]
              return(as.integer(py_to_r(first_obs$shape[1])))
            }
          }
          # Fallback to stored value
          if (is.numeric(values$lazy_data$n_cells)) {
            return(as.integer(values$lazy_data$n_cells))
          }
          # Last resort fallback
          return(1000L)
        }, error = function(e) {
          cat("⚠️ Error getting n_cells, using fallback: 1000\n")
          return(1000L)
        })
      }
      
      n_cells_default <- n_cells_safe()
      
      if (input$heatmap_gene_group_mode_panel2 == "pathways") {
        req(input$heatmap_pathway_select_panel2)
        pathways <- values$lazy_data$pathways
        
        # Handle "All" option
        if ("ALL" %in% input$heatmap_pathway_select_panel2) {
          gene_groups <- pathways
        } else {
          gene_groups <- pathways[input$heatmap_pathway_select_panel2]
        }
        
        # Safe name processing with null checks
        gene_group_names <- sapply(names(gene_groups), function(x) {
          if (is.null(x) || is.na(x)) {
            "Unknown"
          } else if (nchar(x) > 40) {
            paste0(substr(x, 1, 37), "...")
          } else {
            x
          }
        })
        
      } else if (input$heatmap_gene_group_mode_panel2 == "custom_gmt") {
        req(input$heatmap_gmt_sets_panel2)
        all_gmt <- gmt_data_panel2()
        
        # Handle "All" option for GMT
        if ("ALL" %in% input$heatmap_gmt_sets_panel2) {
          gene_groups <- all_gmt
        } else {
          gene_groups <- all_gmt[input$heatmap_gmt_sets_panel2]
        }
        
        # Safe name processing with null checks
        gene_group_names <- sapply(names(gene_groups), function(x) {
          if (is.null(x) || is.na(x)) {
            "Unknown"
          } else if (nchar(x) > 40) {
            paste0(substr(x, 1, 37), "...")
          } else {
            x
          }
        })
        
      } else if (input$heatmap_gene_group_mode_panel2 == "manual") {
        req(input$selected_manual_genesets)
        all_manual <- manual_genesets()
        gene_groups <- all_manual[input$selected_manual_genesets]
        gene_group_names <- names(gene_groups)
      }
      
      req(length(gene_groups) > 0, "No gene sets selected for Panel 2")
      incProgress(0.3)

      # Calculate gene set scores - OPTIMIZED for "All" option
      if (input$heatmap_gene_group_mode_panel2 == "pathways" && "ALL" %in% input$heatmap_pathway_select_panel2) {
        # Ultra-fast batch processing for all pathways
        cat("🚀 Processing ALL pathways with batch method...\n")
        
        # Get all unique genes from all pathways
        all_pathway_genes <- unique(unlist(gene_groups))
        valid_genes <- intersect(toupper(all_pathway_genes), toupper(values$lazy_data$genes))
        
        if (length(valid_genes) > 0) {
          # Get expression for all genes at once
          all_expr_mat <- values$lazy_data$get_genes_expression_batch(valid_genes, layer = matrix_for_viz(), cells = subset_indices)
          
          if (!is.null(all_expr_mat) && nrow(all_expr_mat) > 0) {
            # Create gene lookup for fast indexing
            gene_lookup <- setNames(1:ncol(all_expr_mat), toupper(colnames(all_expr_mat)))
            n_cells_actual <- as.integer(nrow(all_expr_mat))
            
            # Calculate pathway scores using pre-loaded expression matrix
            mat <- sapply(seq_along(gene_groups), function(i) {
              pathway_genes <- intersect(toupper(gene_groups[[i]]), names(gene_lookup))
              if (length(pathway_genes) > 0) {
                # Get indices of pathway genes in the expression matrix
                gene_indices <- gene_lookup[pathway_genes]
                gene_indices <- gene_indices[!is.na(gene_indices)]
                
                if (length(gene_indices) > 0) {
                  rowMeans(all_expr_mat[, gene_indices, drop = FALSE], na.rm = TRUE)
                } else {
                  rep(0, n_cells_actual)
                }
              } else {
                rep(0, n_cells_actual)
              }
            })
          } else {
            # Fallback if matrix loading failed
            mat <- matrix(0, nrow = n_cells_default, ncol = length(gene_groups))
          }
        } else {
          # No valid genes found
          mat <- matrix(0, nrow = n_cells_default, ncol = length(gene_groups))
        }
      } else if (input$heatmap_gene_group_mode_panel2 == "custom_gmt" && "ALL" %in% input$heatmap_gmt_sets_panel2) {
        # Ultra-fast batch processing for all GMT gene sets
        cat("🚀 Processing ALL GMT gene sets with batch method...\n")
        
        # Get all unique genes from all gene sets
        all_geneset_genes <- unique(unlist(lapply(gene_groups, function(x) {
          if (is.list(x) && "genes" %in% names(x)) x$genes else x
        })))
        valid_genes <- intersect(toupper(all_geneset_genes), toupper(values$lazy_data$genes))
        
        if (length(valid_genes) > 0) {
          # Get expression for all genes at once
          all_expr_mat <- values$lazy_data$get_genes_expression_batch(valid_genes, layer = matrix_for_viz(), cells = subset_indices)
          
          if (!is.null(all_expr_mat) && nrow(all_expr_mat) > 0) {
            # Create gene lookup for fast indexing
            gene_lookup <- setNames(1:ncol(all_expr_mat), toupper(colnames(all_expr_mat)))
            n_cells_actual <- as.integer(nrow(all_expr_mat))
            
            # Calculate gene set scores using pre-loaded expression matrix
            mat <- sapply(seq_along(gene_groups), function(i) {
              if (is.list(gene_groups[[i]]) && "genes" %in% names(gene_groups[[i]])) {
                geneset_genes <- intersect(toupper(gene_groups[[i]]$genes), names(gene_lookup))
              } else {
                geneset_genes <- intersect(toupper(gene_groups[[i]]), names(gene_lookup))
              }
              
              if (length(geneset_genes) > 0) {
                # Get indices of gene set genes in the expression matrix
                gene_indices <- gene_lookup[geneset_genes]
                gene_indices <- gene_indices[!is.na(gene_indices)]
                
                if (length(gene_indices) > 0) {
                  rowMeans(all_expr_mat[, gene_indices, drop = FALSE], na.rm = TRUE)
                } else {
                  rep(0, n_cells_actual)
                }
              } else {
                rep(0, n_cells_actual)
              }
            })
          } else {
            # Fallback if matrix loading failed
            mat <- matrix(0, nrow = n_cells_default, ncol = length(gene_groups))
          }
        } else {
          # No valid genes found
          mat <- matrix(0, nrow = n_cells_default, ncol = length(gene_groups))
        }
      } else {
        # Standard processing for individual selections
        mat <- sapply(seq_along(gene_groups), function(i) {
          if (is.list(gene_groups[[i]]) && "genes" %in% names(gene_groups[[i]])) {
            valid_genes <- intersect(gene_groups[[i]]$genes, toupper(values$lazy_data$genes))
          } else {
            valid_genes <- intersect(toupper(gene_groups[[i]]), toupper(values$lazy_data$genes))
          }
          
          if (length(valid_genes) > 0) {
            expr_mat <- values$lazy_data$get_genes_expression_batch(valid_genes, layer = matrix_for_viz(), cells = subset_indices)
            if (!is.null(expr_mat)) {
              rowMeans(expr_mat, na.rm = TRUE)
            } else {
              rep(0, n_cells_default)
            }
          } else {
            rep(0, n_cells_default)
          }
        })
      }
      colnames(mat) <- gene_group_names
      
      incProgress(0.5)
      process_heatmap_data(mat, input$heatmap_score_panel2, input$heatmap_cell_group_mode_panel2, input$heatmap_cluster_by_panel2, subset_indices)
    })
  }, ignoreNULL = FALSE, ignoreInit = FALSE)

  # Helper function to process heatmap data
  process_heatmap_data <- function(mat, score_type, cell_group_mode, cluster_by, cells) {
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
      group_var <- values$lazy_data$get_obs_column(values$lazy_data$zarr_obj, cluster_by, cells)
      req(length(group_var) == nrow(mat))
      
      group_levels <- sort(unique(group_var))
      group_var <- factor(group_var, levels = group_levels)
      
      mat <- sapply(colnames(mat), function(g) tapply(mat[, g], group_var, mean, na.rm = TRUE))
      mat <- t(mat)
    } else if (cell_group_mode == "none") {
      req(cluster_by)
      group_var <- values$lazy_data$get_obs_column(values$lazy_data$zarr_obj, cluster_by, cells)
      req(length(group_var) == nrow(mat))
      
      group_levels <- sort(unique(group_var))
      group_var <- factor(group_var, levels = group_levels)
      o <- order(group_var)
      
      mat <- t(mat[o, , drop = FALSE])
      group_info <- group_var[o]
    }

    list(mat = mat, group_info = group_info)
  }

  # Panel 1
  heatmapPlot1_obj <- reactive({
    req(heatmap_data_panel1())
    data <- heatmap_data_panel1()
    p <- render_heatmap(data, NULL)
    plot_storage$heatmap1 <- p
    p
  })

  output$heatmapPlot1 <- renderPlot(
    {
      req(heatmapPlot1_obj())
      heatmapPlot1_obj()
    },
    height = function() {
      n_rows <- nrow(heatmap_data_panel1()$mat)
      calc_heatmap_height(n_rows)$px  # <-- use px
    }
  )

  # Panel 2
  output$heatmapPlot2 <- renderPlot(
    {
      req(heatmap_data_panel2())
      data <- heatmap_data_panel2()
      p <- render_heatmap(data, NULL)
      plot_storage$heatmap2 <- p
      p
    },
    height = function() {
      n_rows <- nrow(heatmap_data_panel2()$mat)
      calc_heatmap_height(n_rows)$px
    }
  )

  calc_heatmap_height <- function(n_rows, row_cm = 0.5, row_px = 20,
                                  min_cm = 6, max_cm = 40,
                                  min_px = 300, max_px = 2000) {
    list(
      cm = max(min(n_rows * row_cm, max_cm), min_cm),
      px = max(min(n_rows * row_px, max_px), min_px)
    )
  }

  calc_heatmap_size <- function(n_rows, n_cols,
                                row_cm = 0.5, col_cm = 0.3,
                                min_cm = 6, max_cm = 40,
                                min_w = 6, max_w = 50) {
    list(
      height_cm = max(min(n_rows * row_cm, max_cm), min_cm),
      width_cm  = max(min(n_cols * col_cm, max_w), min_w)
    )
  }

  render_heatmap <- function(data, title) {
    mat <- data$mat
    n_rows <- nrow(mat)
    n_cols <- ncol(mat)
    sz <- calc_heatmap_size(n_rows, n_cols)

    library(ComplexHeatmap)
    library(RColorBrewer)

    if (!is.null(data$group_info)) {
      groups <- levels(data$group_info)
      max_colors <- 9
      base_pal <- brewer.pal(max_colors, "Set1")
      n <- length(groups)

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

      ht <- Heatmap(
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
      ht <- Heatmap(
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

    # Only fix width, let Shiny control height
    draw(ht, heatmap_width = unit(sz$width_cm, "cm"))
  }


  # Gene Tab

  # UI: Gene search for violin_plot_genes tab
  output$gene_selected_ui <- renderUI({
    message("Rendering gene_selected_ui") # Debug
    if (is.null(values$lazy_data)) {
      selectizeInput("geneSearch_tab", "🔍 Search Gene:",
                    choices = list("Initialize dataset first..." = "loading"),
                    multiple = TRUE
      )
    } else {
      selectizeInput("geneSearch_tab", "🔍 Search Gene:",
                    choices = NULL, # Rely on server-side updates
                    multiple = TRUE,
                    options = list(
                      placeholder = "Type gene names ...",
                      openOnFocus = FALSE,
                      closeAfterSelect = TRUE,
                      plugins = list("remove_button"),
                      maxItems = 15,
                      onInitialize = I("function() { console.log('geneSearch_tab initialized'); Shiny.setInputValue('geneSearch_tab_initialized', true, {priority: 'event'}); }") # Trigger update
                    )
      )
    }
  })

  # Store selections in reactive values
  values$geneSearch_selection <- NULL
  values$geneSearch_tab_selection <- NULL

  # Update stored selections when user changes them
  observe({
    values$geneSearch_selection <- input$geneSearch
  }) %>% bindEvent(input$geneSearch, ignoreNULL = TRUE, ignoreInit = TRUE)

  observe({
    values$geneSearch_tab_selection <- input$geneSearch_tab
  }) %>% bindEvent(input$geneSearch_tab, ignoreNULL = TRUE, ignoreInit = TRUE)

  # Reset gene selections and cache when dataset changes
  observe({
    message("Dataset changed, resetting gene selections and cache") # Debug
    values$geneSearch_selection <- NULL
    values$geneSearch_tab_selection <- NULL
    gene_data_cache$data <- list()  # Clear entire cache
    displayed_genes$genes <- character(0)  # Reset displayed genes
  }) %>% bindEvent(input$dataset, ignoreNULL = TRUE, ignoreInit = TRUE)

  # Update gene choices when data is loaded, resetting selections
  observe({
    req(values$lazy_data) # Ensure dataset is loaded
    message("Updating selectizeInputs for geneSearch. Gene list length: ", length(values$lazy_data$genes)) # Debug

    updateSelectizeInput(
      session,
      "geneSearch",
      choices = values$lazy_data$genes,
      selected = NULL, # Reset selections
      server = TRUE
    )
  }) %>% bindEvent(values$lazy_data, ignoreNULL = TRUE, ignoreInit = FALSE)

  # Update gene choices for geneSearch_tab when data is loaded or tab is initialized
  observe({
    req(values$lazy_data, input$geneSearch_tab_initialized) # Trigger when geneSearch_tab is initialized or data changes
    message("Updating selectizeInputs for geneSearch_tab. Gene list length: ", length(values$lazy_data$genes)) # Debug

    updateSelectizeInput(
      session,
      "geneSearch_tab",
      choices = values$lazy_data$genes,
      selected = NULL, # Reset selections
      server = TRUE
    )
  }) %>% bindEvent(list(values$lazy_data, input$geneSearch_tab_initialized), ignoreNULL = TRUE, ignoreInit = FALSE)

  # Initialize reactiveValues for caching gene-specific data and displayed genes
  gene_data_cache <- reactiveValues(data = list())
  displayed_genes <- reactiveValues(genes = character(0))  # Track genes to display in UI

  # Group by UI (unchanged)
  output$gene_group_by_ui <- renderUI({
    choices <- categorical_vars()
    if (length(choices) == 0) {
      return(p("No categorical variables available for grouping."))
    }
    selectInput("gene_group_by", "Group by:",
                choices = choices,
                selected = choices[1])
  })

  # Update displayed genes when button is clicked
  observeEvent(input$update_plots, {
    req(input$geneSearch_tab)
    displayed_genes$genes <- input$geneSearch_tab  # Update genes to display
  })

  # Reactive data triggered by button
  gene_expression_data <- eventReactive(input$update_plots, {
    req(input$geneSearch_tab, input$gene_group_by)
    genes <- input$geneSearch_tab

    # Get current subset
    current_subset <- if (!is.null(subset_indices_reactive())) {
      subset_indices_reactive()$combined
    } else {
      NULL
    }

    # Update cache for new or changed genes
    for (gene in genes) {
      cached <- gene_data_cache$data[[gene]]
      invalidate_cache <- is.null(cached) || 
                          !identical(input$gene_group_by, cached$group_by) ||
                          !identical(current_subset, cached$subset_indices)

      if (invalidate_cache) {
        # Get expression for this gene
        
        expr <- values$lazy_data$get_gene_expression(gene_name = gene, layer = matrix_for_viz(), cells = current_subset)
        
        # Get grouping variable
        group_var <- values$lazy_data$get_obs_column(values$lazy_data$z, input$gene_group_by, cells = current_subset)
        req(group_var)
        group_var <- as.factor(group_var)
        
        # Build long-format data frame for this gene
        df <- data.frame(group = group_var, expression = expr, gene = gene)
        
        # Cache the data
        gene_data_cache$data[[gene]] <- list(
          data = df,
          group_by = input$gene_group_by,
          subset_indices = current_subset  # Store subset for invalidation
        )
      }
    }
    
    # Remove cache entries for genes no longer selected
    cached_genes <- names(gene_data_cache$data)
    removed_genes <- setdiff(cached_genes, genes)
    for (gene in removed_genes) {
      gene_data_cache$data[[gene]] <- NULL
      plot_storage[[paste0("gene_", gene)]] <- NULL  # Clear plot storage
    }
    
    # Combine cached data for all selected genes
    long_df <- do.call(rbind, lapply(genes, function(g) {
      gene_data_cache$data[[g]]$data
    }))
    return(long_df)
  })

  # Dynamic UI for plots, triggered by button
  output$gene_plots_ui <- renderUI({
    req(input$update_plots, displayed_genes$genes)  # Depend on button and displayed genes
    genes <- displayed_genes$genes
    
    if (length(genes) == 0) {
      return(p("No genes selected. Click 'Update Plots' to display."))
    }
    
    # Create plotOutput for each gene
    plot_outputs <- lapply(genes, function(gene) {
      plotOutput(paste0("gene_plot_", gene), height = "300px")
    })
    
    do.call(tagList, plot_outputs)
  })

  # Render individual plots when button is clicked
  observeEvent(input$update_plots, {
    req(input$geneSearch_tab, input$gene_group_by, gene_expression_data())
    genes <- input$geneSearch_tab
    
    # Render plots for each gene
    lapply(genes, function(gene) {
      output[[paste0("gene_plot_", gene)]] <- renderPlot({
        # Get cached data for this gene
        gene_df <- gene_data_cache$data[[gene]]$data
        
        # Create plot
        p <- ggplot(gene_df, aes(x = group, y = expression, fill = group)) +
          labs(x = "Group", y = "Expression Level", title = paste(gene)) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1),
                legend.position = "none")
        
        # Add violin or boxplot
        if (input$gene_plot_type == "violin") {
          p <- p + geom_violin(scale = "width", trim = TRUE)
        } else {
          p <- p + geom_boxplot(outlier.size = 0.5)
        }
        
        # Add points if selected
        if (input$gene_show_points) {
          p <- p + geom_jitter(size = input$gene_point_size, width = 0.2, alpha = 0.5)
        }
        
        # Apply log scale if selected
        if (input$gene_log_scale) {
          p <- p + scale_y_continuous(trans = "log1p")
        }
        
        # Store plot
        plot_storage[[paste0("gene_", gene)]] <- p
        p
      })
    })
  })
    
  # NEW REACTIVE EXPRESSION for metadata
  qc_metadata_data <- reactive({
      req(values$lazy_data)
      
      metadata_vals <- NULL

      if (!is.null(subset_indices_reactive())) {
        subset_indices <- subset_indices_reactive()$combined
      } else {
        subset_indices <- NULL
      }
      
      # Robust check for function existence before calling
      if (is.function(values$lazy_data$get_obs_keys)) {
          obs_keys <- values$lazy_data$get_obs_keys(values$lazy_data$z)
          metadata_vars <- c("celltype", "Annotation")
          
          for (var in metadata_vars) {
              if (var %in% obs_keys) {
                  # Return the entire metadata column
                  metadata_vals <- as.character(values$lazy_data$get_obs_column(values$lazy_data$z, var, cells = subset_indices))
                  return(metadata_vals) # Return immediately once found
              }
          }
      }
      
      # If no metadata variable is found, return NULL or a vector of the correct size
      # to avoid row length mismatches. We will handle the length later.
      return(NULL)
  })


  # MODIFIED qc_main_data reactive (with robust length matching)
  qc_main_data <- reactive({
      req(values$lazy_data, input$qc_plot_type)
      

      if (!is.null(subset_indices_reactive())) {
        subset_indices <- subset_indices_reactive()$combined
      } else {
        subset_indices <- NULL
      }

      n_obs <- if (is.null(subset_indices)) nrow(values$lazy_data$df) else length(subset_indices)

      if (input$qc_plot_type == "scatter") {
          req(input$qc_x_metric, input$qc_y_metric, n_obs > 0)  # Guard against empty data

          x <- values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_x_metric, cells = subset_indices)
          y <- values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_y_metric, cells = subset_indices)

          if (is.null(x) || length(x) != n_obs) x <- rep(NA_real_, n_obs)
          if (is.null(y) || length(y) != n_obs) y <- rep(NA_real_, n_obs)

          x <- as.numeric(x)
          y <- as.numeric(y)

          # Group
          group <- NULL
          if (!is.null(input$qc_group_by) && input$qc_group_by != "None") {
            group <- values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_group_by, cells = subset_indices)
          }
          if (is.null(group) || length(group) != n_obs) {
            group <- rep(NA_character_, n_obs)
          }
          group <- as.factor(group)

          # Color
          color_vals <- rep(0L, n_obs)
          if (!is.null(input$qc_color_by) && input$qc_color_by != "None") {
            color_factors <- values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_color_by, cells = subset_indices)
            if (!is.null(color_factors) && length(color_factors) == n_obs) {
              color_vals <- as.integer(as.factor(color_factors)) - 1L
            }
          }

          # Metadata
          metadata_vals <- qc_metadata_data()
          if (is.null(metadata_vals) || length(metadata_vals) != n_obs) {
            metadata_vals <- rep("", n_obs)
          }

          df <- data.frame(
            X = x,
            Y = y,
            Color = color_vals,
            Metadata = metadata_vals,
            Group = group,
            stringsAsFactors = FALSE
          )
          return(df)

      } else if (input$qc_plot_type == "violin") {
          req(input$qc_metric, input$qc_group_by, n_obs > 0)
          metric <- as.numeric(values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_metric, cells = subset_indices))
          group <- as.factor(values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_group_by, cells = subset_indices))

          # check NAs indices in both metric and group and remove them
          na_indices <- is.na(metric) | is.na(group)
          metric <- metric[!na_indices]
          group <- group[!na_indices]

          if (length(metric) != length(group)) {
              stop("Row length mismatch in violin data")
          }
          
          return(data.frame(Metric = metric, Group = group, stringsAsFactors = FALSE))
      }
      
      # Fallback for invalid type
      return(data.frame())
  })


  output$qc_main_plot <- renderUI({
    if (input$qc_plot_type == "violin") {
      tags$div(
        class = "violin_container",
        style = "position: relative; width:800px; height:400px;",
        tags$canvas(id = "violin_plot", width = 800, height = 400, style = "position: absolute; top: 0; left: 0;"),
        tags$svg(id = "violin_overlay", width = 800, height = 400, style = "position: absolute; top: 0; left: 0;")
      )
    } else if (input$qc_plot_type == "scatter") {
      tags$div(id = "scatterplot", style = "width: 100%; height: 500px; position: relative;")
    }
  })

  # Consolidated observer to handle both violin and scatter plots
  observeEvent(list(
      qc_main_data(),
      input$qc_plot_type,
      input$qc_x_metric,
      input$qc_y_metric,
      input$qc_metric,
      input$qc_group_by,
      input$qc_color_by,
      input$qc_show_points,
      input$qc_log_scale,
      input$qc_point_size,
      input$qc_legend_position
  ), {
      # This `req()` is critical. It ensures that all inputs are available
      # and that qc_main_data has successfully computed a value.
      req(input$qc_plot_type, qc_main_data())
      
      dat <- qc_main_data()

      if (!is.null(subset_indices_reactive())) {
        subset_indices <- subset_indices_reactive()$combined
      } else {
        subset_indices <- NULL
      }

      if (input$qc_plot_type == "violin") {
        # Prepare and send data for the violin plot
        datasets <- lapply(unique(dat$Group), function(g) {
            group_data <- dat[dat$Group == g, ]
            list(
                values = group_data$Metric,
                metadata = list(name = g, n = nrow(group_data))
            )
        })
        
        message_to_send <- c(
            list(datasets = datasets),
            list(
                showPoints = input$qc_show_points,
                logScale = input$qc_log_scale,
                pointSize = input$qc_point_size,
                legendPosition = input$qc_legend_position
            )
        )
        
        session$sendCustomMessage("drawViolins", message_to_send)
        plot_storage$qc_main <- TRUE

      } else if (input$qc_plot_type == "scatter") {
        # Prepare and send data for the scatterplot
        scatterplot_matrix <- as.matrix(dat[, c("X", "Y", "Color")])
        
        # New: Get color levels and metadata directly here, after `req`
        color_levels <- NULL
        metadata_vals <- NULL
        
        # Check if the color_by variable exists and is not "None"
        if (!is.null(input$qc_color_by) && input$qc_color_by != "None") {
            # This is now safe because dat (from qc_main_data) is ready
            color_levels <- levels(as.factor(values$lazy_data$get_obs_column(values$lazy_data$z, input$qc_color_by, cells = subset_indices)))
        }

        # Retrieve metadata directly here
        metadata_vars <- c("celltype", "Annotation")
        for (var in metadata_vars) {
            if (var %in% categorical_vars()) { # Your fix is now safe to use here
                metadata_vals <- as.character(values$lazy_data$get_obs_column(values$lazy_data$z, var, cells = subset_indices))
                break
            }
        }
        
        session$sendCustomMessage("scatterplotData", list(
            data = scatterplot_matrix,
            x_label = input$qc_x_metric,
            y_label = input$qc_y_metric,
            color_by = input$qc_color_by,
            color_levels = color_levels,
            point_size = input$qc_point_size,
            metadata = metadata_vals
        ))

        plot_storage$qc_main <- TRUE
      }
  })

  # Keep the separate observer for live violin updates to prevent lag
  observeEvent(list(input$qc_show_points, input$qc_log_scale, input$qc_point_size, input$qc_legend_position), {
    if (input$qc_plot_type == "violin") {
        session$sendCustomMessage("updateViolinInputs", list(
            showPoints = input$qc_show_points,
            logScale = input$qc_log_scale,
            pointSize = input$qc_point_size,
            legendPosition = input$qc_legend_position
        ))
    }
  })


  output$qc_stats_table <- renderTable({
    req(values$lazy_data)

    if (!is.null(subset_indices_reactive())) {
      subset_indices <- subset_indices_reactive()$combined
    } else {
      subset_indices <- NULL
    }

    n_genes_col <- get_first_available_column(values$lazy_data, c("n_genes", "nGenes", "genes_count"))[subset_indices]
    total_counts_col <- get_first_available_column(values$lazy_data, c("total_counts", "totalUMIs"))[subset_indices]
    pct_mito_col <- get_first_available_column(values$lazy_data, c("pct_mito", "percent_mito", "mitochondrial_percent"))[subset_indices]
    pct_ribo_col <- get_first_available_column(values$lazy_data, c("pct_ribo", "percent_ribo", "ribosomal_percent"))[subset_indices]

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

  # Server logic for the enhanced filtering
  observe({
    req(values$lazy_data)

    subset_indices <- isolate({
      if (!is.null(subset_indices_reactive())) {
        subset_indices_reactive()$combined
      } else {
        seq_len(nrow(values$lazy_data$df))
      }
    })
    
    # Get available columns from the dataset
    available_cols <- numeric_obs_keys()
    
    # Update dropdown choices
    updateSelectInput(session, "mito_var_select",
                    choices = c("Choose variable" = "", available_cols),
                    selected = get_first_available_column(values$lazy_data, 
                                                              c("pct_mito", "percent_mito", "mitochondrial_percent", "percent.mt"))[subset_indices])
    
    updateSelectInput(session, "counts_var_select",
                    choices = c("Choose variable" = "", available_cols),
                    selected = get_first_available_column(values$lazy_data, 
                                                              c("total_counts", "totalUMIs", "nCount_RNA", "nUMI"))[subset_indices])
    
    updateSelectInput(session, "genes_var_select",
                    choices = c("Choose variable" = "", available_cols),
                    selected = get_first_available_column(values$lazy_data, 
                                                              c("n_genes", "nGenes", "genes_count", "nFeature_RNA"))[subset_indices])
  })

  # Enhanced filter preview output
  output$filter_preview_enhanced <- renderText({
    req(values$lazy_data)

    if (!is.null(subset_indices_reactive())) {
      subset_indices <- subset_indices_reactive()$combined
    } else {
      subset_indices <- NULL
    }
    
    # Get total number of cells
    total_cells <- if (is.null(subset_indices)) nrow(values$lazy_data$df) else length(subset_indices)
    
    # Initialize passing filter
    passing <- rep(TRUE, total_cells)
    filters_applied <- c()
    
    # Apply mitochondrial filter
    if (!is.null(input$mito_var_select) && input$mito_var_select != "") {
      mito_data <- values$lazy_data$get_obs_column(values$lazy_data$z, input$mito_var_select, cells = subset_indices)
      if (!is.null(mito_data) && !is.null(input$mito_min) && !is.null(input$mito_max)) {
        passing <- passing & (mito_data >= input$mito_min & mito_data <= input$mito_max)
        filters_applied <- c(filters_applied, 
                            paste0(input$mito_var_select, ": ", input$mito_min, "-", input$mito_max, "%"))
      }
    }
    
    # Apply counts filter
    if (!is.null(input$counts_var_select) && input$counts_var_select != "") {
      counts_data <- values$lazy_data$get_obs_column(values$lazy_data$z, input$counts_var_select, cells = subset_indices)
      if (!is.null(counts_data)) {
        counts_filter <- rep(TRUE, length(counts_data))
        
        if (!is.null(input$counts_min) && !is.na(input$counts_min)) {
          counts_filter <- counts_filter & (counts_data >= input$counts_min)
        }
        if (!is.null(input$counts_max) && !is.na(input$counts_max)) {
          counts_filter <- counts_filter & (counts_data <= input$counts_max)
        }
        
        passing <- passing & counts_filter
        
        filter_text <- paste0(input$counts_var_select, ": ")
        if (!is.null(input$counts_min) && !is.na(input$counts_min)) {
          filter_text <- paste0(filter_text, "≥", input$counts_min)
        }
        if (!is.null(input$counts_max) && !is.na(input$counts_max)) {
          if (!is.null(input$counts_min) && !is.na(input$counts_min)) {
            filter_text <- paste0(filter_text, " & ≤", input$counts_max)
          } else {
            filter_text <- paste0(filter_text, "≤", input$counts_max)
          }
        }
        filters_applied <- c(filters_applied, filter_text)
      }
    }
    
    # Apply genes filter
    if (!is.null(input$genes_var_select) && input$genes_var_select != "") {
      genes_data <- values$lazy_data$get_obs_column(values$lazy_data$z, input$genes_var_select, cells = subset_indices)
      if (!is.null(genes_data)) {
        genes_filter <- rep(TRUE, length(genes_data))
        
        if (!is.null(input$genes_min) && !is.na(input$genes_min)) {
          genes_filter <- genes_filter & (genes_data >= input$genes_min)
        }
        if (!is.null(input$genes_max) && !is.na(input$genes_max)) {
          genes_filter <- genes_filter & (genes_data <= input$genes_max)
        }
        
        passing <- passing & genes_filter
        
        filter_text <- paste0(input$genes_var_select, ": ")
        if (!is.null(input$genes_min) && !is.na(input$genes_min)) {
          filter_text <- paste0(filter_text, "≥", input$genes_min)
        }
        if (!is.null(input$genes_max) && !is.na(input$genes_max)) {
          if (!is.null(input$genes_min) && !is.na(input$genes_min)) {
            filter_text <- paste0(filter_text, " & ≤", input$genes_max)
          } else {
            filter_text <- paste0(filter_text, "≤", input$genes_max)
          }
        }
        filters_applied <- c(filters_applied, filter_text)
      }
    }

      
    passing_cells <- sum(passing, na.rm = TRUE)
    
    if (length(filters_applied) == 0) {
      return("No filters applied")
    }
    
    paste0(passing_cells, " / ", total_cells, " cells (", 
          round(passing_cells / total_cells * 100, 1), "%) pass filters:\n",
          paste(filters_applied, collapse = "\n"))
  })

  

  output$qc_detailed_table <- DT::renderDataTable({
    req(values$lazy_data)

    if (!is.null(subset_indices_reactive())) {
      subset_indices <- subset_indices_reactive()$combined
    } else {
      subset_indices <- NULL
    }
    
    # Get all obs column names
    all_cols <- values$lazy_data$obs_keys
    
    # Load each column lazily
    obs_list <- lapply(all_cols, function(colname) {
      values$lazy_data$get_obs_column(values$lazy_data$z, colname, cells = subset_indices)
    })
    names(obs_list) <- all_cols
    
    # Convert to data.frame
    obs_df <- as.data.frame(obs_list)
    
    DT::datatable(
      obs_df,
      options = list(pageLength = 10, scrollX = TRUE, lengthChange = TRUE),
      filter = 'top',
      rownames = FALSE
    )
  })

  canvas_data_ready <- reactiveVal(FALSE)

  # Reactive value to track report generation status
  report_generation_status <- reactiveVal("idle") # Can be "idle", "capturing", "ready"

  # New button to start the report generation process
  observeEvent(input$generate_report_button, {
    # Show a modal dialog to the user
    showModal(modalDialog(
      id = "report_modal",
      title = "Generating Report...",
      "Preparing plot data. Please wait...",
      footer = NULL,
      easyClose = FALSE
    ))

    # Reset the status and trigger the canvas capture (fallback; can remove if unneeded)
    report_generation_status("capturing")
    print("Triggering canvas capture for report (fallback)...")
    session$sendCustomMessage("captureCanvases", list(debug = "modal_trigger"))
  })

  # Observe canvas data sent from JavaScript (fallback; can remove if unneeded)
  observeEvent(input$canvas_data, {
    # This event will now update the plot_storage and the status
    print("Received canvas_data from JavaScript")
    # plot_storage$canvases <- input$canvas_data
    print("Updated plot_storage$canvases")

    # Set the status to "ready" once the canvas data is received
    # This is the line that was missing or in the wrong place.
    report_generation_status("ready")
    removeModal() # Remove the "Please wait" modal
  })

  # Dynamic UI for the download button
  output$download_button_ui <- renderUI({
    if (report_generation_status() == "ready") {
      downloadButton("download_report", "Download Report")
    } else {
      actionButton("generate_report_button", "Generate Report", class = "btn btn-primary btn-sm")
    }
  })

  output$download_report <- downloadHandler(
    filename = function() {
      paste("analysis_report_", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      withProgress(message = 'Generating Report...', {

        # The downloadHandler now executes only AFTER the canvas data is ready.
        # The waiting logic is no longer needed here.
        # if (is.null(plot_storage$canvases) || length(plot_storage$canvases) == 0) {
        #   showNotification("Error: Canvas images were not captured.", type = "error")
        #   return()
        # }

        # Create temporary directory for plots
        temp_dir <- file.path(tempdir(), "plots")
        dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
        plot_dir <- temp_dir

        # Save static plots as image files
        saved_plots <- list()

        # ... (All your plot-saving code) ...
        if (!is.null(plot_storage$heatmap1)) {
          hm1_file <- file.path(plot_dir, "heatmap1.png")
          png(hm1_file, width = 10, height = 6, units = "in", res = 300)
          draw(plot_storage$heatmap1)
          dev.off()
          saved_plots$heatmap1 <- hm1_file
        }

        if (!is.null(plot_storage$heatmap2)) {
          hm2_file <- file.path(plot_dir, "heatmap2.png")
          png(hm2_file, width = 10, height = 6, units = "in", res = 300)
          draw(plot_storage$heatmap2)
          dev.off()
          saved_plots$heatmap2 <- hm2_file
        }

        for (gene in displayed_genes$genes) {
          gene_plot <- plot_storage[[paste0("gene_", gene)]]
          if (!is.null(gene_plot)) {
            gene_file <- file.path(plot_dir, paste0("gene_", gene, ".png"))
            ggsave(gene_file, plot = gene_plot, width = 8, height = 4, dpi = 300)
            saved_plots[[paste0("gene_", gene)]] <- gene_file
          }
        }


        # qc_main plot saving (if applicable)
        qc_plot_data <- NULL

        if (!is.null(plot_storage$qc_main) && !is.null(qc_main_data())) {
          qc_plot_type <- input$qc_plot_type

          dat <- qc_main_data()
          if (qc_plot_type == "violin") {
            # change dat column names
            colnames(dat) <- c(input$qc_metric, input$qc_group_by)
            print(head(dat))
            qc_plot_data <- list(
              data = as.data.frame(dat)
            )

          } else if (qc_plot_type == "scatter") {

            colnames(dat) <- c("X", "Y", "Color", "Metadata", input$qc_group_by)

            qc_plot_data <- list(
              data = as.data.frame(dat),
              color_by = input$qc_color_by
            )
          }
          saved_plots$qc_main <- TRUE  # Indicate that QC plot data is available
        }
        
        all_scatter_data <- list()

        subset_indices <- isolate({
          if (!is.null(subset_indices_reactive())) {
            subset_indices_reactive()$combined
          } else {
            seq_len(nrow(values$lazy_data$df))
          }
        })

        print("Preparing scatter data for report...")
        if (!is.null(values) && !is.null(values$lazy_data) && !is.null(values$lazy_data$df)) {
          print("Lazy data and df found")
          df <- values$lazy_data$df[subset_indices, , drop = FALSE]
          if ("x" %in% colnames(df) && "y" %in% colnames(df)) {
            
            n_obs <- nrow(df)

            # Annotations for main plot
            annotations <- NULL
            if (!is.null(input$colorBy) && !is.null(values$lazy_data$z)) {
              annotations <- values$lazy_data$get_obs_column(values$lazy_data$z, input$colorBy, cells = subset_indices)
            }
            # Fallback or fix length mismatch
            if (is.null(annotations) || length(annotations) != n_obs) {
              annotations <- rep("default", n_obs)
            }

            main_data <- data.frame(
              x = df$x,
              y = df$y,
              annotations = annotations,
              stringsAsFactors = FALSE
            )
            all_scatter_data[['main']] <- main_data
            print(paste("Prepared main scatter data with", nrow(main_data), "rows"))

            # Gene plots
            gene_keys <- names(plot_storage$canvases)[names(plot_storage$canvases) != "main"]
            if (length(gene_keys) > 0) {
              for (gene in gene_keys) {
                if (!is.null(plot_storage$canvases[[gene]])) {
                  gene_annotations <- extract_gene_data_enhanced(gene, values$lazy_data, only_vector = TRUE)[subset_indices]

                  # Force correct length
                  if (is.null(gene_annotations) || length(gene_annotations) != n_obs) {
                    gene_annotations <- rep(NA, n_obs)
                  }

                  gene_data <- data.frame(
                    x = df$x,
                    y = df$y,
                    annotations = gene_annotations,
                    stringsAsFactors = FALSE
                  )
                  all_scatter_data[[gene]] <- gene_data
                  print(paste("Prepared gene scatter data for", gene, "with", nrow(gene_data), "rows"))
                }
              }
            } else {
              print("No gene plots found in plot_storage$canvases")
            }
          } else {
            print("Warning: 'x' or 'y' columns missing in lazy_data$df")
          }
        } else {
          print("Warning: values$lazy_data or df is null")
        }

        # Reset status for the next run
        report_generation_status("idle")

        # Copy the R Markdown template
        tempReport <- file.path(temp_dir, "report_template.Rmd")
        file.copy("report_template.Rmd", tempReport, overwrite = TRUE)

        # Prepare parameters (now includes scatter_data instead of canvases)
        params <- list(
          app_data = list(
            saved_plots = saved_plots,
            # canvas_info = plot_storage$canvases,  # No longer needed
            violin_genes = displayed_genes$genes,
            timestamp = Sys.time(),
            user_inputs = reactiveValuesToList(input),
            plot_dir = plot_dir,
            scatter_data = all_scatter_data,  # NEW: Pass data for recreating scatterplot in Rmd
            qc_plot_data = qc_plot_data  # NEW: Pass data for recreating QC plot in Rmd
          )
        )

        # Render the document
        print("Starting R Markdown rendering...")
        rmarkdown::render(
          input = tempReport,
          output_file = file,
          params = params,
          envir = new.env(parent = globalenv())
        )
        print("Report rendering complete.")
      })
    }
  )

  plot_storage <- reactiveValues(
    heatmap1 = NULL,
    heatmap2 = NULL,
    gene_main = NULL,
    qc_main = NULL,
    canvases = list()  # Keep for fallback if needed
  )

  output$save_changes_ui <- renderUI({
    req(values$lazy_data)
    actionButton("save_changes", "Save Version", class = "btn btn-secondary btn-sm")
  })

  # Define the R client function for the new endpoint
  check_version_async <- function(zarr_url, version_prefix) {
      VERIFY_URL <- "http://localhost:8080/check_version_status"
      
      response <- httr::GET(
          url = VERIFY_URL,
          query = list(
              zarr_url = zarr_url,
              version_prefix = version_prefix
          )
      )
      
      httr::stop_for_status(response, task = "check version status")
      
      return(httr::content(response, as = "parsed"))
  }

  observeEvent(input$save_changes, {
    req(values$lazy_data)

    # ---- open the modal -------------------------------------------------
    showModal(show_prefix_modal())
  })

  # ---- when the user clicks OK in the modal -----------------------------
  observeEvent(input$prefix_ok, {
    # Grab what the user typed (trim whitespace)
    user_prefix <- trimws(input$version_prefix_input %||% "")

    # If empty → default to today's date
    version_prefix <- if (nzchar(user_prefix)) user_prefix else format(Sys.Date(), "%Y-%m-%d")

    # Close the modal
    removeModal()

    N_TOTAL_CELLS <- nrow(values$lazy_data$df)

    # 🆕 GET CURRENT FILTERED STATE (includes numeric filters)
    subset_info <- subset_indices_reactive()
    subset_indices <- subset_info$combined  # These are the cells passing ALL filters
    
    # Create full cell mask
    full_cell_mask <- rep(FALSE, N_TOTAL_CELLS)
    if (length(subset_indices) > 0) full_cell_mask[subset_indices] <- TRUE

    vis_cat <- unlist(visible_categories()$visibleIndices %||% integer(0))
    sel_points <- unlist(selected_points()$selectedIndices %||% integer(0))
    
    # 🆕 CAPTURE NUMERIC FILTERS STATE
    numeric_filters_state <- input$applyFilters  # This is the JSON string from JS
    
    user_id <- Sys.getenv("SHINYPROXY_USERNAME")
    user_id <- trimws(user_id)
    if (is.na(user_id) || user_id == "") {
      user_id <- "jack"
    }

    print(paste("Using user_id:", user_id))
    version_prefix <- paste0(user_id, "_", version_prefix)

    # ⭐ NEW PIPELINE METADATA
    # Calculate which cells pass each stage
    n_total <- nrow(values$lazy_data$df)

    # Stage 1: Category filter
    category_mask <- rep(TRUE, n_total)
    if (!is.null(vis_cat) && length(vis_cat) > 0 && !is.null(input$colorBy)) {
      annotation_data <- values$lazy_data$get_obs_column(values$lazy_data$z, input$colorBy)
      category_names <- sort(unique(annotation_data))
      visible_category_names <- category_names[vis_cat + 1]
      category_mask <- annotation_data %in% visible_category_names
    }
    n_after_category <- sum(category_mask)

    # Stage 2: Numeric filter
    numeric_mask <- rep(TRUE, n_total)
    if (!is.null(filtered_indices()) && length(filtered_indices()) > 0) {
      numeric_mask <- rep(FALSE, n_total)
      numeric_mask[filtered_indices()] <- TRUE
    }
    n_after_numeric <- sum(category_mask & numeric_mask)

    # Stage 3: Lasso selection (final)
    has_lasso <- !is.null(sel_points) && length(sel_points) > 0

    current_payload_data <- list(
      zarr_url = normalize_url(input$dataset),
      user_id  = user_id,
      version_prefix = version_prefix,
      metadata = list(
        filter_criterion = "Multi-stage filtering pipeline",
        n_total_cells = n_total,
        n_after_category = n_after_category,
        n_after_numeric = n_after_numeric,
        n_filtered_cells = length(subset_indices),
        has_lasso_selection = has_lasso,
        analysis_date = format(Sys.Date(), "%Y-%m-%d"),
        umap_colorby = input$colorBy,
        
        # ⭐ PIPELINE STAGE INFO
        filter_pipeline = list(
          stage1_category = list(
            active = length(vis_cat) < length(unique(annotation_data)),
            annotation = input$colorBy,
            visible_indices = vis_cat
          ),
          stage2_numeric = list(
            active = !is.null(numeric_filters_state) && 
                    numeric_filters_state != "" && 
                    numeric_filters_state != "[]",
            filters = if (!is.null(numeric_filters_state) && 
                        numeric_filters_state != "" && 
                        numeric_filters_state != "[]") {
              jsonlite::fromJSON(numeric_filters_state)
            } else {
              list()
            }
          ),
          stage3_lasso = list(
            active = has_lasso,
            n_selected = length(sel_points)
          )
        )
      ),
      zarr_arrays = list(
        cell_mask = I(full_cell_mask),
        visible_categories = I(vis_cat),
        selected_points = I(sel_points),
        category_mask = I(category_mask),  # ⭐ NEW
        numeric_mask = I(numeric_mask)     # ⭐ NEW
      ),
      credential_id = "default"
    )
    
    id <- showNotification("Saving data and verifying version...",
                          duration = NULL, closeButton = FALSE, type = "default")

    initial_promise <- future({
      resp <- httr::POST(
        url = "http://localhost:8080/add_version",
        body = jsonlite::toJSON(current_payload_data, auto_unbox = TRUE),
        encode = "json",
        httr::content_type_json(),
        httr::timeout(60)
      )
      httr::stop_for_status(resp, task = "save version")
      httr::content(resp, as = "parsed")
    })

    verification_promise <- initial_promise %>%
      promises::then(~{
        version_prefix <- .x$version_prefix
        zarr_url       <- .x$zarr_url

        showNotification(
          paste0("Version ", version_prefix, " created. Starting verification..."),
          duration = 3, type = "message"
        )

        future({ check_version_async(zarr_url, version_prefix) }) %>%
          promises::then(~ list(post = .x, check = .x))
      })

    verification_promise %>%
      promises::then(~{
        check <- .x$check
        removeNotification(id)

        if (check$status == "verified") {
          showNotification(
            paste0("Success! Version ", check$version_prefix, " saved and verified."),
            duration = 5, type = "message"
          )
          version_refresh_trigger(version_refresh_trigger() + 1)
        } else {
          showNotification(
            paste0("Warning: Data saved, but verification failed: ", check$message),
            duration = NULL, type = "warning"
          )
        }
      }) %>%
      promises::catch(~{
        removeNotification(id)
        showNotification(
          paste0("Save Failed: ", .x$message),
          duration = NULL, type = "error"
        )
      })
  })

  output$dge_group_by_ui <- renderUI({
    req(values$lazy_data)
    selectInput("dge_group_by", "Group By:", choices = c("None", categorical_vars()), selected = "celltype")
  })

  output$dge_condition_1_ui <- renderUI({
    req(values$lazy_data)
    selectInput("dge_condition_1", "Condition A:", choices = NULL)
  })

  output$dge_condition_2_ui <- renderUI({
    req(values$lazy_data)
    selectInput("dge_condition_2", "Condition B:", choices = NULL)
  })

  observeEvent(input$dge_group_by, {
    req(values$lazy_data, input$dge_group_by)
    
    if (input$dge_group_by == "None") {
      updateSelectInput(session, "dge_condition_1", choices = NULL)
      updateSelectInput(session, "dge_condition_2", choices = NULL)
    } else {
      group_values <- unique(values$lazy_data$get_obs_column(values$lazy_data$z, input$dge_group_by))
      group_values <- sort(na.omit(group_values))
      
      updateSelectInput(session, "dge_condition_1", choices = group_values, selected = group_values[1])
      updateSelectInput(session, "dge_condition_2", choices = group_values, selected = group_values[2])
    }
  })

  observeEvent(input$run_dge, {
    req(values$lazy_data, input$dge_group_by, input$dge_condition_1, input$dge_condition_2)
    
    # Capture input values
    group_by <- input$dge_group_by
    condition_1 <- input$dge_condition_1
    condition_2 <- input$dge_condition_2
    
    # Validate inputs
    if (condition_1 == condition_2) {
      status("Error: Condition A and Condition B must be different.")
      showNotification("Condition A and Condition B must be different.", 
                      type = "error", duration = 5)
      return()
    }
    
    # Update status
    status("Submitting differential expression job...")
    showNotification("Submitting differential expression job...", 
                    type = "message", duration = 3)
    
    # Start Redis worker
    # startLocalWorkers(n = 1, queue = "scanpy_jobs2")
    # config <- redis_config(host = "127.0.0.1", port = 6379) # Adjust as needed
    # plan(redis, queue = "scanpy_jobs2", config = config)  
    plan(multisession, workers = 1)
    
    zarr_url <- input$dataset
    
    # Create a temporary file for parquet output
    output_parquet <- tempfile(fileext = ".parquet")
    
    # Run the future with PARQUET OUTPUT
    f <- future({
      library(reticulate)
      # use_condaenv("sc_rna_env_python2", required = TRUE)
      use_condaenv("shiny_app_env", conda = "/opt/conda/bin/conda", required = TRUE)
      source_python("./scripts/rank_genes_zarr.py")
      
      message("Running fast DGE analysis...")
      
      parquet_path <- rank_genes_zarr_to_parquet(
        s3_path = zarr_url,
        groupby = group_by,
        groups = list(condition_1, condition_2),
        output_path = output_parquet,
        method = "t-test_overestim_var",
        verbose = TRUE
      )
      
      message("DGE analysis complete, parquet written")
      parquet_path
    }, 
    seed = TRUE, 
    globals = list(
      group_by = group_by,
      condition_1 = condition_1,
      condition_2 = condition_2,
      zarr_url = zarr_url,
      output_parquet = output_parquet
    ))
    
    # Get results with error handling
    res <- tryCatch({
      withProgress(message = "Running fast DGE analysis...", 
                  value = 0, {
        parquet_path <- value(f)
        incProgress(0.5)
        
        message("Reading parquet results...")
        res <- arrow::read_parquet(parquet_path)
        incProgress(0.5)
        res
      })
    }, error = function(e) {
      status(paste("Error in DGE analysis:", e$message))
      showNotification(
        paste("Error in DGE analysis:", e$message),
        type = "error",
        duration = 5
      )
      NULL
    })
    
    # Reset plan
    plan(sequential)
    
    # Store results and update status
    if (!is.null(res)) {
      results(res) # Keep this for temporary view
      
      # 🆕 CREATE COMPREHENSIVE METADATA
      analysis_metadata <- list(
        tool = "DGE",
        group_by = group_by,
        condition_1 = condition_1,
        condition_2 = condition_2,
        analysis_date = format(Sys.Date(), "%Y-%m-%d"),
        analysis_time = format(Sys.time(), "%H:%M:%S"),
        n_genes = nrow(res),
        n_significant = sum(res$pvals_adj < 0.05, na.rm = TRUE),
        pct_significant = round(100 * mean(res$pvals_adj < 0.05, na.rm = TRUE), 1),
        method = "t-test_overestim_var",
        dataset = basename(input$dataset),
        user_id = trimws(Sys.getenv("SHINYPROXY_USERNAME", "jack"))
      )
      
      status(paste0("✓ Analysis completed successfully (~4-7s total) - ", 
                  analysis_metadata$n_significant, " significant genes"))
      
      showNotification(
        sprintf("✅ %s vs %s: %d genes analyzed, %d significant (%.1f%%)",
                condition_1, condition_2, nrow(res), 
                analysis_metadata$n_significant, 
                analysis_metadata$pct_significant),
        type = "message",
        duration = 5
      )
      
      # 🆕 CREATE UNIQUE LABEL AND STORE
      timestamp_str <- format(Sys.time(), "%Y%m%d_%H%M%S")
      version_prefix <- isolate(input$version_prefix %||% "local_default")
      label <- paste0("DGE_", group_by, "_", condition_1, "_vs_", condition_2, "_", timestamp_str)
      
      results_store$all[[label]] <- list(
        tool = "DGE",
        data = res,
        cloud = FALSE,
        label = label,
        version_prefix = version_prefix,
        source = "local",
        metadata = analysis_metadata,  # 🆕 STORE METADATA
        parquet_key = timestamp_str     # 🆕 STORE TIMESTAMP FOR UPLOAD
      )
      
      # 🆕 TRIGGER DROPDOWN UPDATE (This ensures new result appears immediately)
      # Force reactive to re-evaluate
      isolate({
        results_store$trigger <- (results_store$trigger %||% 0) + 1
      })
      
      # Optionally clean up temp file
      unlink(output_parquet)
    } else {
      results(NULL)
      status("Analysis failed")
    }
  })
  
  # Central results store (for history)
  results_store <- reactiveValues(all = list())

  # Track cloud-saved results
  cloud_results <- reactiveValues(labels = character(0))

  # Define fsspec_filesystem for S3 access
  fsspec_filesystem <- function() {
    fsspec <- import("fsspec")
    creds <- jsonlite::fromJSON("credentials.json")$default
    fsspec$filesystem("s3", key = creds$aws_access_key_id, secret = creds$aws_secret_access_key)
  }

  fetch_cloud_result_data <- function(zarr_url, version_prefix, parquet_key) {
    flask_url <- Sys.getenv("FLASK_API_URL", "http://127.0.0.1:8080")
    
    tryCatch({
      # First, try to get the list of results for this version
      resp_list <- httr::GET(
        url = paste0(flask_url, "/list_results"),
        query = list(
          zarr_url = zarr_url,
          version_prefix = version_prefix
        )
      )
      
      if (httr::status_code(resp_list) == 200) {
        content <- httr::content(resp_list, as = "parsed")
        
        # Find the matching result by parquet_key
        if (!is.null(content$results)) {
          matching_result <- Find(function(x) {
            !is.null(x$parquet_key) && x$parquet_key == parquet_key
          }, content$results)
          
          if (!is.null(matching_result) && !is.null(matching_result$data)) {
            return(as.data.frame(matching_result$data))
          }
        }
      }
      
      # If that didn't work, try a direct fetch endpoint
      resp_direct <- httr::GET(
        url = paste0(flask_url, "/get_result"),
        query = list(
          zarr_url = zarr_url,
          version_prefix = version_prefix,
          parquet_key = parquet_key
        )
      )
      
      if (httr::status_code(resp_direct) == 200) {
        content <- httr::content(resp_direct, as = "parsed")
        if (!is.null(content$data)) {
          return(as.data.frame(content$data))
        }
      }
      
      return(NULL)
    }, error = function(e) {
      message("Error fetching cloud result data: ", e$message)
      return(NULL)
    })
  }

  # Helper to fetch cloud results
  fetch_cloud_results <- function(zarr_url, flask_url = Sys.getenv("FLASK_API_URL", "http://127.0.0.1:8080")) {
    resp <- httr::GET(
      url = paste0(flask_url, "/list_results"),
      query = list(zarr_url = zarr_url)
    )
    if (httr::status_code(resp) != 200) {
      message("Failed to fetch cloud results. Status: ", httr::status_code(resp))
      return(list())
    }
    parsed <- httr::content(resp, as = "parsed")
    if (!is.list(parsed) || is.null(parsed$results)) {
      message("Invalid response format: no 'results' field")
      return(list())
    }
    results <- parsed$results
    message("Parsed results: ", capture.output(str(results)))
    
    cloud_results_list <- list()
    for (r in results) {
      if (!is.list(r) || is.null(r$version_prefix) || !is.list(r$parquet_files)) {
        message("Skipping invalid version entry: ", capture.output(str(r)))
        next
      }
      for (p in r$parquet_files) {
        if (!is.character(p) || !grepl("\\.parquet$", p)) {
          message("Skipping invalid parquet file: ", capture.output(str(p)))
          next
        }
        label <- paste0("DGE_", r$version_prefix, "_", sub("\\.parquet$", "", p))
        cloud_results_list[[length(cloud_results_list) + 1]] <- list(
          label = label,
          version_prefix = r$version_prefix,
          parquet_file = p,
          zarr_url = zarr_url,
          tool = "DGE",
          cloud = TRUE
        )
      }
    }
    message("Cloud results list: ", capture.output(str(cloud_results_list)))
    
    return(cloud_results_list)
  }

  extract_comparison_info <- function(result_obj, zarr_url = NULL, credential_id = "default") {
    tryCatch({
      # PRIORITY 1: Check stored metadata (for LOCAL results)
      if (!is.null(result_obj$metadata)) {
        return(list(
          group1 = result_obj$metadata$condition_1 %||% "Group1",
          group2 = result_obj$metadata$condition_2 %||% "Group2",
          date = as.Date(result_obj$metadata$analysis_date %||% Sys.Date()),
          n_genes = result_obj$metadata$n_genes %||% nrow(result_obj$data),
          n_significant = result_obj$metadata$n_significant %||% NA,
          pct_significant = result_obj$metadata$pct_significant %||% NA,
          groupby = result_obj$metadata$group_by %||% "unknown"
        ))
      }
      
      # PRIORITY 2: Fetch from cloud (for CLOUD results)
      if (result_obj$cloud && !is.null(result_obj$version_prefix) && !is.null(zarr_url)) {
        flask_url <- Sys.getenv("FLASK_API_URL", "http://127.0.0.1:8080")
        
        # Extract parquet key from label or parquet_file
        parquet_key <- if (!is.null(result_obj$parquet_key)) {
          result_obj$parquet_key
        } else if (!is.null(result_obj$parquet_file)) {
          sub("\\.parquet$", "", result_obj$parquet_file)
        } else {
          sub("^.*_(\\d{8}_\\d{6})$", "\\1", result_obj$label)
        }
        
        resp <- httr::GET(
          url = paste0(flask_url, "/get_result_metadata"),
          query = list(
            zarr_url = zarr_url,
            version_prefix = result_obj$version_prefix,
            parquet_key = parquet_key
          ),
          httr::timeout(10)
        )
        
        if (httr::status_code(resp) == 200) {
          metadata <- httr::content(resp, as = "parsed")
          return(list(
            group1 = metadata$condition_1 %||% "Group1",
            group2 = metadata$condition_2 %||% "Group2",
            date = as.Date(metadata$analysis_date %||% Sys.Date()),
            n_genes = metadata$n_genes %||% NA,
            n_significant = metadata$n_significant %||% NA,
            pct_significant = metadata$pct_significant %||% NA,
            groupby = metadata$group_by %||% "unknown"
          ))
        }
      }
      
      NULL
    }, error = function(e) {
      message("Could not extract comparison info: ", e$message)
      NULL
    })
  }

  # 🆕 ENHANCED available_dge_results() with better display names
  available_dge_results <- reactive({
    req(input$dataset, input$version_prefix)
    
    # Add dependency on trigger to force re-evaluation when new results added
    results_store$trigger
    
    zarr_url <- normalize_url(input$dataset)
    version_prefix <- input$version_prefix
    
    message("\n🔄 Refreshing available DGE results...")
    message("   Dataset: ", basename(input$dataset))
    message("   Version prefix: ", version_prefix)
    
    # 🆕 CHECK IF DATASET IS CLOUD-BASED
    is_cloud_dataset <- grepl("^s3://|^https?://", zarr_url)
    message("   Is cloud dataset: ", is_cloud_dataset)
    
    # 1. Fetch cloud results ONLY if dataset is cloud-based
    cloud_filtered <- list()
    if (is_cloud_dataset) {
      cloud_results_list <- tryCatch({
        fetch_cloud_results(zarr_url)
      }, error = function(e) {
        message("   ⚠️ Failed to fetch cloud results: ", e$message)
        list()
      })
      message("   Fetched ", length(cloud_results_list), " cloud results from API")
      
      # 2. Filter cloud results by version_prefix
      if (length(cloud_results_list) > 0) {
        cloud_filtered <- Filter(function(x) {
          identical(x$version_prefix, version_prefix)
        }, cloud_results_list)
      }
      message("   After version filtering: ", length(cloud_filtered), " cloud results")
    } else {
      message("   Skipping cloud fetch for local dataset")
    }
    
    # 3. Get local results (all versions)
    local_all <- Filter(function(x) x$tool == "DGE", results_store$all)
    message("   Found ", length(local_all), " local results in store")
    
    # 4. Filter local results by version_prefix AND ensure they're actually local
    local_filtered <- list()
    if (length(local_all) > 0) {
      local_filtered <- Filter(function(x) {
        # 🆕 MUST be local source AND match version (or have no version)
        is_local <- is.null(x$source) || x$source == "local"
        matches_version <- is.null(x$version_prefix) || identical(x$version_prefix, version_prefix)
        is_local && matches_version
      }, local_all)
    }
    message("   After filtering: ", length(local_filtered), " local results")
    
    # 5. Build final list WITH ENHANCED DISPLAY NAMES
    results_list <- list()
    
    # Helper function to create display name
    create_display_name <- function(comparison_info, compact = TRUE) {
      if (is.null(comparison_info)) {
        return(NULL)
      }
      
      max_length <- if (compact) 10 else 15
      
      group1 <- if (nchar(comparison_info$group1) > max_length) {
        paste0(substr(comparison_info$group1, 1, max_length - 2), "...")
      } else {
        comparison_info$group1
      }
      
      group2 <- if (nchar(comparison_info$group2) > max_length) {
        paste0(substr(comparison_info$group2, 1, max_length - 2), "...")
      } else {
        comparison_info$group2
      }
      
      if (compact) {
        groupby_part <- if (!is.null(comparison_info$groupby) && comparison_info$groupby != "unknown") {
          paste0(substr(comparison_info$groupby, 1, 8), ": ")
        } else {
          ""
        }
        
        display_name <- sprintf(
          "%s%s vs %s (%s)",
          groupby_part,
          group1,
          group2,
          format(comparison_info$date, "%m/%d")
        )
        
        if (!is.na(comparison_info$n_genes)) {
          genes_k <- if (comparison_info$n_genes >= 1000) {
            paste0(round(comparison_info$n_genes / 1000, 1), "K")
          } else {
            as.character(comparison_info$n_genes)
          }
          display_name <- paste0(display_name, " • ", genes_k)
        }
        
        if (!is.na(comparison_info$n_significant)) {
          display_name <- paste0(display_name, " • ", comparison_info$n_significant, "↑")
        }
        
      } else {
        groupby_prefix <- if (!is.null(comparison_info$groupby) && comparison_info$groupby != "unknown") {
          paste0("[", comparison_info$groupby, "] ")
        } else {
          ""
        }
        
        display_name <- sprintf(
          "%s%s vs %s (%s)",
          groupby_prefix,
          group1,
          group2,
          format(comparison_info$date, "%b %d")
        )
        
        if (!is.na(comparison_info$n_genes)) {
          display_name <- paste0(display_name, sprintf(" • %s genes", 
                                                      format(comparison_info$n_genes, big.mark = ",")))
        }
        
        if (!is.na(comparison_info$n_significant)) {
          display_name <- paste0(display_name, sprintf(" • %d sig", comparison_info$n_significant))
        }
      }
      
      display_name
    }
    
    # Cloud results (only if cloud dataset)
    if (is_cloud_dataset && length(cloud_filtered) > 0) {
      for (res in cloud_filtered) {
        label <- res$label
        comparison_info <- extract_comparison_info(res, zarr_url, credential_id = "default")
        display_name <- create_display_name(comparison_info, compact = TRUE) %||% label
        
        results_list[[label]] <- c(res, list(
          source = "cloud",
          display_name = display_name
        ))
      }
    }
    
    # Local results (should not duplicate cloud results anymore)
    cloud_labels <- names(results_list)
    if (length(local_filtered) > 0) {
      for (label in names(local_filtered)) {
        # Skip if this label is already in cloud results
        if (label %in% cloud_labels) {
          message("   ⚠️ Skipping duplicate: ", label, " (already in cloud)")
          next
        }
        
        res <- local_filtered[[label]]
        comparison_info <- extract_comparison_info(res, zarr_url = NULL)
        display_name <- create_display_name(comparison_info, compact = TRUE) %||% label
        
        results_list[[label]] <- c(res, list(
          source = "local",
          display_name = display_name
        ))
      }
    }
    
    # 🆕 FIX: Safely count cloud vs local
    if (length(results_list) > 0) {
      # Use unlist() to ensure we get a logical vector
      n_cloud <- sum(unlist(lapply(results_list, function(x) {
        identical(x$source, "cloud")
      })))
      n_local <- sum(unlist(lapply(results_list, function(x) {
        identical(x$source, "local")
      })))
    } else {
      n_cloud <- 0
      n_local <- 0
    }
    
    message("✅ Final results: ", length(results_list), " total")
    message("   Cloud: ", n_cloud)
    message("   Local: ", n_local)
    
    results_list
  })

  

  # ===================================================================
  # UPDATE DROPDOWN whenever dataset, version_prefix, or new results appear
  # ===================================================================
  # 🆕 UPDATE DROPDOWN with enhanced display names and no duplicates
  observe({
    results_list <- available_dge_results()
    
    if (length(results_list) == 0) {
      updateSelectInput(session, "previous_result", choices = character(0), selected = character(0))
      return()
    }
    
    # Sort by timestamp (newest first)
    labels <- names(results_list)
    timestamps <- sapply(labels, function(label) {
      # Extract timestamp: DGE_..._20251109_151317 -> 20251109151317
      ts_str <- sub("^.*_(\\d{8}_\\d{6})$", "\\1", label)
      if (grepl("^\\d{8}_\\d{6}$", ts_str)) {
        as.numeric(gsub("_", "", ts_str))
      } else {
        0  # Put malformed labels at the end
      }
    })
    
    order_idx <- order(timestamps, decreasing = TRUE)
    labels <- labels[order_idx]
    
    # Get display names and sources
    display_names <- sapply(labels, function(label) {
      results_list[[label]]$display_name %||% label
    })
    
    sources <- sapply(labels, function(label) {
      results_list[[label]]$source
    })
    
    # 🆕 Add source prefix with better emojis and spacing
    display_labels <- ifelse(
      sources == "cloud",
      paste0("☁️  ", display_names),  # Extra space for readability
      paste0("💻  ", display_names)
    )
    
    # Preserve current selection if it still exists
    current_selection <- input$previous_result
    
    # 🆕 If current selection was just uploaded, it's now cloud - keep it selected
    if (!is.null(current_selection)) {
      # Strip any emoji prefix for comparison
      clean_selection <- sub("^☁️  |^💻  ", "", current_selection)
      if (clean_selection %in% labels) {
        current_selection <- clean_selection
      } else if (!current_selection %in% labels) {
        current_selection <- character(0)
      }
    } else {
      current_selection <- character(0)
    }
    
    updateSelectInput(
      session,
      "previous_result",
      choices = setNames(labels, display_labels),
      selected = current_selection
    )
    
    cat("📋 Dropdown updated with", length(labels), "results\n")
    cat("   Sources: Cloud =", sum(sources == "cloud"), ", Local =", sum(sources == "local"), "\n")
  }) |> 
    bindEvent(
      input$dataset,
      input$version_prefix,
      available_dge_results(),
      ignoreNULL = FALSE
    )

  # Save to cloud
  observeEvent(input$save_to_cloud, {
    tool <- input$result_tool
    req(tool == "DGE", results(), input$version_prefix)
    
    # CHECK IF DATASET IS CLOUD-BASED
    zarr_url <- normalize_url(input$dataset)
    is_cloud_dataset <- grepl("^s3://|^https?://", zarr_url)
    
    if (!is_cloud_dataset) {
      showNotification(
        "⚠️ Cannot save to cloud: Dataset is local. Only cloud datasets support versioning.",
        type = "warning",
        duration = 5
      )
      return()
    }
    
    user_id <- trimws(Sys.getenv("SHINYPROXY_USERNAME"))
    if (is.na(user_id) || user_id == "") {
      user_id <- "jack"
    }
    
    version_prefix <- input$version_prefix
    
    # 🆕 HANDLE "ORIGINAL" VERSION - MUST CREATE A NEW VERSION
    if (version_prefix == "original") {
      showModal(modalDialog(
        title = "Create Version First",
        div(
          p("You're currently on the 'original' version. To save analysis results to the cloud, you need to create a version first."),
          p(strong("This version will capture your current cell selection/filtering state.")),
          hr(),
          textInput("new_version_name", "Version Name:", 
                  value = paste0("analysis_", format(Sys.Date(), "%Y%m%d")),
                  placeholder = "e.g., filtered_bcells")
        ),
        footer = tagList(
          modalButton("Cancel"),
          actionButton("create_and_save", "Create Version & Save Results", class = "btn-primary")
        ),
        easyClose = FALSE
      ))
      return()
    }
    
    # For non-original versions, continue with normal save flow
    selected_label <- sub("^☁️  |^💻  ", "", input$previous_result)
    
    # Get the result object to access metadata
    result_obj <- results_store$all[[selected_label]]
    req(result_obj)
    
    parquet_key <- result_obj$parquet_key %||% sub("^.*_(\\d{8}_\\d{6})$", "\\1", selected_label)
    
    if (parquet_key == selected_label || !grepl("^\\d{8}_\\d{6}$", parquet_key)) {
      showNotification("Error: Could not parse timestamp from result label.", type = "error")
      return()
    }
    
    # Validate version_prefix
    flask_url <- Sys.getenv("FLASK_API_URL", "http://127.0.0.1:8080")
  
    # Write parquet locally first
    local_parquet <- tempfile(fileext = ".parquet")
    arrow::write_parquet(results(), local_parquet)
    
    # Verify file size
    file_size <- file.info(local_parquet)$size / 1024  # KB
    cat("📦 Local parquet size:", round(file_size, 2), "KB\n")
    
    # Prepare metadata to send with upload
    metadata_to_send <- result_obj$metadata %||% list(
      tool = "DGE",
      analysis_date = format(Sys.Date(), "%Y-%m-%d"),
      n_genes = nrow(results()),
      user_id = user_id
    )
    
    id <- showNotification("Uploading DGE results...", duration = NULL, closeButton = FALSE, type = "default")
    
    # Upload via multipart with metadata
    future({
      resp <- httr::POST(
        url = paste0(flask_url, "/upload_parquet"),
        body = list(
          zarr_url = zarr_url,
          user_id = user_id,
          version_prefix = version_prefix,
          parquet_key = parquet_key,
          metadata = jsonlite::toJSON(metadata_to_send, auto_unbox = TRUE),
          file = httr::upload_file(local_parquet)
        ),
        encode = "multipart",
        httr::timeout(300)
      )
      httr::stop_for_status(resp, task = "upload parquet")
      httr::content(resp, as = "parsed")
    }) %...>%
      (function(response) {
        removeNotification(id)
        
        # Verify uploaded size
        if (!is.null(response$uploaded_size_kb)) {
          cat("☁️ Cloud parquet size:", response$uploaded_size_kb, "KB\n")
          if (abs(response$uploaded_size_kb - file_size) > 1) {
            showNotification("⚠️ Size mismatch detected!", type = "warning", duration = 5)
          }
        }
        
        # REMOVE from local store (it's now in cloud)
        if (selected_label %in% names(results_store$all)) {
          results_store$all[[selected_label]] <- NULL
          cat("✅ Removed", selected_label, "from local store (now in cloud)\n")
        }
        
        showNotification(
          sprintf("✅ %s saved to cloud (version: %s)",
                  result_obj$metadata$condition_1 %||% "Results",
                  version_prefix),
          duration = 5, type = "message"
        )
        
        # Force dropdown refresh
        isolate({
          results_store$trigger <- (results_store$trigger %||% 0) + 1
        })
      }) %...!%
      (function(error) {
        removeNotification(id)
        showNotification(
          paste0("❌ Upload Failed: ", error$message),
          duration = NULL, type = "error"
        )
      })
    
    # Clean up temp file
    later::later(~ unlink(local_parquet), delay = 5)
  })

  # 🆕 HANDLE CREATE VERSION THEN SAVE
  observeEvent(input$create_and_save, {
    req(input$new_version_name)
    
    new_version_name <- trimws(input$new_version_name)
    
    if (new_version_name == "" || new_version_name == "original") {
      showNotification("Please provide a valid version name (not 'original')", 
                      type = "error", duration = 3)
      return()
    }
    
    removeModal()
    
    # Step 1: Create the version
    N_TOTAL_CELLS <- nrow(values$lazy_data$df)
    
    subset_info <- subset_indices_reactive()
    subset_indices <- subset_info$combined
    
    full_cell_mask <- rep(FALSE, N_TOTAL_CELLS)
    if (length(subset_indices) > 0) full_cell_mask[subset_indices] <- TRUE
    
    vis_cat <- unlist(visible_categories()$visibleIndices %||% integer(0))
    sel_points <- unlist(selected_points()$selectedIndices %||% integer(0))
    numeric_filters_state <- input$applyFilters
    
    user_id <- trimws(Sys.getenv("SHINYPROXY_USERNAME"))
    if (is.na(user_id) || user_id == "") user_id <- "jack"
    
    version_prefix <- paste0(user_id, "_", new_version_name)
    zarr_url <- normalize_url(input$dataset)
    
    create_payload <- list(
      zarr_url = zarr_url,
      user_id = user_id,
      version_prefix = version_prefix,
      metadata = list(
        filter_criterion = "Created for analysis results",
        n_filtered_cells = length(subset_indices),
        analysis_date = format(Sys.Date(), "%Y-%m-%d"),
        umap_colorby = input$colorBy,
        has_numeric_filters = !is.null(numeric_filters_state) && 
                            numeric_filters_state != "" && 
                            numeric_filters_state != "[]",
        numeric_filters = if (!is.null(numeric_filters_state) && 
                            numeric_filters_state != "" && 
                            numeric_filters_state != "[]") {
          jsonlite::fromJSON(numeric_filters_state)
        } else {
          list()
        }
      ),
      zarr_arrays = list(
        cell_mask = I(full_cell_mask),
        visible_categories = I(vis_cat),
        selected_points = I(sel_points)
      ),
      credential_id = "default"
    )
    
    id <- showNotification("Creating version and saving results...",
                          duration = NULL, closeButton = FALSE, type = "default")
    
    # Step 2: Create version, then save results
    future({
      # Create version
      resp1 <- httr::POST(
        url = "http://localhost:8080/add_version",
        body = jsonlite::toJSON(create_payload, auto_unbox = TRUE),
        encode = "json",
        httr::content_type_json(),
        httr::timeout(60)
      )
      httr::stop_for_status(resp1, task = "create version")
      
      version_info <- httr::content(resp1, as = "parsed")
      
      # Now upload results to this new version
      selected_label <- sub("^☁️  |^💻  ", "", input$previous_result)
      result_obj <- results_store$all[[selected_label]]
      parquet_key <- result_obj$parquet_key %||% sub("^.*_(\\d{8}_\\d{6})$", "\\1", selected_label)
      
      local_parquet <- tempfile(fileext = ".parquet")
      arrow::write_parquet(results(), local_parquet)
      
      metadata_to_send <- result_obj$metadata %||% list(
        tool = "DGE",
        analysis_date = format(Sys.Date(), "%Y-%m-%d"),
        n_genes = nrow(results()),
        user_id = user_id
      )
      
      resp2 <- httr::POST(
        url = "http://localhost:8080/upload_parquet",
        body = list(
          zarr_url = zarr_url,
          user_id = user_id,
          version_prefix = version_prefix,
          parquet_key = parquet_key,
          metadata = jsonlite::toJSON(metadata_to_send, auto_unbox = TRUE),
          file = httr::upload_file(local_parquet)
        ),
        encode = "multipart",
        httr::timeout(300)
      )
      httr::stop_for_status(resp2, task = "upload parquet")
      
      unlink(local_parquet)
      
      list(
        version = httr::content(resp1, as = "parsed"),
        upload = httr::content(resp2, as = "parsed")
      )
    }) %...>%
      (function(response) {
        removeNotification(id)
        
        # Remove from local store
        selected_label <- sub("^☁️  |^💻  ", "", input$previous_result)
        if (selected_label %in% names(results_store$all)) {
          results_store$all[[selected_label]] <- NULL
        }
        
        showNotification(
          sprintf("✅ Version '%s' created and results saved!", new_version_name),
          duration = 5, type = "message"
        )
        
        # Update version dropdown to include new version
        version_refresh_trigger(version_refresh_trigger() + 1)
        
        # Select the new version
        updateSelectInput(session, "version_prefix", selected = version_prefix)
        
        # Force results dropdown refresh
        isolate({
          results_store$trigger <- (results_store$trigger %||% 0) + 1
        })
      }) %...!%
      (function(error) {
        removeNotification(id)
        showNotification(
          paste0("❌ Failed: ", error$message),
          duration = NULL, type = "error"
        )
      })
  })

  output$analysis_display <- renderUI({
    req(input$result_tool)
    tagList(
      h4(paste("Results for", input$result_tool)),
      
      # conditionalPanel(
      #   condition = "output.has_results == true",
      #   div(
      #     style = "margin-bottom: 10px;",
      #     downloadButton("download_result", "⬇️ Download CSV", class = "btn-secondary btn-sm"),
      #     actionButton("save_to_cloud", "☁️ Save to Cloud", class = "btn-success btn-sm")
      #   )
      # ),
      
      # VOLCANO PLOT + TABLE FOR DGE (OPTIMIZED)
      conditionalPanel(
        condition = "input.result_tool == 'DGE' && output.has_results == true",
        fluidRow(
          # VOLCANO PLOT (left - 7 columns for more space)
          column(7,
            h5("🌋 Volcano Plot"),
            p("Click or brush to select genes | Click background to deselect", 
              style = "font-size: 0.85em; color: #666; margin-bottom: 10px;"),
            my_scatterplotOutput("volcano_plot", height = "700px", width = "100%")  # Taller for better view
          ),
          
          # GENES TABLE (right - 5 columns)
          column(5,
            h5("📋 Selected Genes"),
            div(
              style = "margin-bottom: 10px;",
              radioButtons("volcano_table_mode", NULL, inline = TRUE,
                          choices = list(
                            "Top 1000 + Selected" = "all",
                            "Selected only" = "selected"
                          ),
                          selected = "all"),
              # Info text
              uiOutput("volcano_table_info")
            ),
            DT::dataTableOutput("volcano_genes_table", height = "650px")
          )
        ),
        hr()
      ),
      
      # DT::dataTableOutput("results_table")
    )
  })

  output$save_to_cloud_button_ui <- renderUI({
    req(input$dataset)
    
    zarr_url <- normalize_url(input$dataset)
    is_cloud_dataset <- grepl("^s3://|^https?://", zarr_url)
    
    if (is_cloud_dataset) {
      div(
        style = "margin-top: 25px;",
        actionButton("save_to_cloud", "☁️ Save", 
                    class = "btn-success btn-sm",
                    style = "width: 100%;")
      )
    } else {
      div(
        style = "margin-top: 25px;",
        tags$span(
          class = "label label-default",
          style = "display: block; padding: 6px; text-align: center; font-size: 11px;",
          "💻 Local dataset"
        )
      )
    }
  })

  
  # ===================================================================
  # REACTIVE: Currently selected DGE result (local OR cloud)
  # ===================================================================
  current_dge_result <- reactive({
    req(input$previous_result, input$version_prefix)
    
    selected_label <- sub("^(Cloud|Local) ", "", input$previous_result)
    message("Loading public result: ", selected_label)
    
    # LOCAL FIRST
    local <- results_store$all[[selected_label]]
    if (!is.null(local) && local$tool == "DGE") {
      message("Local cache hit: ", nrow(local$data), " rows")
      return(local)
    }

    # === PUBLIC HTTPS DIRECT READ ===
    base_url <- input$dataset  # ← ALREADY PUBLIC HTTPS URL
    # Example: https://scrnaseq-browser.s3.us-east-2.amazonaws.com/corrected_....zarr
    
    # Extract timestamp from label: DGE_jack_test_20251109_133651 → 20251109_133651
    timestamp <- sub("^.*_(\\d{8}_\\d{6})$", "\\1", selected_label)
    if (timestamp == selected_label) {
      showNotification("Cannot parse timestamp", type = "error")
      return(NULL)
    }
    
    parquet_file <- paste0(timestamp, ".parquet")
    public_url <- file.path(
      base_url, 
      "versions", 
      input$version_prefix, 
      "results", 
      parquet_file
    )
    
    message("PUBLIC DIRECT LOAD: ", public_url)
    
    pt <- proc.time()
    df <- tryCatch({
      arrow::read_parquet(public_url)
    }, error = function(e) {
      showNotification(paste("Load failed:", e$message), type = "error")
      message("ERROR: ", e$message)
      NULL
    })
    
    req(df, nrow(df) > 0)
    
    message("PUBLIC LOAD SUCCESS: ", nrow(df), " rows in ", 
      round((proc.time() - pt)[3], 3), "s")
    
    list(
      tool = "DGE",
      data = df,
      cloud = TRUE,
      label = selected_label,
      version_prefix = input$version_prefix,
      source = "public"
    )
  })

  # ===================================================================
  # VOLCANO PLOT - Now works for BOTH local and cloud
  # ===================================================================
  output$volcano_plot <- renderMy_scatterplot({
    req(input$result_tool == "DGE")
    
    result_data <- current_dge_result()
    req(result_data, result_data$data)
    
    df <- result_data$data
    
    # === [REST OF YOUR VOLCANO CODE - UNCHANGED] ===
    # Gene naming
    if ("names" %in% colnames(df)) {
      df$gene <- df$names
    } else if ("gene" %in% colnames(df)) {
      df$gene <- df$gene
    } else if (!is.null(rownames(df)) && !all(rownames(df) == as.character(1:nrow(df)))) {
      df$gene <- rownames(df)
    } else {
      df$gene <- paste0("Gene_", seq_len(nrow(df)))
    }
    
    df$neg_log10_padj <- -log10(pmax(df$pvals_adj, 1e-300))
    max_finite <- max(df$neg_log10_padj[is.finite(df$neg_log10_padj)], na.rm = TRUE)
    df$neg_log10_padj[is.infinite(df$neg_log10_padj)] <- max_finite * 1.2
    
    df$significant <- ifelse(
      df$pvals_adj < 0.05 & abs(df$logfoldchanges) > 1, 
      "Significant", "Not Significant"
    )
    
    valid_idx <- !is.na(df$logfoldchanges) & !is.na(df$neg_log10_padj)
    df <- df[valid_idx, ]
    
    x_axis_limit <- 10
    y_axis_limit <- 50
    
    df$logfoldchanges_capped <- pmax(-x_axis_limit, pmin(x_axis_limit, df$logfoldchanges))
    df$neg_log10_padj_capped <- pmin(y_axis_limit, df$neg_log10_padj)
    df$is_capped <- (abs(df$logfoldchanges) > x_axis_limit) | (df$neg_log10_padj > y_axis_limit)
    
    # Store for table
    current_checksum <- digest::digest(df$gene)
    if (is.null(volcano_table_data()) || digest::digest(volcano_table_data()$gene) != current_checksum) {
      volcano_table_data(df)
      volcano_selected_indices(integer(0))
    }
    
    df <- volcano_table_data()
    req(df)
    selected <- volcano_selected_indices()
    
    categories <- rep("Not Significant", nrow(df))
    categories[df$significant == "Significant"] <- "Significant"
    capped_idx <- which(df$is_capped)
    if (length(capped_idx) > 0) {
      categories[capped_idx] <- paste0(categories[capped_idx], "_Capped")
    }
    if (length(selected) > 0) {
      for (idx in selected) {
        base <- if (df$significant[idx] == "Significant") "Selected_Significant" else "Selected"
        categories[idx] <- if (df$is_capped[idx]) paste0(base, "_Capped") else base
      }
    }
    
    color_palette <- c(
      "Not Significant" = "#CCCCCC",
      "Not Significant_Capped" = "#999999",
      "Significant" = "#E74C3C",
      "Significant_Capped" = "#C0392B",
      "Selected" = "#F1C40F",
      "Selected_Capped" = "#F39C12",
      "Selected_Significant" = "#D35400",
      "Selected_Significant_Capped" = "#A04000"
    )
    
    point_size <- if (nrow(df) > 10000) 2 else if (nrow(df) > 5000) 3 else 4
    
    data_version <- digest::digest(list(df$gene, df$logfoldchanges_capped, df$neg_log10_padj_capped, categories))
    unique_id <- paste0("volcano_", substr(data_version, 1, 8))
    
    my_scatterplot(
      x = df$logfoldchanges_capped,
      y = df$neg_log10_padj_capped,
      colorBy = categories,
      custom_palette = color_palette,
      xlab = "log2 Fold Change",
      ylab = "-log10(adjusted p-value)",
      xrange = c(-10, 10),
      yrange = c(0, 50),
      showAxes = TRUE,
      showTooltip = TRUE,
      opacity = 0.6,
      size = point_size,
      backgroundColor = "#FFFFFF",
      legend_title = "Significance",
      enableDownload = TRUE,
      gene_names = df$gene,
      plotId = "volcano_plot",
      elementId = unique_id,
      dataVersion = data_version
    )
  })

  # --- INFO TEXT ABOUT TABLE CONTENT ---
  output$volcano_table_info <- renderUI({
    req(volcano_table_data())
    
    n_total <- nrow(volcano_table_data())
    n_selected <- length(volcano_selected_indices())
    
    if (input$volcano_table_mode == "all") {
      p(paste0("Showing top ", min(1000, n_total), " genes (including ", 
              n_selected, " selected)."),
        style = "font-size: 0.8em; color: #666; margin: 5px 0;")
    } else { # mode == "selected"
      p(paste0("Showing ", n_selected, " selected genes."),
        style = "font-size: 0.8em; color: #666; margin: 5px 0;")
    }
  })

  # --- INITIAL TABLE RENDER ---
  output$volcano_genes_table <- DT::renderDataTable({
    req(volcano_table_data())
    
    df <- volcano_table_data()
    result <- prepare_volcano_table(df, integer(0), "all", max_rows = 1000)
    
    cat("📋 Table: rendering", nrow(result$display), "rows\n")
    
    DT::datatable(
      result$display,
      selection = list(mode = 'multiple', selected = result$selected_rows),
      options = list(
        pageLength = 50,
        scrollY = "580px", 
        scrollCollapse = TRUE,
        dom = 'ftp', 
        ordering = TRUE, 
        order = list(list(1, 'desc'), list(4, 'desc')),
        searching = TRUE,
        deferRender = TRUE,
        # PROPER way to hide column - use columnDefs
        columnDefs = list(
          list(targets = 1, visible = FALSE)  # Hide Sel column (0-indexed, so column 1)
        ),
        rowCallback = DT::JS(
          "function(row, data, index) {",
          "  // data[5] is 'Status' column",
          "  if (data[5] === 'Significant') {",
          "    $(row).addClass('significant-row');",
          "  }",
          "  // data[1] is the 'Sel' column (hidden but still in data)",
          "  if (data[1] === true) {",
          "    $(row).addClass('selected-row-highlight');",
          "  }",
          "}"
        )
      ),
      rownames = FALSE,
      filter = 'none'
    ) %>%
      DT::formatStyle(
        'Status',
        color = DT::styleEqual('Significant', '#d32f2f'),
        fontWeight = DT::styleEqual('Significant', 'bold')
      )
  })

  # --- OBSERVER: Update table when selection or mode changes (OPTIMIZED) ---
  observe({
    req(volcano_table_data())
    req(input$previous_result)  # Add this line to watch for result changes
    
    selected <- volcano_selected_indices()
    mode <- input$volcano_table_mode
    df <- volcano_table_data()
    
    # Generate the data frame for the current view
    result <- prepare_volcano_table(df, selected, mode, max_rows = 1000)
    
    cat("📋 Table update:", nrow(result$display), "rows,", length(selected), "selected\n")
    
    # Use replaceData for speed
    proxy <- DT::dataTableProxy("volcano_genes_table")
    # Reset paging if mode changed, otherwise keep it
    reset_paging <- !is.null(isolate(input$volcano_table_mode)) && 
                    isolate(input$volcano_table_mode) != mode
                    
    proxy %>% DT::replaceData(result$display, resetPaging = reset_paging, rownames = FALSE)
    
    # Update row selection
    if (length(result$selected_rows) > 0) {
      shinyjs::delay(50, {
        proxy %>% DT::selectRows(result$selected_rows)
      })
    } else {
      shinyjs::delay(50, {
        proxy %>% DT::selectRows(NULL)
      })
    }
  }) %>% bindEvent(volcano_selected_indices(), input$volcano_table_mode, input$previous_result, ignoreNULL = FALSE)


  # --- OBSERVER: Handle Plot Selection ---
  observeEvent(input$volcano_plot_selected, {
    if (volcano_updating_from_table()) return()
    
    volcano_updating_from_plot(TRUE)
    on.exit(volcano_updating_from_plot(FALSE))
    
    indices_raw <- input$volcano_plot_selected$indices
    
    # Handle list or atomic vector from reglScatterplot
    if (is.list(indices_raw)) {
      indices_0based <- unlist(indices_raw, use.names = FALSE)
    } else {
      indices_0based <- as.numeric(indices_raw)
    }
    
    if (length(indices_0based) > 0) {
      # Convert to 1-based indices for R data frames
      full_indices_1based <- sort(unique(indices_0based + 1))
      volcano_selected_indices(full_indices_1based)
    } else {
      volcano_selected_indices(integer(0))
    }
  })

  # --- OBSERVER: Handle Table Row Selection ---
  observeEvent(input$volcano_genes_table_rows_selected, {
    # If the selection update is coming *from the plot*, ignore it here
    if (volcano_updating_from_plot()) return()
    
    volcano_updating_from_table(TRUE)
    on.exit(volcano_updating_from_table(FALSE))
    
    table_rows <- input$volcano_genes_table_rows_selected
    
    if (length(table_rows) > 0) {
      req(volcano_table_data())
      
      df <- volcano_table_data()
      selected_current <- volcano_selected_indices()
      
      # Re-calculate the current table view to get the correct original_indices for the selected rows
      current_table_result <- prepare_volcano_table(df, selected_current, input$volcano_table_mode, max_rows = 1000)
      
      # Map the selected table rows back to the full dataset indices
      selected_genes <- current_table_result$display$Gene[table_rows]
      full_indices <- which(df$gene %in% selected_genes)
      
      full_indices <- full_indices[!is.na(full_indices)]
      volcano_selected_indices(sort(full_indices))
      
      # Send to plot
      plot_indices_0based <- full_indices - 1
      session$sendCustomMessage("select_plot_points", list(
        indices = plot_indices_0based
      ))
    } else {
      volcano_selected_indices(integer(0))
      session$sendCustomMessage("clear_plot_selection", list())
    }
  })

  # Download handler
  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("volcano_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      showNotification("Use the widget's download button to export the plot", type = "info")
    }
  )


  
  # Global cache for already-loaded cloud results
  cloud_cache <- new.env(parent = emptyenv())

  # Reactive flag for whether results exist
  output$has_results <- reactive({
    !is.null(input$previous_result) && input$previous_result != ""
  })
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)

  # Download currently shown result
  output$download_result <- downloadHandler(
    filename = function() paste0(input$result_tool, "_results_", Sys.Date(), ".csv"),
    content = function(file) {
      data <- switch(input$result_tool,
        "DGE" = {
          selected_label <- input$previous_result
          if (!is.null(selected_label) && nzchar(selected_label)) {
            if (selected_label %in% names(results_store$all)) {
              results_store$all[[selected_label]]$data
            } else if (grepl("^DGE_", selected_label)) {
              version_prefix <- sub("^DGE_([^_]+)_.*", "\\1", selected_label)
              parquet_key <- sub("^DGE_[a-zA-Z0-9_-]+_", "", selected_label)
              zarr_url <- normalize_url(input$dataset)
              fs <- fsspec_filesystem()
              parquet_path <- sprintf("%s/versions/%s/results/%s.parquet", zarr_url, version_prefix, parquet_key)
              if (fs$exists(parquet_path)) {
                arrow::read_parquet(parquet_path)
              } else {
                NULL
              }
            } else {
              NULL
            }
          } else {
            results()
          }
        },
        "Clustering" = NULL,
        "Enrichment" = NULL
      )
      req(data)
      write.csv(data, file, row.names = FALSE)
    }
  )

  # Send data to JS
  observe({
    req(numeric_obs_keys())
    print("Sending numeric vars to JS")
    session$sendCustomMessage("updateNumericVars", numeric_obs_keys())
  })

  # Handle variable data requests from JS
  observeEvent(input$requestVarData, {
    req(input$requestVarData$var)
    var_name <- input$requestVarData$var
    
    # Fetch data using your function
    data <- values$lazy_data$get_obs_column(values$lazy_data$z, var_name)
    
    # Send to filter.js (for histogram)
    session$sendCustomMessage("updateVarData", list(
      var = var_name, 
      data = data
    ))
    
    # Send to scatterplot.js (for filtering)
    session$sendCustomMessage("updateNumericData", list(
      columnName = var_name,
      values = as.numeric(data)  # Ensure it's numeric
    ))
  })
}

# Clean up plan on app stop
# onStop(function() {
#   stopLocalWorkers(queue = "scanpy_jobs")
#   plan(sequential)
# })

# shinyApp(ui, server)
# shinyApp(ui = ui, server = server)