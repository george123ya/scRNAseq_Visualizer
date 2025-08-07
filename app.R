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

if (!file.exists("processed_data/relaxed_epdsc_annotated_data.h5")) {
  dir.create("processed_data", showWarnings = FALSE)
  message("📦 .h5 file not found, downloading from Google Drive...")
  
  # Use gdown with the FILE ID
  system("gdown https://drive.google.com/uc?id=1DwykyMuvohpo1aFQsQerlGsmLl-IiHpy -O processed_data/relaxed_epdsc_annotated_data.h5")
}

# Load required libraries
suppressPackageStartupMessages({
  library(shiny)
  library(shinyjs)
  library(shinydashboard)
  library(DT)
  library(plotly)
  library(base64enc)
  library(viridisLite)
  library(colorspace)
  library(jsonlite)
  library(reticulate)
  library(dplyr)
})

# options(shiny.host = "0.0.0.0")
# options(shiny.port = 3838)

# Source global configurations and helper functions
source("global.R")

# Configure Python environment for scanpy/anndata

# use_condaenv("sc_rna_env_python2", required = TRUE)
use_condaenv("shiny_app_env", conda = "/opt/conda/bin/conda", required = TRUE)


# ==============================================================================
# PYTHON INTERFACE SETUP
# ==============================================================================

# Import Python libraries for h5ad file handling
anndata <- import("anndata", delay_load = TRUE)
scanpy <- import("scanpy", delay_load = TRUE)

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

#' Load and process h5ad files
#' 
#' @param file_path Path to the .h5ad file
#' @return List containing processed data frames and metadata
#' @description Loads AnnData objects and extracts UMAP coordinates, 
#' cell annotations, gene expression matrices, and optional MAGIC data
load_h5ad <- function(file_path) {
  cat("🔄 Loading .h5ad file:", file_path, "\n")
  
  # Load AnnData object
  ad_obj <- anndata$read_h5ad(file_path)
  
  # Extract UMAP coordinates with fallback options
  umap_keys <- c("X_umap_normal", "X_umap", "X_umap_magic")
  umap <- NULL
  
  for (k in umap_keys) {
    if (k %in% ad_obj$obsm_keys()) {
      umap <- ad_obj$obsm[k]
      cat("📍 Found UMAP coordinates:", k, "\n")
      break
    }
  }
  
  if (is.null(umap)) {
    stop("❌ No UMAP coordinates found. Expected keys: ", paste(umap_keys, collapse = ", "))
  }
  
  # Create base dataframe with coordinates
  df <- as.data.frame(umap)
  colnames(df) <- c("x", "y")
  
  # Helper function to extract observation columns with fallback names
  get_obs_col <- function(possible_names) {
    for (name in possible_names) {
      if (name %in% names(ad_obj$obs)) {
        cat("📊 Found annotation:", name, "\n")
        return(as.factor(ad_obj$obs[[name]]))
      }
    }
    return(NULL)
  }
  
  # Extract common cell annotations
  df$cluster <- get_obs_col(c("pheno_100", "leiden", "cluster_feature", "Cluster_Fine"))
  df$cell_type <- get_obs_col(c("celltype", "Annotation"))
  df$condition <- get_obs_col(c("Condition", "Sample"))
  df$time <- get_obs_col(c("Time"))
  
  # Extract MAGIC UMAP coordinates if available
  if ("X_umap_magic" %in% ad_obj$obsm_keys()) {
    umap_magic <- ad_obj$obsm["X_umap_magic"]
    df$magic_x <- umap_magic[, 1]
    df$magic_y <- umap_magic[, 2]
    cat("✨ MAGIC UMAP coordinates available\n")
  } else {
    df$magic_x <- NA
    df$magic_y <- NA
    cat("⚠️  MAGIC UMAP coordinates not found\n")
  }
  
  # Get available expression layers
  available_layers <- reticulate::py_to_r(
    reticulate::import_builtins()$list(ad_obj$layers$keys())
  )
  cat("📋 Available layers:", paste(available_layers, collapse = ", "), "\n")
  
  # Extract expression data with priority order
  layers_priority <- c("log1p_data", "raw", "lognorm_pseudocount.1")
  expr <- NULL
  
  for (layer in layers_priority) {
    if (layer %in% available_layers) {
      expr <- ad_obj$layers[layer]
      cat("📈 Using expression layer:", layer, "\n")
      break
    }
  }
  
  # Extract MAGIC imputed expression if available
  expr_magic <- NULL
  if ("MAGIC_imputed_data" %in% available_layers) {
    expr_magic <- ad_obj$layers["MAGIC_imputed_data"]
    cat("✨ MAGIC imputed expression available\n")
  }
  
  # Extract gene names
  genes <- as.character(py_to_r(ad_obj$var_names$tolist()))
  cat("🧬 Total genes:", length(genes), "\n")
  
  # Extract SEACells summary if available
  seacell_df <- NULL
  if ("SEACells_summary" %in% names(ad_obj$uns)) {
    seacell_df <- as.data.frame(ad_obj$uns[["SEACells_summary"]]$obsm["X_umap"])
    colnames(seacell_df) <- c("UMAP_1", "UMAP_2")
    seacell_df$SEACell <- rownames(seacell_df)
    seacell_df$cell_type <- ad_obj$uns[["SEACells_summary"]]$obs[["cluster_feature"]]
    cat("🔬 SEACells metacells available:", nrow(seacell_df), "\n")
  }
  
  cat("✅ Successfully loaded:", nrow(df), "cells\n\n")
  
  return(list(
    df = df,
    expr = expr,
    expr_magic = expr_magic,
    genes = genes,
    ad = ad_obj,
    seacell_df = seacell_df
  ))
}

# ==============================================================================
# GENE EXPRESSION EXTRACTION
# ==============================================================================

#' Extract and encode gene expression data for visualization
#' 
#' @param gene_names Character vector of gene names to extract
#' @param h5_data List containing loaded h5ad data
#' @return List with encoded expression data and statistics
#' @description Efficiently extracts gene expression data from sparse matrices
#' and encodes it for JavaScript visualization
extract_gene_data <- function(gene_names, h5_data) {
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
    # Match requested genes to available genes
    gene_indices <- match(gene_names, h5_data$genes)
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
    
    valid_genes <- gene_names[valid_mask]
    valid_indices <- gene_indices[valid_mask]
    
    cat("🔍 Extracting", length(valid_genes), "genes:", paste(valid_genes, collapse = ", "), "\n")
    
    # Efficient sparse matrix extraction
    if (inherits(h5_data$expr, "dgRMatrix")) {
      cat("⚡ Fast dgRMatrix extraction\n")
      gene_rows <- h5_data$expr[valid_indices, , drop = FALSE]
      expr_matrix <- as.matrix(Matrix::t(gene_rows))
    } else if (inherits(h5_data$expr, "dgCMatrix")) {
      cat("⚡ Fast dgCMatrix extraction\n")
      expr_matrix <- as.matrix(h5_data$expr[, valid_indices, drop = FALSE])
    } else {
      expr_matrix <- as.matrix(h5_data$expr[, valid_indices, drop = FALSE])
    }
    
    cat("📊 Expression matrix:", nrow(expr_matrix), "cells ×", ncol(expr_matrix), "genes\n")
    
    # Calculate gene expression statistics
    gene_ranges <- list()
    if (ncol(expr_matrix) > 1) {
      col_mins <- apply(expr_matrix, 2, min, na.rm = TRUE)
      col_maxs <- apply(expr_matrix, 2, max, na.rm = TRUE)
      col_means <- colMeans(expr_matrix, na.rm = TRUE)
      
      for (i in seq_along(valid_genes)) {
        gene_ranges[[i]] <- list(
          min = col_mins[i],
          max = col_maxs[i],
          mean = col_means[i]
        )
      }
    } else {
      # Single gene case
      vals <- expr_matrix[, 1]
      gene_ranges[[1]] <- list(
        min = min(vals, na.rm = TRUE),
        max = max(vals, na.rm = TRUE),
        mean = mean(vals, na.rm = TRUE)
      )
    }
    names(gene_ranges) <- valid_genes
    
    # Encode expression data for JavaScript
    gene_data_vector <- as.vector(expr_matrix)
    binary_data <- writeBin(as.numeric(gene_data_vector), raw())
    encoded_data <- base64enc::base64encode(binary_data)
    
    cat("💾 Encoded data size:", nchar(encoded_data), "characters\n")
    
    return(list(
      genes = valid_genes,
      data = encoded_data,
      magic_data = "",  # Add MAGIC data extraction if needed
      ranges = gene_ranges,
      magic_ranges = list(),
      nrows = nrow(expr_matrix),
      ncols = ncol(expr_matrix)
    ))
    
  }, error = function(e) {
    cat("❌ Error extracting gene data:", e$message, "\n")
    return(list(
      genes = character(0),
      data = character(0),
      magic_data = character(0),
      ranges = list(),
      magic_ranges = list()
    ))
  })
}

# ==============================================================================
# DEMO DATA GENERATION
# ==============================================================================

#' Generate realistic demo scRNA-seq data
#' 
#' @param n Number of cells to generate
#' @return Data frame with simulated single-cell data
#' @description Creates biologically plausible scRNA-seq data for testing
#' and demonstration purposes
generate_demo_data <- function(n) {
  cat("🎭 Generating demo data with", n, "cells\n")
  
  # Define realistic cluster centers in UMAP space
  centers <- matrix(c(
    -1.2, -0.8,  # T cells
    -0.5, -1.0,  # B cells  
     0.3, -0.7,  # NK cells
     1.0, -0.2,  # Monocytes
     0.8,  0.8,  # Dendritic cells
    -0.8,  0.5,  # Macrophages
     0.0,  0.9,  # Neutrophils
    -1.5,  0.2   # Stem cells
  ), ncol = 2, byrow = TRUE)
  
  k <- nrow(centers)
  
  # Realistic cluster proportions
  cluster_probs <- c(0.20, 0.15, 0.12, 0.18, 0.10, 0.08, 0.12, 0.05)
  cluster_ids <- sample(1:k, n, replace = TRUE, prob = cluster_probs)
  
  # Generate coordinates with cluster-specific noise
  xs <- rnorm(n, mean = centers[cluster_ids, 1], sd = 0.12)
  ys <- rnorm(n, mean = centers[cluster_ids, 2], sd = 0.12)
  
  # Create biological annotations
  cell_type_map <- c("T_cells", "B_cells", "NK_cells", "Monocytes", 
                     "Dendritic", "Macrophages", "Neutrophils", "Stem_cells")
  cell_types <- cell_type_map[cluster_ids]
  
  # Add biological mixing (10% of cells change type)
  mix_indices <- sample(1:n, round(n * 0.1))
  cell_types[mix_indices] <- sample(cell_type_map, length(mix_indices), replace = TRUE)
  
  # Generate experimental metadata
  samples <- paste0("Sample_", sample(LETTERS[1:4], n, replace = TRUE, 
                                      prob = c(0.3, 0.25, 0.25, 0.2)))
  conditions <- sample(c("Control", "Treated", "Vehicle"), n, replace = TRUE, 
                      prob = c(0.4, 0.35, 0.25))
  timepoints <- sample(c("0h", "6h", "12h", "24h", "48h"), n, replace = TRUE, 
                      prob = c(0.25, 0.2, 0.2, 0.2, 0.15))
  
  # Generate realistic gene expression with dropout
  cluster_expr_means <- c(0.1, 0.2, 0.4, 0.7, 0.9, 0.3, 0.8, 0.05)
  cluster_expr_sds <- c(0.15, 0.18, 0.20, 0.15, 0.12, 0.20, 0.15, 0.25)
  
  base_expr <- rnorm(n, 
                     mean = cluster_expr_means[cluster_ids], 
                     sd = cluster_expr_sds[cluster_ids])
  
  # Add spatial gradients
  spatial_gradient <- 0.3 * ((xs + 2) / 4) + 0.2 * ((ys + 1.5) / 3)
  
  # Simulate dropout (common in scRNA-seq)
  dropout_prob <- pmax(0.05, pmin(0.6, 0.4 - base_expr))
  dropout_mask <- runif(n) > dropout_prob
  
  # Generate multiple gene expressions
  gene_expr1 <- pmax(0, pmin(10, (base_expr + 0.3 * spatial_gradient + 
                                  rnorm(n, 0, 0.1)) * 10))
  gene_expr2 <- pmax(0, pmin(10, (base_expr + 0.4 * spatial_gradient + 
                                  rnorm(n, 0, 0.1)) * 10))
  gene_expr3 <- pmax(0, pmin(10, (base_expr + 0.5 * spatial_gradient + 
                                  rnorm(n, 0, 0.1)) * 10))
  
  # Apply dropout
  gene_expr1[!dropout_mask] <- 0
  gene_expr2[!dropout_mask] <- 0
  gene_expr3[!dropout_mask] <- 0
  
  # Create MAGIC coordinates (slightly transformed)
  magic_x <- xs + rnorm(n, mean = 0.1, sd = 0.05) + 0.6 * sin(ys * 3)
  magic_y <- ys + rnorm(n, mean = 0.1, sd = 0.05) + 0.6 * cos(xs * 3)
  
  # Combine all data
  demo_data <- data.frame(
    x = xs,
    y = ys,
    magic_x = magic_x,
    magic_y = magic_y,
    cluster = paste0("Cluster_", cluster_ids),
    cell_type = cell_types,
    condition = conditions,
    time = timepoints,
    sample_id = samples,
    gene_expression1 = gene_expr1,
    gene_expression2 = gene_expr2,
    gene_expression3 = gene_expr3,
    stringsAsFactors = FALSE
  )
  
  cat("✅ Generated demo data:", nrow(demo_data), "cells with", 
      length(unique(demo_data$cluster)), "clusters\n")
  
  return(demo_data)
}

# ==============================================================================
# USER INTERFACE
# ==============================================================================

ui <- fluidPage(
  theme = bslib::bs_theme(version = 4, bootswatch = "flatly"),
  useShinyjs(),
  
  # Custom CSS styling
  tags$head(
    tags$script(src = "scatterplot.js"),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
  ),
  
  # Main Title Section
  div(class = "main-title",
    h1("🧬 scRNA-seq Visualizer"),
    div(class = "subtitle", 
        "Interactive single-cell RNA sequencing data exploration platform")
  ),
  
  # Main Layout
  fluidRow(
    # Sidebar Panel
    column(3,
      # Dataset Selection
      div(class = "well",
        div(class = "section-header", "📁 Dataset"),
        uiOutput("datasetUI"),
        tags$small(class = "text-muted", 
          "Select your .h5ad file or use demo data for exploration")
      ),
      
      # Gene Search and Visualization
      div(class = "well",
        div(class = "section-header", "🔬 Gene Expression"),
        uiOutput("geneSearchUI"),
        uiOutput("colorByUI"),
        div(class = "mt-2",
          checkboxInput("activateMAGIC", 
            HTML("<strong>🪄 MAGIC Imputation</strong>"), 
            value = FALSE),
          tags$small(class = "text-muted", 
            "Note: For visualization only, not for statistical analysis")
        )
      ),
      
      # Status and Data Info
      div(class = "info-card",
        div(id = "status-container",
          textOutput("status")
        ),
        uiOutput("dataInfo")
      ),
      
      # Advanced Tools
      div(class = "well",
        div(class = "section-header", "🛠️ Analysis Tools"),
        
        # Gene Set Signature
        div(class = "tool-section",
          checkboxInput("show_gene_set", 
            HTML("<strong>📊 Gene Set Signature</strong>"), FALSE),
          conditionalPanel(
            condition = "input.show_gene_set == true",
            selectizeInput("gene_set_input", "Select Genes:", 
              choices = NULL, multiple = TRUE,
              options = list(maxItems = 10, placeholder = "Choose genes...")),
            actionButton("calc_gene_score", "Calculate Signature", 
              class = "btn-primary btn-sm")
          )
        ),
        
        # SEACell Toggle
        div(class = "tool-section",
          checkboxInput("show_seacell_toggle", 
            HTML("<strong>🔬 SEACell Metacells</strong>"), FALSE),
          conditionalPanel(
            condition = "input.show_seacell_toggle == true",
            radioButtons("cell_level", "Display Level:", 
              choices = c("Single Cells" = "single", "SEACells" = "meta"),
              selected = "single", inline = TRUE)
          )
        ),
        
        # Differential Expression
        div(class = "tool-section",
          checkboxInput("show_dge", 
            HTML("<strong>📈 Differential Expression</strong>"), FALSE),
          conditionalPanel(
            condition = "input.show_dge == true",
            selectInput("dge_group_by", "Group by:", 
              choices = c("Cluster" = "cluster", "Cell Type" = "cell_type")),
            selectInput("dge_condition_1", "Condition A:", choices = NULL),
            selectInput("dge_condition_2", "Condition B:", choices = NULL),
            actionButton("run_dge", "Run DE Test", class = "btn-primary btn-sm"),
            tags$small(class = "text-muted", 
              "Click a gene in results to view violin plot.")
          )
        ),

        # Differential Expression with Milo
        div(class = "tool-section",
          checkboxInput("show_milo", 
            HTML("<strong>📊 Differential Expression with Milo</strong>"), FALSE),
          conditionalPanel(
            condition = "input.show_milo == true",
            selectInput("milo_group_by", "Group by:", 
              choices = c("Cluster", "Cell Type")),
            selectInput("milo_condition_1", "Condition A:", choices = NULL),
            selectInput("milo_condition_2", "Condition B:", choices = NULL),
            actionButton("run_milo", "Run Milo Analysis", class = "btn-primary btn-sm"),
            tags$small(class = "text-muted", 
              "Milo requires a pre-computed graph. Ensure you have run the necessary steps.")
          )
        ),
      )
    ),

    mainPanel(
      div(style = "position: relative; width: 1100px; height: 800px;",
          tags$div(id = "loadingOverlay", 
                   tags$div(class = "spinner"), 
                   tags$span("Processing data...")),
          tags$canvas(id = "scatterplot_canvas", width = 800, height = 400, 
                     style = "border:1px solid #ccc; border-radius: 4px;")
      ),
      # Legends container for your dynamic legends
      tags$div(id = "legendsContainer", style = "margin-top: 15px;"),
      verbatimTextOutput("selected_points")
    )
  ),

  # Inject viridis colors for gene expression palette to JS
  tags$script(HTML(sprintf("
    const viridisColors = %s;
    let scatterplot;
    let clusters = [];
    let geneExprRange = [];
    let clusterColors = [];
    let points = [];
    let filteredPoints = [];
    let annotations = {};
    let currentColorBy = '';
    let numAnnotations = 0;
    let useTransition = false;

  ", toJSON(viridis_hex, auto_unbox = TRUE))))
)

# Fixed server function with proper initialization

server <- function(input, output, session) {
  # Reactive values to store processed data
  values <- reactiveValues(
    processed_data = NULL,
    annotation_names = NULL,
    h5_data = NULL,  # Store h5_data in reactive values instead of global
    current_dataset = NULL
  )
  
  observeEvent(input$activateMAGIC, {
    session$sendCustomMessage("toggleMAGIC", input$activateMAGIC)
  })
  
  # FIXED: Move the data generation/processing logic into a separate observer
  # that only runs when h5_data is available
  observeEvent(values$h5_data, {
    req(values$h5_data)  # Ensure h5_data exists
    
    session$sendCustomMessage("showSpinner", TRUE)
    output$status <- renderText("Processing data...")

    # Use the loaded h5_data
    raw_data <- values$h5_data$df

    print(head(raw_data))  # Debugging: print first few rows
    
    # Process the dataframe
    processed_result <- process_dataframe(raw_data)
    values$processed_data <- processed_result
    
    # Create color by choices
    annotation_names <- names(processed_result$annotations)
    values$annotation_names <- annotation_names

    print(head(processed_result$data))
    
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
    
    output$status <- renderText(paste0("Rendered ", nrow(raw_data), " cells with ", 
                                      length(annotation_names), " annotations"))
  })

  # FIXED: Enhanced dataset loading with proper error handling
  observeEvent(input$dataset, {
    req(input$dataset)
    
    # Skip if placeholder values
    if (input$dataset %in% c("", "❌ Folder not found", "❌ No .h5ad files found")) {
      return()
    }
    
    # Check if file actually exists
    if (!file.exists(input$dataset)) {
      showNotification(
        paste("File not found:", input$dataset), 
        type = "error", 
        duration = 5
      )
      return()
    }
    
    # Show loading notification
    showNotification(
      paste("Loading dataset:", basename(input$dataset)), 
      type = "message", 
      duration = 3
    )
    
    tryCatch({
      # Load the dataset and store in reactive values
      values$current_dataset <- input$dataset
      values$h5_data <- load_h5ad(input$dataset)  # Store in reactive values
      
      # Success notification
      showNotification(
        paste("✅ Successfully loaded:", basename(input$dataset)), 
        type = "message", 
        duration = 3
      )
      
      cat("✅ Loaded dataset:", basename(input$dataset), "\n")
      cat("   Path:", input$dataset, "\n")
      cat("   Cells:", nrow(values$h5_data$df), "\n")
      cat("   Genes:", length(values$h5_data$genes), "\n")
      
    }, error = function(e) {
      # Error notification
      showNotification(
        paste("❌ Error loading file:", e$message), 
        type = "error", 
        duration = 10
      )
      cat("❌ Error loading dataset:", e$message, "\n")
      
      # Clear any partial data
      values$h5_data <- NULL
      values$processed_data <- NULL
    })
  })

  # UI components
  output$datasetUI <- renderUI({
    folder_path <- "./processed_data"
    
    # Check if folder exists
    if (!dir.exists(folder_path)) {
      return(
        div(
          selectInput("dataset", "Select Dataset:", 
                     choices = list("❌ Folder not found" = ""), 
                     selected = ""),
          helpText(paste("Path not found:", folder_path), style = "color: red;")
        )
      )
    }
    
    # Get all .h5ad files recursively
    files <- list.files(folder_path, full.names = TRUE, recursive = TRUE)
    
    if (length(files) == 0) {
      return(
        div(
          selectInput("dataset", "Select Dataset:", 
                     choices = list("❌ No .h5ad files found" = ""), 
                     selected = ""),
          helpText(paste("No .h5ad files found in:", folder_path), style = "color: orange;")
        )
      )
    }
    
    # Create user-friendly display names
    display_names <- sapply(files, function(file) {
      rel_path <- gsub(paste0("^", folder_path, "/?"), "", file)
      clean_name <- gsub("\\.h5ad$", "", basename(rel_path))
      
      if (dirname(rel_path) != ".") {
        folder_name <- basename(dirname(rel_path))
        return(paste0(folder_name, " / ", clean_name))
      } else {
        return(clean_name)
      }
    })
    
    choices <- setNames(files, display_names)
    choices <- choices[order(names(choices))]
    
    file_info_text <- ""
    if (length(files) > 0) {
      first_file <- files[1]
      file_size <- file.size(first_file)
      file_size_mb <- round(file_size / (1024^2), 2)
      file_info_text <- paste0("Files found: ", length(files), " • First file size: ", file_size_mb, " MB")
    }
    
    return(
      div(
        selectInput("dataset", "Select Dataset:", 
                   choices = choices, 
                   selected = 'relaxed_epdsc_annotated_data.h5',
                   width = "100%"),
        helpText(file_info_text, style = "color: #666; font-size: 11px;")
      )
    )
  })

  # UI: Gene search (wait for data to be loaded)
  output$geneSearchUI <- renderUI({
    if (is.null(values$h5_data)) {
      selectizeInput("geneSearch", "🔍 Search Gene:",
                    choices = list("Load dataset first..." = "loading"),
                    multiple = TRUE)
    } else {
      selectizeInput(
        "geneSearch",
        "🔍 Search Gene:",
        choices = NULL,
        multiple = TRUE,
        options = list(
          placeholder = "Type or click to select genes...",
          openOnFocus = FALSE,
          closeAfterSelect = TRUE,
          plugins = list("remove_button")
        )
      )
    }
  })

  # FIXED: Gene search observer - check if h5_data exists
  observeEvent(input$geneSearch, {
    req(values$h5_data)  # Ensure h5_data is loaded
    
    if (!is.null(input$geneSearch) && length(input$geneSearch) > 0) {
      gene_list <- as.character(input$geneSearch)
      
      print(paste("🔍 Searching for genes:", paste(gene_list, collapse = ", ")))
      
      # Extract both normal and MAGIC gene expression data
      gene_data <- extract_gene_data(gene_list, values$h5_data)
      
      print(paste("📊 Extracted data summary:"))
      print(paste("  - Genes found:", length(gene_data$genes)))
      print(paste("  - Normal data size:", nchar(gene_data$data)))
      print(paste("  - MAGIC data size:", nchar(gene_data$magic_data)))
      
      session$sendCustomMessage("geneSearchChange", list(
        genes = gene_list,
        expression_data = gene_data
      ))
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
    }
  }, ignoreNULL = FALSE)

  # FIXED: Update gene choices when data is loaded
  observeEvent(values$h5_data, {
    req(values$h5_data)
    
    updateSelectizeInput(
      session,
      "geneSearch",
      choices = values$h5_data$genes,
      server = TRUE
    )
  })
  
  # Dynamic UI for color selection
  output$colorByUI <- renderUI({
    if (is.null(values$annotation_names)) {
      selectInput("colorBy", "Color cells by:", choices = list("Load dataset first..." = "loading"))
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

  # FIXED: Test observer - only run when h5_data is available
  observe({
    req(values$h5_data)  # Only run when h5_data exists
    
    # Test with a few common genes if they exist
    test_genes <- intersect(c("CD3E", "CD4", "CD8A", "IL2", "IFNG"), values$h5_data$genes)
    if (length(test_genes) > 0) {
      print(paste("🧪 Testing with genes:", paste(test_genes[1:min(2, length(test_genes))], collapse = ", ")))
    }
  })
  
  # Data info output
  output$dataInfo <- renderUI({
    req(values$h5_data, values$processed_data)
    
    # Check if MAGIC data is available
    has_magic_coords <- !all(is.na(values$h5_data$df$magic_x))
    has_magic_expr <- !is.null(values$h5_data$expr_magic)
    
    result <- values$processed_data
    info_text <- tagList(
      h5("📋 Dataset Info:"),
      p(paste("🔹 Total cells:", nrow(values$h5_data$df))),
      p(paste("🔹 Total genes:", length(values$h5_data$genes))),
      p(paste("🔹 Annotations:", length(result$annotations))),
      p(paste("🔹 MAGIC coordinates:", ifelse(has_magic_coords, "✓ Available", "❌ Not available"))),
      p(paste("🔹 MAGIC expression:", ifelse(has_magic_expr, "✓ Available", "❌ Not available"))),
      if (length(result$annotations) > 0) {
        div(
          strong("Available annotations:"),
          br(),
          paste(names(result$annotations), collapse = ", ")
        )
      }
    )
    
    div(style = "background: #f8f9fa; padding: 10px; border-radius: 5px; font-size: 12px;", info_text)
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
}

# shinyApp(ui, server)
shinyApp(ui = ui, server = server)