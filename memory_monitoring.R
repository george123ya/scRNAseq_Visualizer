# ==============================================================================
# RAM Usage Monitoring Functions for Shiny scRNA-seq App
# ==============================================================================

# Function to get current memory usage
get_memory_usage <- function() {
  tryCatch({
    # Get R session memory usage
    mem_info <- gc(verbose = FALSE)
    
    # Calculate total memory used by R (in MB)
    r_memory_mb <- sum(mem_info[, 2]) * 8 / 1024  # Convert from 8-byte units to MB
    
    # Get system memory info (Linux/Unix specific)
    system_memory <- NULL
    if (file.exists("/proc/meminfo")) {
      meminfo <- readLines("/proc/meminfo", n = 10)
      total_line <- grep("^MemTotal:", meminfo, value = TRUE)
      available_line <- grep("^MemAvailable:", meminfo, value = TRUE)
      
      if (length(total_line) > 0 && length(available_line) > 0) {
        total_kb <- as.numeric(gsub(".*?([0-9]+).*", "\\1", total_line))
        available_kb <- as.numeric(gsub(".*?([0-9]+).*", "\\1", available_line))
        
        system_memory <- list(
          total_gb = round(total_kb / 1024 / 1024, 2),
          available_gb = round(available_kb / 1024 / 1024, 2),
          used_gb = round((total_kb - available_kb) / 1024 / 1024, 2),
          used_percent = round((total_kb - available_kb) / total_kb * 100, 1)
        )
      }
    }
    
    # Alternative for cross-platform (less detailed)
    if (is.null(system_memory)) {
      # Use R's memory.limit() on Windows or estimate
      if (.Platform$OS.type == "windows") {
        mem_limit_mb <- memory.limit()
        system_memory <- list(
          total_gb = round(mem_limit_mb / 1024, 2),
          available_gb = NA,
          used_gb = round(r_memory_mb / 1024, 2),
          used_percent = round(r_memory_mb / mem_limit_mb * 100, 1)
        )
      }
    }
    
    return(list(
      r_session_mb = round(r_memory_mb, 2),
      r_session_gb = round(r_memory_mb / 1024, 2),
      system_memory = system_memory,
      timestamp = Sys.time()
    ))
    
  }, error = function(e) {
    return(list(
      r_session_mb = NA,
      r_session_gb = NA,
      system_memory = NULL,
      error = e$message,
      timestamp = Sys.time()
    ))
  })
}

# Function to format memory usage for display
format_memory_display <- function(mem_info) {
  if (is.na(mem_info$r_session_mb)) {
    return("Memory info unavailable")
  }
  
  r_mem_text <- paste0("R Session: ", mem_info$r_session_gb, " GB")
  
  if (!is.null(mem_info$system_memory)) {
    sys_mem <- mem_info$system_memory
    system_text <- paste0(
      "System: ", sys_mem$used_gb, "/", sys_mem$total_gb, " GB", 
      " (", sys_mem$used_percent, "%)"
    )
    return(paste(r_mem_text, system_text, sep = " | "))
  } else {
    return(r_mem_text)
  }
}

# Function to get Python process memory (if using reticulate)
get_python_memory <- function() {
  tryCatch({
    # Check if reticulate is loaded and Python is available
    if (requireNamespace("reticulate", quietly = TRUE) && reticulate::py_available()) {
      # Import psutil if available
      psutil <- tryCatch(
        reticulate::import("psutil"),
        error = function(e) NULL
      )
      
      if (!is.null(psutil)) {
        # Get current process
        process <- psutil$Process()
        mem_info <- process$memory_info()
        
        return(list(
          python_memory_mb = round(mem_info$rss / 1024 / 1024, 2),
          python_memory_gb = round(mem_info$rss / 1024 / 1024 / 1024, 2)
        ))
      }
    }
    
    return(NULL)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Function to estimate AnnData object memory usage
estimate_anndata_memory <- function(lazy_data) {
  if (is.null(lazy_data)) return(0)
  
  tryCatch({
    # Estimate based on data types and dimensions
    n_cells <- nrow(lazy_data$df)
    n_genes <- length(lazy_data$genes)
    
    # Rough estimates (in MB)
    coords_mb <- (n_cells * 2 * 8) / 1024 / 1024  # x,y coordinates
    annotations_mb <- (n_cells * length(lazy_data$df) * 4) / 1024 / 1024  # annotations
    metadata_mb <- 1  # gene names, etc.
    
    if (lazy_data$memory_strategy == "true_lazy") {
      # Minimal footprint for true lazy loading
      estimated_mb <- coords_mb + annotations_mb + metadata_mb
    } else {
      # Add some overhead for selective loading
      estimated_mb <- (coords_mb + annotations_mb + metadata_mb) * 1.5
    }
    
    return(round(estimated_mb, 2))
    
  }, error = function(e) {
    return(0)
  })
}

# Memory monitoring reactive function for Shiny
create_memory_monitor <- function(session, update_interval = 5000) {
  # Create a reactive timer
  timer <- reactiveTimer(update_interval)  # Update every 5 seconds
  
  # Reactive expression for memory info
  memory_info <- reactive({
    timer()  # Depend on timer
    
    # Get current memory usage
    mem_info <- get_memory_usage()
    python_mem <- get_python_memory()
    
    # Combine information
    if (!is.null(python_mem)) {
      mem_info$python_memory_mb <- python_mem$python_memory_mb
      mem_info$python_memory_gb <- python_mem$python_memory_gb
    }
    
    return(mem_info)
  })
  
  return(memory_info)
}

# UI component for memory display
memory_display_ui <- function() {
  div(
    style = "background: #f8f9fa; padding: 8px 12px; margin: 5px 0; border-radius: 4px; font-family: monospace; font-size: 11px;",
    div(
      style = "display: flex; justify-content: space-between; align-items: center;",
      span("💾 Memory Usage:", style = "font-weight: bold; color: #495057;"),
      span(id = "memory-usage-text", "Loading...", style = "color: #007bff;")
    ),
    div(
      id = "memory-details",
      style = "margin-top: 4px; font-size: 10px; color: #6c757d;"
    )
  )
}

# JavaScript code to update memory display
memory_display_js <- "
Shiny.addCustomMessageHandler('updateMemoryDisplay', function(message) {
  const memoryText = document.getElementById('memory-usage-text');
  const memoryDetails = document.getElementById('memory-details');
  
  if (memoryText) {
    memoryText.textContent = message.main_text;
  }
  
  if (memoryDetails && message.details) {
    memoryDetails.innerHTML = message.details;
  }
});
"