#!/usr/bin/env Rscript
# ==============================================================================
# SHINY APP LAUNCHER - Run this instead of shiny::runApp()
# ==============================================================================

cat("🚀 Starting Shiny Application...\n\n")

# Load required libraries
library(shiny)

# Source UI and Server
cat("📦 Loading UI...\n")
source("ui.R", local = TRUE)

cat("📦 Loading Server...\n")
source("server.R", local = TRUE)

# Verify both exist
if (!exists("ui")) {
  stop("❌ UI not defined! Check ui.R")
}

if (!exists("server")) {
  stop("❌ Server not defined! Check server.R")
}

cat("✅ UI loaded\n")
cat("✅ Server loaded\n\n")

# Create and run app
cat("🎨 Creating Shiny app...\n")
app <- shinyApp(ui = ui, server = server)

cat("✅ App created successfully!\n")
cat("🌐 Launching on http://127.0.0.1:8631\n\n")

# Run the app
shiny::runApp(app, port = 8631, launch.browser = TRUE)