# ==============================================================================
# DEBUGGING SCRIPT - Run this to find the server function error
# ==============================================================================

cat("🔍 SHINY SERVER DEBUGGING\n")
cat("=" , rep("=", 79), "\n\n")

# Step 1: Check if server.R can be sourced
cat("Step 1: Checking syntax of server.R...\n")
tryCatch({
  source("server.R", echo = FALSE)
  cat("✅ server.R sourced successfully\n\n")
}, error = function(e) {
  cat("❌ ERROR in server.R:\n")
  cat("   ", e$message, "\n\n")
  cat("Line with error (look for syntax issues):\n")
  print(e)
})

# Step 2: Check if server is defined
cat("Step 2: Checking if 'server' is defined...\n")
if (exists("server")) {
  cat("✅ 'server' exists\n")
  cat("   Class:", class(server), "\n")
  cat("   Type: Function?\n\n")
} else {
  cat("❌ 'server' does not exist!\n")
  cat("   Make sure you define: server <- function(input, output, session) { ... }\n\n")
}

# Step 3: Check if ui is defined
cat("Step 3: Checking if 'ui' is defined...\n")
if (exists("ui")) {
  cat("✅ 'ui' exists\n\n")
} else {
  cat("❌ 'ui' does not exist!\n")
  cat("   Make sure you source('ui.R') or define ui\n\n")
}

# Step 4: Try to create minimal app
cat("Step 4: Testing minimal Shiny app...\n")
tryCatch({
  if (exists("ui") && exists("server")) {
    app <- shiny::shinyApp(ui = ui, server = server)
    cat("✅ Shiny app created successfully!\n")
    cat("   Try: shiny::runApp('.', port = 8630)\n\n")
  } else {
    cat("❌ Cannot create app - missing ui or server\n\n")
  }
}, error = function(e) {
  cat("❌ ERROR creating Shiny app:\n")
  cat("   ", e$message, "\n\n")
})

# Step 5: Check for common issues
cat("Step 5: Checking for common issues...\n")

# Check for unclosed braces
server_text <- tryCatch({
  readLines("server.R")
}, error = function(e) NULL)

if (!is.null(server_text)) {
  opening_braces <- sum(grepl("\\{", server_text))
  closing_braces <- sum(grepl("\\}", server_text))
  
  cat("   Opening braces {:", opening_braces, "\n")
  cat("   Closing braces }:", closing_braces, "\n")
  
  if (opening_braces != closing_braces) {
    cat("   ⚠️  WARNING: Brace mismatch! (", 
        opening_braces - closing_braces, "unclosed)\n")
  } else {
    cat("   ✅ Braces match\n")
  }
}

cat("\n" , rep("=", 80), "\n")
cat("NEXT STEPS:\n")
cat("1. Fix any errors listed above\n")
cat("2. Try: shiny::runApp('.', port = 8630)\n")
cat("3. If still broken, share the error message\n")
cat("=" , rep("=", 80), "\n")