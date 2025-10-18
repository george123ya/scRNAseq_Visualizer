# ==============================================================================
# USER INTERFACE
# ==============================================================================

library(shinyjs)
library(shiny)
# library(shinydashboard)
# library(DT)
library(plotly)
library(jsonlite)
library(fontawesome)
library(reglScatterplot)

# ==============================================================================
# Enhanced scRNA-seq Interactive Visualizer UI with QC Tab
# ==============================================================================

ui <- fluidPage(
  theme = bslib::bs_theme(version = 4, bootswatch = "flatly"),
  useShinyjs(),
  
  # Custom CSS styling
  tags$head(
    tags$script(src = "scatterplot.js"),
    tags$script(src = "qc_plots.js"),
    tags$script(src = "https://d3js.org/d3.v7.min.js"),
    tags$script(src = 'https://cdnjs.cloudflare.com/ajax/libs/html2canvas/1.4.1/html2canvas.min.js'),
    tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
    tags$style(HTML("
      .nav-tabs .nav-link.active {
        background-color: #007bff;
        color: white;
        border-color: #007bff;
      }
      .qc-metric-card {
        background: #f8f9fa;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        padding: 15px;
        margin-bottom: 15px;
      }
      .filter-row {
        margin-bottom: 15px;
      }
      
      .filter-label {
        margin-bottom: 8px;
        font-size: 14px;
        color: #495057;
      }
      
      .filter-preview {
        border-left: 3px solid #007bff;
      }
      .qc-plot-container {
        margin-bottom: 20px;
        padding: 15px;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        background-color: #ffffff;
      }
      .metric-value {
        font-size: 24px;
        font-weight: bold;
        color: #007bff;
      }
      .metric-label {
        font-size: 14px;
        color: #6c757d;
        margin-top: 5px;
      }
      
      /* Collapsible section styles */
      .collapsible-section {
        margin-bottom: 15px;
        border-radius: 8px;
        background: #f8f9fa;
      }
      
      .collapsible-header {
        padding: 0px 15px;
        background: linear-gradient(135deg, #3a8dde, #007bff); /* smooth gradient */
        color: white;
        cursor: pointer;
        border-radius: 10px 10px 0 0;
        user-select: none;
        font-weight: 600;
        font-size: 15px;
        display: inline-flex;
        align-items: center;
        gap: 8px;  /* space between icon and text */
        box-shadow: 0 2px 5px rgba(0, 123, 255, 0.3); /* subtle shadow */
        transition: background 0.3s ease, box-shadow 0.3s ease;
        }

        .collapsible-header:hover {
        background: linear-gradient(135deg, #2e7ad1, #0056b3);
        box-shadow: 0 4px 10px rgba(0, 86, 179, 0.5);
        }

        .collapse-icon {
        font-size: 16px;
        transition: transform 0.3s ease;
        display: inline-block;
        color: white;
        }

        .collapse-icon.collapsed {
        transform: rotate(-90deg);
        }
      
        .collapsible-content {
            padding: 15px;
            border-radius: 0px 15px 15px 15px;
            border: 1px solid #dee2e6;
        }
        
        .section-header {
            margin: 0;
            font-size: 14px;
            font-weight: 500;
        }

        .gene-plots {
          # max-height: 600px;  /* Adjust as needed */
          # overflow-y: auto;   /* Enable vertical scrolling */
          padding: 10px;
        }
    "))
  ),
  
  # Main Title Section
  div(class = "main-title",
    h1("🧬 scRNA-seq Browser"),
  ),
  
    # Main Layout with Tabs
    fluidRow(
    # Sidebar Panel (always visible)
    column(3,
      # titlePanel("RAM usage over time"),
      # plotlyOutput("ram_plot", height = "300px"),

      # Dataset Selection - Always Expanded
      div(class = "collapsible-section",
        div(class = "collapsible-header", onclick = "toggleCollapse('dataset-section')",
          span(class = "section-header", "📁 Dataset"),
          span(class = "collapse-icon", "▼")
        ),
        div(id = "dataset-section", class = "collapsible-content",
          uiOutput("datasetUI"),
          tags$small(class = "text-muted", 
            "Select your .h5ad file or use demo data for exploration")
        ),
        # actionButton("test_capture", "Prepare Download"),
        # downloadButton("download_report", "Download Full Report")
        uiOutput("download_button_ui")

      ),
      
      # Status and Data Info - Collapsible
      div(class = "collapsible-section",
        div(class = "collapsible-header", onclick = "toggleCollapse('status-section')",
          span(class = "section-header", "ℹ️ Data Information"),
          span(class = "collapse-icon", "▼")
        ),
        div(id = "status-section", class = "collapsible-content",
          div(id = "status-container",
            textOutput("status")
          ),
          uiOutput("dataInfo"),
          uiOutput("selectedInfo")
        )
      ),
      
      # Tab-specific controls (conditionally shown)
      conditionalPanel(
        condition = "input.main_tabs == 'visualization'",
        # Gene Search and Visualization - Collapsible
        div(class = "collapsible-section",
          div(class = "collapsible-header", onclick = "toggleCollapse('gene-section')",
            span(class = "section-header", "🔬 Gene Expression"),
            span(class = "collapse-icon", "▼")
          ),
          div(id = "gene-section", class = "collapsible-content",
            uiOutput("geneSearchUI"),
            # selectizeInput(
            #   "geneSearch",
            #   "🔍 Search Gene:",
            #   choices = NULL,
            #   multiple = TRUE,
            #   options = list(
            #     placeholder = "Type gene names...",
            #     create = TRUE,            # <--- this makes Enter immediate
            #     openOnFocus = FALSE,
            #     closeAfterSelect = TRUE,
            #     plugins = list("remove_button"),
            #     maxItems = 15
            #   )
            # ),
            uiOutput("colorByUI"),
            # div(class = "mt-2",
            #   checkboxInput("activateMAGIC", 
            #     HTML("<strong>🪄 MAGIC Imputation</strong>"), 
            #     value = FALSE),
            #   tags$small(class = "text-muted", 
            #     "Note: For visualization only, not for statistical analysis")
            # ),
            uiOutput("chooseMatrixUI")
          )
        )
      ),

      # Gene Expression Tab Controls
      conditionalPanel(
        condition = "input.main_tabs == 'violin_plot_genes'",
        div(class = "collapsible-section",
          div(class = "collapsible-header", onclick = "toggleCollapse('gene-params-section')",
            span(class = "section-header", "🔍 Gene Plot Parameters"),
            span(class = "collapse-icon", "▼")
          ),
          div(id = "gene-params-section", class = "collapsible-content",
            
            conditionalPanel(
              condition = "input.qc_plot_type != 'scatter'",
              uiOutput("gene_selected_ui"),
              uiOutput("gene_group_by_ui"),  # Add the Group by selector
              actionButton("update_plots", "Update Plots", class = "btn btn-primary btn-sm"),  # Add button
            ),
            
            # Plot type selection
            radioButtons("gene_plot_type", "Plot Type:",
              choices = list(
                "Violin Plot" = "violin",
                "Box Plot" = "box"
              ),
              selected = "violin",
              inline = TRUE
            ),

            checkboxInput("gene_show_points", "Show individual points", value = TRUE),
            checkboxInput("gene_log_scale", "Log scale y-axis", value = FALSE),
            sliderInput("gene_point_size", "Point size:", 
              min = 0.1, max = 2, value = 0.5, step = 0.1)
          )
        )
      ),

      # Heatmap

      conditionalPanel(
        condition = "input.main_tabs == 'heatmap_genes'",

        # --- PANEL 1: INDIVIDUAL GENES ---
        div(class = "collapsible-section",
          div(
            class = "collapsible-header",
            onclick = "toggleCollapse('heatmap-panel1')",
            span(class = "section-header", "Panel 1: Individual Genes"),
            span(class = "collapse-icon", "▼")
          ),
          div(
            id = "heatmap-panel1", 
            class = "collapsible-content",

            textInput(
              "heatmap_genes_panel1",
              label = "Genes to visualize (comma-separated):",
              placeholder = "e.g., KRT10, BRCA2"
            ),

            # Cell grouping settings (shared between panels)
            h5("Cell Grouping Settings"),
            radioButtons(
              "heatmap_cell_group_mode_panel1",
              label = "Group cells by:",
              choices = c(
                "Categorical variable (cluster cells, show means)" = "categorical",
                "None (show all cells, color by group)" = "none"
              ),
              selected = "categorical"
            ),
            uiOutput("heatmap_cluster_by_ui_panel1"),

            radioButtons(
              "heatmap_score_panel1",
              label = "Scale values by:",
              choices = c("Z-score" = "z_score", "No scale (X layer)" = "log"),
              selected = "z_score",
              inline = TRUE
            ),

            actionButton(
              "update_heatmap_panel1",
              "Update Panel 1",
              class = "btn btn-primary btn-sm",
              style = "margin-top: 8px; width: 100%;"
            )
          )
        ),

        # --- PANEL 2: GENE SETS ---
        div(class = "collapsible-section",
          div(
            class = "collapsible-header",
            onclick = "toggleCollapse('heatmap-panel2')",
            span(class = "section-header", "Panel 2: Gene Sets / Pathways"),
            span(class = "collapse-icon", "▼")
          ),
          div(
            id = "heatmap-panel2", 
            class = "collapsible-content",

            radioButtons(
              "heatmap_gene_group_mode_panel2",
              label = "Gene set source:",
              choices = c(
                "Custom gene sets (from .gmt)" = "custom_gmt",
                "Manual gene sets" = "manual"
              ),
              selected = "custom_gmt"
            ),

            # Pathway selector
            conditionalPanel(
              condition = "input.heatmap_gene_group_mode_panel2 == 'pathways'",
              uiOutput("heatmap_pathway_select_panel2")
            ),

            # File upload for .gmt file
            conditionalPanel(
              condition = "input.heatmap_gene_group_mode_panel2 == 'custom_gmt'",
              fileInput(
                "gmt_file_panel2",
                label = "Upload .gmt file:",
                accept = c(".gmt"),
                placeholder = "Select a .gmt file"
              ),
              conditionalPanel(
                condition = "output.gmt_loaded_panel2",
                uiOutput("heatmap_gmt_select_panel2")
              )
            ),

            # Manual gene set builder
            conditionalPanel(
              condition = "input.heatmap_gene_group_mode_panel2 == 'manual'",
              h5("Build Custom Gene Sets"),
              wellPanel(
                fluidRow(
                  column(4,
                    textInput("manual_geneset_name", "Gene Set Name:", placeholder = "e.g., My_Genes")
                  ),
                  column(6,
                    textInput("manual_geneset_genes", "Genes (comma-separated):", placeholder = "e.g., GENE1, GENE2, GENE3")
                  ),
                  column(2,
                    br(),
                    actionButton("add_manual_geneset", "Add", class = "btn btn-success btn-sm")
                  )
                )
              ),
              uiOutput("manual_genesets_display"),
              conditionalPanel(
                condition = "output.has_manual_genesets",
                selectizeInput(
                  "selected_manual_genesets",
                  "Select gene sets to plot:",
                  choices = NULL,
                  multiple = TRUE,
                  options = list(placeholder = "Select gene sets...")
                )
              )
            ),

            # Cell grouping settings (shared between panels)
            h5("Cell Grouping Settings"),
            radioButtons(
              "heatmap_cell_group_mode_panel2",
              label = "Group cells by:",
              choices = c(
                "Categorical variable (cluster cells, show means)" = "categorical",
                "None (show all cells, color by group)" = "none"
              ),
              selected = "categorical"
            ),
            uiOutput("heatmap_cluster_by_ui_panel2"),

            radioButtons(
              "heatmap_score_panel2",
              label = "Scale values by:",
              choices = c("Z-score" = "z_score", "No scale (X layer)" = "log"),
              selected = "z_score",
              inline = TRUE
            ),

            actionButton(
              "update_heatmap_panel2",
              "Update Panel 2",
              class = "btn btn-primary btn-sm",
              style = "margin-top: 8px; width: 100%;"
            ),

            # Space
            tags$hr(style = "margin: 10px 0; border-color: #e9ecef;"),

            conditionalPanel(
              condition = "input.heatmap_gene_group_mode_panel2 == 'custom_gmt' && output.gmt_loaded_panel2",

              h5("Pathways / Geneset UMAP"),
              
              uiOutput("umap_gmt_select_panel2"),
              
              selectInput(
                "umap_score_method_panel2",
                "Score aggregation:",
                choices = c("Mean" = "mean", "Sum" = "sum"),
                selected = "mean",
                width = "100%",
                # style = "margin-top: 8px;"
              ),
              
              actionButton(
                "update_umap_panel2",
                "Update UMAP Signature",
                class = "btn btn-primary btn-sm",
                # style = "margin-top: 8px; width: 100%;"
              )
            ),

          )
        )
      ),

      
      # QC Tab Controls
      conditionalPanel(
        condition = "input.main_tabs == 'qc'",
        div(class = "collapsible-section",
          div(class = "collapsible-header", onclick = "toggleCollapse('qc-params-section')",
            span(class = "section-header", "🔍 QC Parameters"),
            span(class = "collapse-icon", "▼")
          ),
          div(id = "qc-params-section", class = "collapsible-content",
            conditionalPanel(
              condition = "input.qc_plot_type != 'scatter'",
              uiOutput("qc_metric_ui"),
              uiOutput("qc_group_by_ui")
            ),

            conditionalPanel(
              condition = "input.qc_plot_type == 'scatter'",
              uiOutput("qc_x_metric_ui"),
              uiOutput("qc_y_metric_ui"),
              uiOutput("qc_color_by_ui")
            ),
            
            # Plot type selection
            radioButtons("qc_plot_type", "Plot Type:",
              choices = list(
                "Violin Plot" = "violin",
                # "Box Plot" = "box",
                # "Histogram" = "histogram",
                "Scatter Plot" = "scatter"
              ),
              selected = "violin",
              inline = TRUE
            ),

            checkboxInput("qc_show_points", "Show individual points", value = TRUE),
            checkboxInput("qc_log_scale", "Log scale (where applicable)", value = FALSE),
            sliderInput("qc_point_size", "Point size:", 
            min = 0.1, max = 2, value = 0.5, step = 0.1),
            radioButtons("qc_legend_position", "Label Position:",
               choices = c("Legend" = "Legend", "Bottom Labels" = "Bottom Labels"),
               selected = "Legend")  
          )
        )
      ),

      # Tools Tab Controls
      conditionalPanel(
        condition = "input.main_tabs == 'results'",
        # Gene Set Signature - Collapsible
        div(class = "collapsible-section",
          div(class = "collapsible-header", onclick = "toggleCollapse('gene-set-section')",
            span(class = "section-header", "📊 Gene Set Signature"),
            span(class = "collapse-icon collapsed", "▼")
          ),
          div(id = "gene-set-section", class = "collapsible-content", style = "display: none;",
            checkboxInput("show_gene_set", 
              HTML("<strong>Enable Gene Set Analysis</strong>"), FALSE),
            conditionalPanel(
              condition = "input.show_gene_set == true",
              selectizeInput("gene_set_input", "Select Genes:", 
                choices = NULL, multiple = TRUE,
                options = list(maxItems = 10, placeholder = "Choose genes...")),
              actionButton("calc_gene_score", "Calculate Signature", 
                class = "btn-primary btn-sm")
            )
          )
        ),
        
        # SEACell Toggle - Collapsible
        div(class = "collapsible-section",
          div(class = "collapsible-header", onclick = "toggleCollapse('seacell-section')",
            span(class = "section-header", "🔬 SEACell Metacells"),
            span(class = "collapse-icon collapsed", "▼")
          ),
          div(id = "seacell-section", class = "collapsible-content", style = "display: none;",
            checkboxInput("show_seacell_toggle", 
              HTML("<strong>Enable SEACell View</strong>"), FALSE),
            conditionalPanel(
              condition = "input.show_seacell_toggle == true",
              radioButtons("cell_level", "Display Level:", 
                choices = c("Single Cells" = "single", "SEACells" = "meta"),
                selected = "single", inline = TRUE)
            )
          )
        ),
        
        # Differential Expression - Collapsible
        div(class = "collapsible-section",
          div(class = "collapsible-header", onclick = "toggleCollapse('dge-section')",
            span(class = "section-header", "📈 Differential Expression"),
            span(class = "collapse-icon collapsed", "▼")
          ),
          div(id = "dge-section", class = "collapsible-content", style = "display: none;",
            checkboxInput("show_dge", 
              HTML("<strong>Enable DE Analysis</strong>"), FALSE),
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
          )
        ),

        # Milo Analysis - Collapsible
        div(class = "collapsible-section",
          div(class = "collapsible-header", onclick = "toggleCollapse('milo-section')",
            span(class = "section-header", "📊 Differential Expression with Milo"),
            span(class = "collapse-icon collapsed", "▼")
          ),
          div(id = "milo-section", class = "collapsible-content", style = "display: none;",
            checkboxInput("show_milo", 
              HTML("<strong>Enable Milo Analysis</strong>"), FALSE),
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
          )
        )
      )
    ),
    
    # Main Panel with Tabs
    column(9,
      tabsetPanel(id = "main_tabs",
        
        # Visualization Tab (original content)
        tabPanel("🎨 Visualization", value = "visualization",
          div(style = "position: relative; width: 100%; height: 800px; margin-top: 15px;",
            tags$div(id = "loadingOverlay", 
                     tags$div(class = "spinner"), 
                     tags$span("Processing data...")),
            tags$canvas(id = "scatterplot_canvas", width = 800, height = 400, 
                       style = "border:1px solid #ccc; border-radius: 4px;")
          ),
          # Legends container for dynamic legends
          tags$div(id = "legendsContainer", style = "margin-top: 15px;"),
          verbatimTextOutput("selected_points")
        ),

        # Gene Expression Tab
        tabPanel("🎻 Violin Plot", value = "violin_plot_genes",
          div(style = "margin-top: 20px;",
            # h3("📈 Gene Expression"),
            fluidRow(
              column(8,
                div(class = "gene-plots",
                  h4("Expression Distribution"),
                  uiOutput("gene_plots_ui")      # Dynamic plot outputs
                )
              )
            )
          )
        ),
        tabPanel("🗺️ Heatmap", value = "heatmap_genes",
          div(style = "margin-top: 20px;",
            # h3("📈 Gene Expression Heatmaps"),

            # First heatmap
            div(class = "gene-plots",
              h4("Panel 1: Genes Heatmap"),
              plotOutput("heatmapPlot1", height = "auto")
            ),

            # Add some spacing between panels
            br(), br(),

            # Second heatmap
            div(class = "gene-plots",
              h4("Panel 2: Pathways/Gene Sets Heatmap"),
              plotOutput("heatmapPlot2", height = "auto")
            ),

            br(), br(),

            # UMAP for gene sets
            my_scatterplotOutput("umapPlot", width = "100%", height = "600px")

          )
        ),
        tabPanel("🔍 Quality Control", value = "qc",
          div(style = "margin-top: 15px;",
            
            # Main QC Plots Row
            fluidRow(
              # QC Violin/Box Plots
              column(8,
                div(class = "qc-plot-container",
                  h4("QC Metrics Distribution"),
                  # plotOutput("qc_main_plot", height = "400px")
                  uiOutput("qc_main_plot")
                )
              ),
              
              # Enhanced filtering UI section
              column(4,
                div(class = "qc-plot-container",
                  h4("QC Statistics"),
                  tableOutput("qc_stats_table")
                ),
                
                # Interactive Filters Section
                div(class = "qc-plot-container",
                  h4("Interactive Filters"),
                  
                  # Mitochondrial Percentage Filter
                  div(class = "filter-row",
                    div(class = "filter-label", 
                      tags$strong("Mitochondrial %")
                    ),
                    fluidRow(
                      column(6,
                        selectInput("mito_var_select", 
                                  label = NULL,
                                  choices = c("Choose variable" = ""),
                                  selected = "",
                                  width = "100%")
                      ),
                      column(3,
                        numericInput("mito_min", 
                                    label = NULL,
                                    value = 0,
                                    min = 0,
                                    max = 100,
                                    step = 0.1,
                                    width = "100%")
                      ),
                      column(3,
                        numericInput("mito_max", 
                                    label = NULL,
                                    value = 5,
                                    min = 0,
                                    max = 100,
                                    step = 0.1,
                                    width = "100%")
                      )
                    )
                  ),
                  
                  tags$hr(style = "margin: 10px 0; border-color: #e9ecef;"),
                  
                  # Total Counts Filter
                  div(class = "filter-row",
                    div(class = "filter-label", 
                      tags$strong("Total Counts")
                    ),
                    fluidRow(
                      column(6,
                        selectInput("counts_var_select", 
                                  label = NULL,
                                  choices = c("Choose variable" = ""),
                                  selected = "",
                                  width = "100%")
                      ),
                      column(3,
                        numericInput("counts_min", 
                                    label = NULL,
                                    value = 1000,
                                    min = 0,
                                    step = 100,
                                    width = "100%")
                      ),
                      column(3,
                        numericInput("counts_max", 
                                    label = NULL,
                                    value = NA,
                                    min = 0,
                                    step = 100,
                                    width = "100%")
                      )
                    )
                  ),
                  
                  tags$hr(style = "margin: 10px 0; border-color: #e9ecef;"),
                  
                  # Number of Genes Filter
                  div(class = "filter-row",
                    div(class = "filter-label", 
                      tags$strong("Number of Genes")
                    ),
                    fluidRow(
                      column(6,
                        selectInput("genes_var_select", 
                                  label = NULL,
                                  choices = c("Choose variable" = ""),
                                  selected = "",
                                  width = "100%")
                      ),
                      column(3,
                        numericInput("genes_min", 
                                    label = NULL,
                                    value = 200,
                                    min = 0,
                                    step = 50,
                                    width = "100%")
                      ),
                      column(3,
                        numericInput("genes_max", 
                                    label = NULL,
                                    value = NA,
                                    min = 0,
                                    step = 50,
                                    width = "100%")
                      )
                    )
                  ),
                  
                  # Filter Preview
                  div(class = "filter-preview",
                    style = "margin-top: 15px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;",
                    h5("Filter Preview", style = "margin-bottom: 10px; color: #495057;"),
                    textOutput("filter_preview_enhanced"),
                    br(),
                    div(style = "font-size: 12px; color: #6c757d;",
                      "Cells passing current filters"
                    )
                  )
                )
              )
            ),
              
            # Detailed Tables Row
            fluidRow(
              column(12,
                div(class = "qc-plot-container",
                  h4("Detailed QC Table"),
                  p("Interactive table showing QC metrics for all cells. Use filters to explore data."),
                  DT::dataTableOutput("qc_detailed_table")
                )
              )
            )
          )
        ),
        
        # Results/Analysis Tab (for DE results, etc.)
        tabPanel("📊 Analysis Results", value = "results",
          div(style = "margin-top: 15px;",
            h3("🧬 Analysis Results"),
            
            # DE Results
            conditionalPanel(
              condition = "input.show_dge == true",
              div(class = "qc-plot-container",
                h4("📈 Differential Expression Results"),
                DT::dataTableOutput("de_results")
              )
            ),
            
            # Gene Set Results  
            conditionalPanel(
              condition = "input.show_gene_set == true",
              div(class = "qc-plot-container",
                h4("📊 Gene Set Signature Results"),
                verbatimTextOutput("gene_set_summary")
              )
            ),
            
            # Milo Results
            conditionalPanel(
              condition = "input.show_milo == true",
              div(class = "qc-plot-container",
                h4("📊 Milo Analysis Results"), 
                verbatimTextOutput("milo_results")
              )
            ),
            
            # Default message when no analysis is active
            conditionalPanel(
              condition = "input.show_dge == false && input.show_gene_set == false && input.show_milo == false",
              div(class = "qc-plot-container text-center",
                h4("🔬 No Analysis Active"),
                p("Enable analysis tools in the Visualization tab to see results here."),
                tags$ul(class = "text-left",
                  tags$li("Gene Set Signature Analysis"),
                  tags$li("Differential Expression Analysis"),
                  tags$li("Milo Differential Abundance Testing")
                )
              )
            )
          )
        )
      )
    )
  ),

  # JavaScript for collapsible functionality
  tags$script(HTML("
    function toggleCollapse(sectionId) {
      const content = document.getElementById(sectionId);
      const header = content.previousElementSibling;
      const icon = header.querySelector('.collapse-icon');
      
      if (content.style.display === 'none') {
        content.style.display = 'block';
        icon.classList.remove('collapsed');
      } else {
        content.style.display = 'none';
        icon.classList.add('collapsed');
      }
    }
    
    // Initialize collapsed state for some sections
    document.addEventListener('DOMContentLoaded', function() {
      // Keep dataset and status sections expanded by default
      // All others start collapsed except the first section in each tab
    });
  ")),

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
