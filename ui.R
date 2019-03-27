
library(plotly)

fluidPage(
  
  fluidRow(
    
    sidebarLayout(
      
      sidebarPanel(
        fluidRow(
          #sliderInput("N_samps", "Number of samples", 1, 8, 1),
          numericInput("N_samps", "Number of samples", value=1),
          selectInput("samp_select", "Sample to mark", choices=c(1, "Water"), selected=1),
          actionButton("load_layout", "Load Layout"),
          actionButton("save_layout", "Save Layout"),
          br("\n"),
          fileInput("pyro_file", "Pyro data file"),
          actionButton("start_analysis", "Start analysis")
        )
      ),
      
      mainPanel(
        div(
          style = "position:relative",
          plotOutput(
            "plate",
            click = "plate_click",
            dblclick = "plate_dblclick",
            brush = brushOpts(id = "plate_brush"),
            hover = hoverOpts(id = "plate_hover", delay = 200, delayType = "debounce")
          ),
          uiOutput("hover_info")
        )
      )
      
    )
    
  ),
  
  fixedRow(
    # uiOutput("assay_tabs")
    column(3,
           actionButton("remove_failed", "Remove failed CpGs"),
           br("\n"),
           selectInput("assay_select", "Assay", NULL),
           #checkboxInput("remove_failed", "Remove failed CpGs", value=TRUE),
           checkboxGroupInput("CpG_select", "CpGs"),
           br("\n\n"),
           downloadButton("export_tsv", "Export data (tsv)"),
           downloadButton("export_xlsx", "Export data (xlsx)")
    ),
    column(3,
           tableOutput("assay_table"),
           style="overflow-y: scroll"),
    column(6,
           plotlyOutput("assay_plot")
    )
  )
  
)