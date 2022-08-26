#'Produces the UI for the Gene Counts by Taxa Module
#' @param id module id
#' @importFrom shiny actionButton HTML numericInput selectInput selectizeInput sliderInput radioButtons downloadButton tagList NS




RA_ui <- function(id){

  ns <- NS(id)
  #help button UI
  tagList(
    actionButton(ns("help"), label = "Help/About this Plot", icon = icon("question-circle", lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Create Barplots of Genes by Taxa</h4><br><br>"),
    actionButton(ns("start"), label = "Load Dataset", style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h5>Choose Data to Plot</h5><br><br>"),
    numericInput(inputId = ns("sig_level"),label = "Select Significance Threshold", value = 0.05),
    actionButton(ns("subset_sig"), "Subset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    selectInput(ns("KO_level"), label = "Select Orthology Metadata Category", choices = c("KO_Class", "KO_Subclass_1", "KO_Subclass_2")),
    selectInput(ns("KO_group"), label = "Select Gene Group", choices = NULL, multiple = FALSE),
    actionButton(ns("subset_ortho"), "Subset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    selectizeInput(ns("Gene"), "Select Gene to Plot", choices = NULL,multiple = FALSE),
    selectInput(ns("taxa_level"), "Select Taxa Level to Group by", choices = c("Phylum", "Class", "Order", "Family", "Genus", "Genus_Species"), multiple = FALSE),
    HTML("<br><br><h5>Adjust Plot Parameters</h5><br><br>"),
    sliderInput(ns("height"), "Plot Height", min = 100, max = 1500, value = 600),
    sliderInput(ns("width"), "Plot Width", min = 100, max = 1500, value = 700),
    sliderInput(ns("font_size"), "Font Size",min = 5, max = 30, value = 14),
    sliderInput(ns("font"), "Font Size", min = 5, max = 25, value = 9),
    sliderInput(ns("border_size"), "Border Size", min = 0.5, max = 5, value = 1.5),
    HTML("<br><br><h5>Download Plot</h5><br><br>"),
    radioButtons(ns("extension"), "Save As:",
                 choices = c("pdf", "png","svg", "pptx"), inline = TRUE),
    numericInput(ns("figure_height"), label = "Figure Height(cm)", value = 5),
    numericInput(ns("figure_width"), label = "Figure Width(cm)", value = 5),
    downloadButton(ns("download_plot"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747")

  )

}

#'Produces the Server for the Gene Counts by Taxa Module
#' @param id module id
#' @importFrom shiny moduleServer observeEvent reactiveValues updateSelectizeInput req updateSelectInput reactive renderUI downloadHandler showModal modalDialog renderPlot
#' @importFrom thematic thematic_shiny
#' @importFrom colourpicker colourInput
#' @importFrom ggplot2 ggsave
#' @importFrom officer read_pptx add_slide ph_with ph_location_type

RA_server <- function(id){

  moduleServer(id, function(input, output, session){
    thematic::thematic_shiny()

    gene_choice <- reactiveValues(df = NULL)

    observeEvent(input$start, {


      gene_choice$df <- readRDS("inst/extdata/HF_Strep_After_vs_Before.rds")


      updateSelectizeInput(inputId = "Gene", choices = unique(gene_choice$df[,"KO_NAME"]), selected = character(0))



    })

    observeEvent(input$subset_sig,{

      sig_level <- input$sig_level

      gene_data <- gene_choice$df

      gene_data<- gene_data[gene_data$padj <= sig_level,]

      gene_data <- gene_data[!is.na(gene_data$NAME),]

      updateSelectizeInput(inputId = "Gene", choices = unique(gene_data[,"KO_NAME"]),selected = character(0))

      gene_choice$df <- gene_data

    })

    #a reactive to update the fill_levels UI based on which fill_variable choice is selected as an input
    ortho_choice <- reactive({
      req(!is.null(gene_choice$df))
      ortho_data <- gene_choice$df
      ortho_variable <- ortho_data[,input$KO_level]
      return(ortho_variable)

    })

    #an observation which dynamically updates the fill_levels input options when the fill_variable input is changed
    observeEvent(ortho_choice(), {
      choices <- unique(ortho_choice())
      updateSelectInput(inputId = "KO_group", choices = choices)
    })

    observeEvent(input$subset_ortho, {

      gene_data <- gene_choice$df
      KO_level <- input$KO_level
      KO_group <- input$KO_group
      gene_data <- subset(gene_data, gene_data[,KO_level]== KO_group)

      updateSelectizeInput(inputId = "Gene", choices = unique(gene_data[,"KO_NAME"]),selected = character(0))

      gene_choice$df <- gene_data

    })


    data <- reactive({

      req(length(gene_choice$df) > 0)
      req(input$Gene)

      sub_genes <- input$Gene

      sub_genes <- gsub(pattern = "-", replacement = "_", sub_genes)

      sub_genes <- gsub(pattern = "+", replacement = "", sub_genes)

      sub_genes <- gsub(pattern = "\\(", replacement = "", sub_genes)

      sub_genes <- gsub(pattern = "\\)", replacement = "", sub_genes)

      sub_genes <- gsub(pattern = '^[^_]*_', replacement = "", sub_genes)

      data <- readRDS("inst/extdata/tidy_countdata_with_taxa_HF_Strep.rds")

      data <- subset(data, NAME %in% sub_genes)

      taxa <- input$taxa_level

      data[,taxa] <- gsub(pattern = "-", replacement = "_", x = data[,taxa])

      data[,taxa] <- gsub(pattern = " ", replacement = "_", x = data[,taxa])





      return(data)
    })

    #dynamic UI that only appears once fill_level inputs are selected
    output$myPanel <- renderUI({
      req(length(input$taxa_level)==1)
      req(data())

      taxa <- input$taxa_level
      df <- data()
      lev <- sort(unique(df[,taxa])) # sorting so that "things" are unambigious
      cols <- gg_fill_hue(length(lev))

      # New IDs "colX1" so that it partly coincide with input$select...
      lapply(seq_along(lev), function(i) {
        colourpicker::colourInput(inputId = NS(id,paste0("col", lev[i])),
                                  label = paste0("Choose Color for ", lev[i]),
                                  value = cbind(cols[i])
        )

      })
    })

    RA_plot <- reactive({

      req(data)
      req(length(input$taxa_level)==1)

      taxa <- input$taxa_level

      cols <- paste0("c(", paste0("input$col", sort(unique(data()[,taxa])), collapse = ", "), ")")
      # print(cols)
      cols <- eval(parse(text = cols))

      plot <- plot_gene_taxa_RAs(metagenome_data = data(), genes = input$Gene, taxa_level = input$taxa_level,
                                 font_size = input$font_size, cols = cols, border = input$border_size)

      return(plot)

    })


    output$RA_plot <- renderPlot(width = function() input$width, height = function() input$height,res = 96,{

      RA_plot()

    })


    output$download_plot <- downloadHandler(
      filename = function() {
        paste("RA_dotplot", input$extension, sep = ".")
      },

      content = function (file) {
        if (grepl(".pptx", file)==TRUE){

          doc <-  read_pptx()
          doc <- add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- ph_with(doc, value = RA_plot(), location = ph_location_type(type = "body"))
          print(doc, file)

        }else {

          ggplot2::ggsave(file, RA_plot(), device = input$extension, width = input$figure_width, height = input$figure_height, units = "cm", dpi = 300,
                 limitsize = FALSE)

        }

      }
    )

    #help button pop up text box
    observeEvent(input$help,{
      showModal(modalDialog(
        title = "Help",
        HTML("
        The data in the plots are read counts showing the amount of reads for the gene by taxonomic groups across samples, and the plot is faceted by treatment group.<br>
       <br>
       Perform the following steps to create a plot.<br>
       <br>
        1. Choose criteria to subset genes by. Options are significance, fold change, and KO groups (Optional).<br>
        <br>
        2. Select genes to plot from the drop down menu. Multiple genes can be plotted at once.<br>
        <br>
        The plot can be downloaded as either a PNG, PDF, SVG, or PPTX file."
        ),easyClose = TRUE))
    })





  })
}
#' Produces Gene Counts by Taxa Plot
#' @param id module id
#' @importFrom shiny NS plotOutput
RAOutput <- function(id){
  ns <- NS(id)

  plotOutput(ns("RA_plot"),width = 500, height = 1000)



}
#' Produces Color Palette Widget
#' @param id module id
#' @importFrom shiny NS uiOutput
RA_timepoint_color <- function(id){

  ns <- NS(id)
  uiOutput(ns("myPanel"))

}




