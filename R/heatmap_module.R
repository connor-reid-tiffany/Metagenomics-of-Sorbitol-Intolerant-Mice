#'Produces heatmap module ui
#' @param id module id
#' @importFrom shiny actionButton HTML numericInput selectInput selectizeInput sliderInput radioButtons downloadButton tagList NS checkboxInput

heatmap_ui <- function(id){

  ns <- NS(id)
  #help button UI
  tagList(
    actionButton(ns("help"), label = "Help/About this Plot", icon = icon("question-circle", lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Create a Heatmap</h4><br><br>"),
    HTML("<br><br><h5>Choose Data to Plot</h5><br><br>"),
    selectInput(ns("fill_variable"), "Select Orthology Metadata Category", choices = c("KO_Class", "KO_Subclass_1", "KO_Subclass_2")),
    #choices are selected server side based off observed event of fill_variable
    selectizeInput(ns("fill_levels"), "Select Genes to Plot", choices = NULL,multiple = FALSE),
    #could reduce this to a function, might be useful for other modules
    actionButton(ns("exclude_Others"), "Create Plot with Selected Genes",icon = icon("chart-bar"),
                 style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h5>Adjust Plot Parameters</h5><br><br>"),
    checkboxInput(ns("checkbox"), "Toggle Gene Labels"),
    checkboxInput(ns("col_dendrogram"), "Toggle Column Dendrogram"),
    #could reduce to a function using a dataframe of values for each input and a functional. sliders for plot dimensions etc.
    sliderInput(ns("height"), "Plot Height", min = 100, max = 10000, value = 1500),
    sliderInput(ns("width"), "Plot Width", min = 100, max = 10000, value = 1500),
    sliderInput(ns("font_size"), "Scale Y Axis Font Size", min = 0.5, max = 3, value = 1, step = 0.1),
    sliderInput(ns("font_size_col"), "Scale X Axis Font Size", min = 0.5, max = 3, value = 1, step = 0.1),
    HTML("<br><br><h5>Download Plot</h5><br><br>"),
    radioButtons(ns("extension"), "Save As:",
                 choices = c("png", "pdf", "svg", "pptx"), inline = TRUE),
    numericInput(ns("figure_height"), label = "Figure Height(cm)", value = 5),
    numericInput(ns("figure_width"), label = "Figure Width(cm)", value = 5),
    downloadButton(ns("download_plot"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747")

  )

}

#' Produces heatmap server
#' @param id module id
#' @importFrom shiny moduleServer reactive observeEvent reactiveValues renderPlot downloadHandler showModal modalDialog req updateSelectInput
#' @importFrom ggplot2 ggsave
#' @importFrom ggplotify as.grob as.ggplot
#' @importFrom gplots heatmap.2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom reshape2 dcast
#' @importFrom officer read_pptx add_slide ph_with ph_location_type
#' @importFrom thematic thematic_shiny
#' @importFrom promises future_promise

heatmap_server <- function(id){

  moduleServer(id, function(input, output, session){
    thematic::thematic_shiny()

    #a reactive to update the data UI based on which data type  is selected as an input
    data <- reactive({


      data_stats <- readRDS("inst/extdata//tidy_normalized_countdata_HF_Strep.rds")

      data_stats$Abundance <- log(1 + data_stats$Abundance)


      data_stats$Sample <- as.character(data_stats$Sample)

      data_stats[,sapply(data_stats,is.character)] <- as.data.frame(lapply(data_stats[,sapply(data_stats, is.character)],
                                                                           function(x) gsub(pattern = " ", replacement = "_", x = x )))
      data_stats$KO_Subclass_2 <- gsub(pattern = "^.+?(?=_).", replacement = "", x = data_stats$KO_Subclass_2, perl = TRUE)
      data_stats$KO_Subclass_2 <- gsub(pattern = "\\[PATH.*", replacement = "", x = data_stats$KO_Subclass_2)

      data_stats$KO_NAME <- paste0(data_stats$ENTRY, "_", data_stats$NAME)

      return(data_stats)

    })


    #a reactive to update the fill_levels UI based on which fill_variable choice is selected as an input
    fill_choice <- reactive({
      fill_data <- data()
      fill_variable <- fill_data[,input$fill_variable]
      return(fill_variable)

    })

    #an observation which dynamically updates the fill_levels input options when the fill_variable input is changed
    observeEvent(fill_choice(), {
      choices <- unique(fill_choice())
      updateSelectInput(inputId = "fill_levels", choices = c(choices))
    })

    #a reactive value meant to change the plot output depending on which create a plot action button is pressed
    #this works by either exluding data that is not in the fill_levels input, or changing exluded data to "Other" to include it
    dat <- reactiveValues(df = NULL)

    observeEvent(input$exclude_Others, {

      dat$df <- data()


      dat$df[,input$fill_variable] <- sapply(dat$df[, input$fill_variable], function(x) replace(x, !x %in% input$fill_levels, NA))
      dat$df <- dat$df[!is.na(dat$df[,input$fill_variable]),]
      dat$df <- dat$df[!is.na(dat$df[,"Abundance"]),]

      dat$df[,"ENTRY_Sample"] <- paste0(dat$df[,"ENTRY"],"_", dat$df[,"Sample"])
      dat$df <- dat$df[!duplicated(dat$df[,c("ENTRY_Sample")]),]






    })


    #reactive which creates the volcano plot. is dependent on data(), and also the reactive value dat(aka plot button being pressed)
    #takes inputs from fill_levels color, font, width, height thus requires these as well
    heatmap_output <- reactive(
      {
        if (is.null(dat$df)) return()

        req(input$fill_levels)
        req(data())


        plot_data <- dat$df[,c("KO_NAME", "Sample", "Abundance")]
        plot_data <- reshape2::dcast(plot_data, Sample ~ KO_NAME)
        rownames(plot_data) <- plot_data$Sample
        plot_data <- plot_data[,-1]


        breaks = seq(0,max(as.matrix(plot_data)),length.out=26)
        coul <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(25)

        if(input$checkbox){

          margins = c(5,0.3)

        } else{

          margins = c(5,13)

        }

        if(input$col_dendrogram){

          dendrogram <- c("none")

        } else{

          dendrogram <- c("col")
        }

        font_size <- input$font_size
        font_size_col <- input$font_size_col



        promises::future_promise({ggplotify::as.ggplot(ggplotify::as.grob(function() gplots::heatmap.2(x = t(plot_data),col = coul,scale = "none", dendrogram = dendrogram, trace = "none", density.info = "none",
                                                                                             breaks = breaks, margins =
                                                                                               margins, cexRow = font_size, cexCol = font_size_col)))
})


      }
    )

    #here so that height and width sliders function properly
    output$heatmap <- renderPlot(width = function() input$width, height = function() input$height,res = 96,{

      heatmap_output()

    })



    #download handler for plot, the boolean after content is required for pptx to be an option, as it is not covered by ggsave
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("heatmap", input$extension, sep = ".")
      },

      content = function (file) {
        if (grepl(".pptx", file)==TRUE){

          doc <-  officer::read_pptx()
          doc <- officer::add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- officer::ph_with(doc, value = heatmap_output(), location = officer::ph_location_type(type = "body"))
          print(doc, file)

        }else {

          ggplot2::ggsave(file, heatmap_output(), device = input$extension, width = input$figure_width, height = input$figure_height, units = "cm", dpi = 300,
                 limitsize = FALSE)

        }

      }
    )

    #help button pop up text box
    observeEvent(input$help,{
      showModal(modalDialog(
        title = "Help",
        HTML("
        The data shown in the heatmap were first transformed using the median of ratios method in DESeq2, and were then transformed by the glog for centering.
        The hierarchical clustering algorithm used to generate gene and sample dendrograms uses euclidean distance with a complete agglomeration method.
        Samples that end in 1 are mice on a HF diet before streptomycin treatment, while samples that end in 2 are mice on a HF diet 4 weeks after streptomycin treatment.<br>
       <br>
       Perform the following steps to create a plot.<br>
       <br>
        1. Choose a KO metadata category to subset the data (optional).<br>
        <br>
        2. Choose a Gene group within that category to subset the data (optional).<br>
        <br>
        3. Click either create a plot button. Create a plot with select genes produces a heatmap of the subsetted data. Create a plot with all genes
        makes a heatmap with all the genes included.<br>
        <br>
        4. The checkbox for gene labels can be toggled if there are too many genes for labels to work.
       The plot can be downloaded as either a PNG, PDF, SVG, or PPTX file."
        ),easyClose = TRUE))
    })




  })
}

#'Produces Heatmap
#' @param id module id
#' @importFrom shiny NS plotOutput

heatmapOutput <- function(id){
  ns <- NS(id)

  plotOutput(ns("heatmap"),width = 1000, height = 500)



}
