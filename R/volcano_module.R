#' Produces UI for Volcano Plot Module
#' @param id module id
#' @importFrom shiny actionButton HTML numericInput selectInput selectizeInput sliderInput radioButtons downloadButton tagList NS

volcano_ui <- function(id){

  ns <- NS(id)
  #help button UI
  tagList(
    actionButton(ns("help"), label = "Help/About this Plot", icon = icon("question-circle", lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Create a Volcano Plot</h4><br><br>"),
    HTML("<br><br><h5>Choose Data to Plot</h5><br><br>"),
    selectInput(ns("fill_variable"), "Select Orthology Metadata Category", choices = c("KO_Class", "KO_Subclass_1", "KO_Subclass_2")),
    #choices are selected server side based off observed event of fill_variable
    selectizeInput(ns("fill_levels"), "Select Genes to Plot", choices = NULL,multiple = TRUE),
    #could reduce this to a function, might be useful for other modules
    actionButton(ns("exclude_Others"), "Create Plot with Selected Genes",icon = icon("chart-bar"),
                 style="color: #fff; background-color: #0694bf; border-color: #013747"),
    actionButton(ns("include_Others"), "Create Plot with All Genes",icon = icon("chart-bar"),
                 style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h5>Adjust Plot Parameters</h5><br><br>"),
    actionButton(ns("rev_x"), "invert foldchange",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    #could reduce to a function using a dataframe of values for each input and a functional. sliders for plot dimensions etc.
    sliderInput(ns("height"), "Plot Height", min = 100, max = 1500, value = 500),
    sliderInput(ns("width"), "Plot Width", min = 100, max = 1500, value = 1000),
    sliderInput(ns("size"), "Point Size", min = 0.5, max = 10, value = 2),
    sliderInput(ns("sig_line_size"), "Significance Line Size", min = 0.5, max = 5, value = 1.5),
    sliderInput(ns("font"), "Font Size", min = 5, max = 25, value = 9),
    sliderInput(ns("border_size"), "Border Size", min = 0.5, max = 5, value = 1.5),
    HTML("<br><br><h5>Download Plot</h5><br><br>"),
    radioButtons(ns("extension"), "Save As:",
                 choices = c("png", "pdf", "svg", "pptx"), inline = TRUE),
    numericInput(ns("figure_height"), label = "Figure Height(cm)", value = 5),
    numericInput(ns("figure_width"), label = "Figure Width(cm)", value = 5),
    downloadButton(ns("download_plot"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747")

  )

}
#' Produces Server for Volcano Plot Module
#' @param id module id
#' @importFrom shiny moduleServer observeEvent reactiveValues updateSelectizeInput req updateSelectInput reactive renderUI downloadHandler showModal modalDialog renderPlot brushedPoints
#' @importFrom thematic thematic_shiny
#' @importFrom DT datatable renderDataTable
#' @importFrom officer read_pptx add_slide ph_with ph_location_type
#' @importFrom ggplot2 ggplot aes geom_point geom_hline theme_bw theme element_blank element_text element_rect xlim scale_color_manual coord_cartesian guides guide_legend ggsave
volcano_server <- function(id){

  moduleServer(id, function(input, output, session){
    thematic::thematic_shiny()

    #a reactive to update the data UI based on which data type  is selected as an input
    data <- reactive({

      data_stats <- readRDS("inst/extdata/HF_Strep_After_vs_Before.rds")
      data_stats[,sapply(data_stats,is.character)] <- as.data.frame(lapply(data_stats[,sapply(data_stats, is.character)],
                                                                           function(x) gsub(pattern = " ", replacement = "_", x = x )))
      data_stats$KO_Subclass_2 <- gsub(pattern = "^.+?(?=_).", replacement = "", x = data_stats$KO_Subclass_2, perl = TRUE)
      data_stats$KO_Subclass_2 <- gsub(pattern = "\\[PATH.*", replacement = "", x = data_stats$KO_Subclass_2)
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
      updateSelectInput(inputId = "fill_levels", choices = c(choices, "Other"))
    })

    #a reactive to make sure the table output from brush is correct
    calc_yaxis <- reactive({
      if (is.null(dat$df)) return()
      calc_data <- dat$df


      calc_data <- calc_data[,!names(calc_data) %in% c("baseMean", "pvalue", "stat","ENTRY",
                                                       "NAME")]

      calc_data <- data.frame(log2FoldChange = calc_data$log2FoldChange, padj = calc_data$padj, KO_NAME = calc_data$KO_NAME,
                              KO_Class = calc_data$KO_Class, KO_Subclass_1 = calc_data$KO_Subclass_1,
                              KO_Subclass_2 = calc_data$KO_Subclass_2, lfcSE = calc_data$lfcSE)

      calc_data$`-log10(padj)` <- -log10(calc_data$padj)

      return(calc_data)

    })

    #a reactive value meant to change the plot output depending on which create a plot action button is pressed
    #this works by either exluding data that is not in the fill_levels input, or changing exluded data to "Other" to include it
    dat <- reactiveValues(df = NULL)


    observeEvent(input$include_Others, {

      dat$df <- data()
      dat$df[,input$fill_variable] <- sapply(dat$df[, input$fill_variable], function(x) replace(x, !x %in% input$fill_levels, "Other"))
      dat$df <- dat$df[!duplicated(dat$df[,c("ENTRY")]),]
      choices <- unique(fill_choice())
      updateSelectInput(inputId = "fill_levels", choices = c(choices, "Other"), selected = c(input$fill_levels, "Other"))

    })

    observeEvent(input$exclude_Others, {

      dat$df <- data()
      dat$df[,input$fill_variable] <- sapply(dat$df[, input$fill_variable], function(x) replace(x, !x %in% input$fill_levels, NA))
      dat$df <- dat$df[!is.na(dat$df[,input$fill_variable]),]
      dat$df <- dat$df[!duplicated(dat$df[,c("ENTRY")]),]

    })

    #observation for the invert fold change action button, for if someone wants to swap numerator and denominator contrasts
    observeEvent(input$rev_x, {

      if (is.null(dat$df)) return()
      dat$df[,"log2FoldChange"] <- dat$df[,"log2FoldChange"] * -1


    })

    #reactive value for plot zoom
    ranges <- reactiveValues(x = NULL, y = NULL)

    #dynamic UI that only appears once fill_level inputs are selected
    output$myPanel <- renderUI({
      lev <- sort(unique(input$fill_levels)) # sorting so that "things" are unambigious
      cols <- gg_fill_hue(length(lev))

      # New IDs "colX1" so that it partly coincide with input$select...
      lapply(seq_along(lev), function(i) {
        colourpicker::colourInput(inputId = NS(id,paste0("col", lev[i])),
                                  label = paste0("Choose Color for ", lev[i]),
                                  value = cbind(cols[i])
        )

      })
    })

    #reactive which creates the volcano plot. is dependent on data(), and also the reactive value dat(aka plot button being pressed)
    #takes inputs from fill_levels color, font, width, height thus requires these as well
    volcano_output <- reactive(
      {
        if (is.null(dat$df)) return()
        cols <- paste0("c(", paste0("input$col", sort(input$fill_levels), collapse = ", "), ")")
        # print(cols)
        cols <- eval(parse(text = cols))
        # print(cols)
        # To prevent errors
        req(input$fill_levels)
        req(length(cols) == length(input$fill_levels))
        req(data())

        #req(input$exclude_Others | input$include_Others)


        variable <- dat$df[,input$fill_variable]

        ggplot2::ggplot(dat$df, ggplot2::aes(x = log2FoldChange, y = -log10(padj), color = unlist(variable)
                           #,color = unlist(variable),
                           #shape = unlist(variable)
        )) +
          ggplot2::geom_point(size = input$size) +
          ggplot2::geom_hline(ggplot2::aes(yintercept = -log10(0.05)),linetype = "dashed", size = input$sig_line_size, color = "black") +
          ggplot2::xlim((0-abs(max(data()[,"log2FoldChange"]))),(0+max(abs(data()[,"log2FoldChange"]))))+
          ggplot2::scale_color_manual(values = cols) +
          ggplot2::coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = TRUE) +
          #scale_color_manual(values = c(rep("black", length(cols))))+
          #scale_shape_manual(values = c(rep(21, length(cols))))+
          ggplot2::theme_bw()+
          ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_rect(size = input$border_size)) +
          ggplot2::theme(legend.title = ggplot2::element_blank()) +
          ggplot2::theme(axis.text = ggplot2::element_text(size = input$font)) +
          ggplot2::theme(axis.title = ggplot2::element_text(size = input$font)) +
          ggplot2::theme(legend.text = ggplot2::element_text(size = input$font), legend.position = "top")+
          ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))

      }
    )

    #here so that height and width sliders function properly
    output$volcano <- renderPlot(width = function() input$width, height = function() input$height,res = 96,{

      volcano_output()

    })

    #observer for plot zoom
    observeEvent(input$plot_dblclick, {
      brush <- input$plot_brush
      if (!is.null(brush)) {
        ranges$x <- c(brush$xmin, brush$xmax)
        ranges$y <- c(brush$ymin, brush$ymax)

      } else {
        ranges$x <- NULL
        ranges$y <- NULL
      }
    })
    #creates table describing selected datapoints on plot
    output$data <- DT::renderDataTable({
      req(input$plot_brush)
      dt <- brushedPoints(calc_yaxis(), input$plot_brush)

      DT::datatable(data = dt, options = list(scrollX = TRUE), width = 1800)
    })

    #download handler for plot, the boolean after content is required for pptx to be an option, as it is not covered by ggsave
    output$download_plot <- downloadHandler(
      filename = function() {
        paste("volcano_plot", input$extension, sep = ".")
      },

      content = function (file) {
        if (grepl(".pptx", file)==TRUE){

          doc <-  officer::read_pptx()
          doc <- officer::add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- officer::ph_with(doc, value = volcano_output(), location = officer::ph_location_type(type = "body"))
          print(doc, file)

        }else {

          ggplot2::ggsave(file, volcano_output(), device = input$extension,width = input$figure_width, height = input$figure_height, units = "cm", dpi = 300,
                 limitsize = FALSE)

        }

      }
    )

    #help button pop up text box
    observeEvent(input$help,{
      showModal(modalDialog(
        title = "Help",
        HTML("
        The data in the plot are derived from a paired DESeq2 wald test.
        One group is HF diet before receiving streptomycin, the other group is HF diet 4 weeks after receiving streptomycin. Samples with a negative log2FoldChange are lower in the
        mice 4 weeks after streptomycin treatment , samples with a positive log2FoldChange are higher after treatment.<br>
       <br>
       Perform the following steps to create a plot.<br>
       <br>
        1. Select a KO metadata category.<br>
        <br>
        2. Select KO group to color points by.<br>
        <br>
        3. Click one of the create a plot buttons to make the plot. Choosing to only use the KO groups selected will only plot genes in those groups.
        Choosing to include all genes will plot those groups and then plot all other genes as 'other', enabling you to select colors.<br>
        <br>
        The plot can be downloaded as either a PNG, PDF, SVG, or PPTX file."
        ),easyClose = TRUE))
    })




  })
}
#' Creates Volcano Plot
#' @param id module id
#' @importFrom shiny NS plotOutput brushOpts
volcanoOutput <- function(id){
  ns <- NS(id)

  plotOutput(ns("volcano"),dblclick = ns("plot_dblclick"), brush = brushOpts(
    id = ns("plot_brush"),
    resetOnNew = TRUE
  ),width = 1000, height = 500)



}
#' Creates Data Table
#' @param id module id
#' @importFrom DT dataTableOutput
#' @importFrom shiny NS
vol_table <- function(id){
  ns <- NS(id)
  DT::dataTableOutput(ns("data"))

}
#' Creates color Palette
#' @param id module id
#' @importFrom shiny NS uiOutput
color_output <- function(id){

  ns <- NS(id)
  uiOutput(ns("myPanel"))

}
