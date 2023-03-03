#' Produces UI for Gene Dotplot Module
#' @param id module id
#' @importFrom shiny actionButton HTML numericInput selectInput selectizeInput sliderInput radioButtons downloadButton tagList NS checkboxInput


dotplot_ui <- function(id){

  ns <- NS(id)
  #help button UI
  tagList(
    actionButton(ns("help"), label = "Help/About this Plot", icon = icon("question-circle", lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h4>Create Gene Dotplots</h4><br><br>"),
    actionButton(ns("start"), label = "Load Dataset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    HTML("<br><br><h5>Choose Data to Plot</h5><br><br>"),
    numericInput(inputId = ns("sig_level"),label = "Select Significance Threshold", value = 0.05),
    actionButton(ns("subset_sig"), "Subset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    selectInput(ns("KO_level"), label = "Select Orthology Metadata Category", choices = c("KO_Class", "KO_Subclass_1", "KO_Subclass_2")),
    selectInput(ns("KO_group"), label = "Select Gene Group", choices = NULL, multiple = FALSE),
    actionButton(ns("subset_ortho"), "Subset",style="color: #fff; background-color: #0694bf; border-color: #013747"),
    selectizeInput(ns("Gene"), "Select Gene to Plot", choices = NULL,multiple = TRUE),
    #could reduce to a function using a dataframe of values for each input and a functional. sliders for plot dimensions etc.
    HTML("<br><br><h5>Adjust Plot Parameters</h5><br><br>"),
    checkboxInput(ns("checkbox"), "Toggle Point Labels"),
    sliderInput(ns("height"), "Plot Height", min = 100, max = 1500, value = 300),
    sliderInput(ns("width"), "Plot Width", min = 100, max = 1500, value = 1500),
    sliderInput(ns("size"), "Point Size", min = 1, max = 10, value = 2, step = 0.25),
    sliderInput(ns("font_size"), "Font Size",min = 5, max = 30, value = 14),
    sliderInput(ns("border_size"), "Border Size", min = 0.5, max = 5, value = 1.5),
    HTML("<br><br><h5>Download Plot</h5><br><br>"),
    radioButtons(ns("extension"), "Save As:",
                 choices = c("pdf", "png","svg", "pptx"), inline = TRUE),
    numericInput(ns("fig_width"), "Base Figure Width", value = 17, min = 1, max = 30),
    numericInput(ns("fig_height"), "Base Figure Height", value = 4, min = 1, max = 30),
    downloadButton(ns("download_plot"), "Save Plot",style="color: #fff; background-color: #0694bf; border-color: #013747")

  )

}

#' Produces Server for Gene Dotplot Module
#' @param id module id
#' @importFrom shiny moduleServer observeEvent reactiveValues updateSelectizeInput req updateSelectInput reactive renderUI downloadHandler showModal modalDialog renderPlot renderDataTable
#' @importFrom thematic thematic_shiny
#' @importFrom stats ave qnorm
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme element_text element_blank element_rect position_jitter stat_summary scale_color_manual ylab coord_flip ggtitle guides geom_errorbar element_line guide_legend
#' @importFrom ggrepel geom_label_repel
#' @importFrom cowplot plot_grid save_plot
#' @importFrom colourpicker colorInput
#' @importFrom DT datatable
#' @importFrom officer read_pptx add_slide ph_with ph_location_type
dotplot_server <- function(id){

  moduleServer(id, function(input, output, session){
    thematic::thematic_shiny()

    gene_choice <- reactiveValues(df = NULL)

    observeEvent(input$start, {


      gene_choice$df <- readRDS("inst/extdata/HF_Strep_After_vs_Before.rds")


      updateSelectizeInput(inputId = "Gene", choices = unique(gene_choice$df[,"KO_NAME"]))



    })

    observeEvent(input$subset_sig,{

      sig_level <- input$sig_level

      gene_data <- gene_choice$df

      gene_data<- gene_data[gene_data$padj <= sig_level,]

      updateSelectizeInput(session,inputId = "Gene", choices = unique(gene_data[,"KO_NAME"]), server = TRUE)

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

      updateSelectizeInput(session,inputId = "Gene", choices = unique(gene_data[,"KO_NAME"]), server = TRUE)

      gene_choice$df <- gene_data

    })


    data <- reactive({

      req(length(input$Gene)>= 1)

      data <- readRDS("inst/extdata/tidy_normalized_countdata_HF_Strep.rds")

      sig_data <- readRDS("inst/extdata/HF_Strep_After_vs_Before.rds")

      sample_data <- readRDS("inst/extdata/sample_data_HF_strep.rds")

      data$KO_NAME <- paste0(data$ENTRY, "_", data$NAME)

      data <- merge(data, sample_data, by = "Sample")

      sub_genes <- input$Gene

      sub_genes <- gsub(pattern = "_.*", replacement = "", sub_genes)

      data <- data[data$ENTRY %in% sub_genes,]

      data[,"ENTRY_Sample"] <- paste0(data[,"ENTRY"],"_", data[,"Sample"])
      data<- data[!duplicated(data[,c("ENTRY_Sample")]),]


      data$padj <- sig_data$padj[match(data$ENTRY, sig_data$ENTRY)]

      return(data)
    })

    #dynamic UI that only appears once fill_level inputs are selected
    output$myPanel <- renderUI({
      lev <- sort(unique(data()$Timepoint)) # sorting so that "things" are unambigious
      cols <- gg_fill_hue(length(lev))

      # New IDs "colX1" so that it partly coincide with input$select...
      lapply(seq_along(lev), function(i) {
        colourpicker::colourInput(inputId = NS(id,paste0("col", lev[i])),
                                  label = paste0("Choose Color for ", lev[i]),
                                  value = cbind(cols[i])
        )

      })
    })

    dot_plot <- reactive({

      req(data)

      df <- data()
      df$Abundance <- log2(df$Abundance)
      df$max_abund <- ave(seq_len(nrow(df)), df$NAME,
                                          FUN = function(i) df$Abundance[i][which.max(df$Abundance[i])])
      cols <- paste0("c(", paste0("input$col", sort(unique(data()$Timepoint)), collapse = ", "), ")")
      # print(cols)
      cols <- eval(parse(text = cols))

      df_2 <- data.frame(NAME = unique(df$NAME))
      df_2$padj <- df$padj[match(df_2$NAME, df$NAME)]
      df_2$ENTRY <- df$ENTRY[match(df_2$NAME, df$NAME)]
      df_2$Abundance <- df$max_abund[match(df_2$ENTRY, df$ENTRY)]
      df_2$Abundance <- df_2$Abundance + 0.5

      df$Timepoint[df$Timepoint=="T1"] <- "HF_before_Str"
      df$Timepoint[df$Timepoint=="T2"] <- "HF_4_weeks_after_Str"

      p1 <- ggplot2::ggplot(df, ggplot2::aes(x=NAME, y=Abundance, color = Timepoint)) +

        ggplot2::geom_point(position=ggplot2::position_jitter(w=0.1,h=0), size = input$size,shape = 21, fill = "white", stroke = 1.5) +
        {if (input$checkbox) ggrepel::geom_label_repel(ggplot2::aes(label = Sample),box.padding   = 0.75,
                                              point.padding = 0.75)} +
        ggplot2::geom_text(data = df_2,ggplot2::aes(label = ifelse(padj<=0.05, "*", ""), group = NAME),color = "black", vjust = .7, size = 10)+
        #stat_summary(fun.y = "mean", geom = "point", size = 5, color = "black") +
        ggplot2::stat_summary(fun.y = "mean", geom = "segment", ggplot2::aes(xend=..x.. -0.45, yend = ..y..),size = 2) +
        ggplot2::stat_summary(fun.y = "mean", geom = "segment", ggplot2::aes(xend=..x.. +0.45, yend = ..y..), size = 2) +
        ggplot2::scale_color_manual(values = cols) +
        ggplot2::ylab("log2 Normalized Abundance") +
        ggplot2::coord_flip() +
        ggplot2::ggtitle("") +
        ggplot2::theme_bw()+
        ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_rect(size = input$border_size), axis.text = ggplot2::element_text(size = input$font_size),
              axis.title = ggplot2::element_text(size = input$font_size),title = ggplot2::element_text(size = input$font_size-1.5),
              legend.text = ggplot2::element_text(size = input$font_size - 2),
              legend.title = ggplot2::element_blank(),
              axis.title.y = ggplot2::element_blank(), legend.position = c(0.3, 1.10),
              legend.key.size = ggplot2::unit(0.2, "cm"),axis.ticks=ggplot2::element_line(size=1.2), axis.ticks.length=ggplot2::unit(.25, "cm"))+
        ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1))

      data <- gene_choice$df

      data$KO_NAME <- paste0(data$ENTRY, "_", data$NAME)

      data <- data[data$KO_NAME %in% input$Gene,]

      p2 <- ggplot2::ggplot(data, ggplot2::aes(x = NAME, y = log2FoldChange)) + ggplot2::geom_point(size = 5, color = "black") +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = log2FoldChange - qnorm(0.025)*lfcSE, ymax = log2FoldChange + qnorm(0.025)*lfcSE), size = 1.25, color = "black") +
        ggplot2::ggtitle("95% Confidence Interval") +
        ggplot2::coord_flip() +
        ggplot2::theme_bw()+
        ggplot2::theme(panel.grid = ggplot2::element_blank(), panel.border = ggplot2::element_rect(size = input$border_size), axis.text = ggplot2::element_text(size = input$font_size),
              axis.title = ggplot2::element_text(size = input$font_size), title = ggplot2::element_text(size = input$font_size - 1.5),
              axis.title.y = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),axis.ticks=ggplot2::element_line(size=1.2),axis.ticks.length=ggplot2::unit(.20, "cm"))

      p_all <- cowplot::plot_grid(p1, p2, rel_widths = c(1,0.5), ncol = 2)

      return(p_all)

      })


    output$dotplot <- renderPlot(width = function() input$width, height = function() input$height,res = 96,{

      dot_plot()

    })

    #creates table describing selected datapoints on plot
    output$dotplot_table <- DT::renderDataTable({
      req(data())

      df <- data()

      DT::datatable(data = df, options = list(scrollX = TRUE), width = 1800)
    })

    output$download_plot <- downloadHandler(
      filename = function() {
        paste("dotplot", input$extension, sep = ".")
      },

      content = function (file) {
        if (grepl(".pptx", file)==TRUE){

          doc <-  officer::read_pptx()
          doc <- officer::add_slide(doc, layout =  'Title and Content', master = 'Office Theme')
          doc <- officer::ph_with(doc, value = dot_plot(), location = officer::ph_location_type(type = "body"))
          print(doc, file)

        }else {

          cowplot::save_plot(file, dot_plot(), device = input$extension, base_width = input$fig_width, base_height = input$fig_height)

        }

      }
    )

    #help button pop up text box
    observeEvent(input$help,{
      showModal(modalDialog(
        title = "Help",
        HTML("The plot
        on the left shows log2 transformed, median of ratio normalized read counts for each gene by sample, with bars denoting group means and an asterisk
        denoting a signficant difference between groups (if applicable) following FDR correction.
        The plot on the right is a 95% confidence interval of the log2FoldChange between groups, calculated with the following equation:
        log2FoldChange + qnorm(0.025)*lfcSE and log2FoldChange - qnorm(0.025)*lfcSE, for upper and lower bounds respectively. <br>
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
#' Creates Dotplot
#' @param id module id
#' @importFrom shiny plotOutput NS

dotplotOutput <- function(id){
  ns <- NS(id)

  plotOutput(ns("dotplot"),width = 1000, height = 500)



}
#' Creates Color Palette
#' @param id module id
#' @importFrom shiny uiOutput NS
timepoint_color <- function(id){

  ns <- NS(id)
  uiOutput(ns("myPanel"))

}
#' Creates Table Output
#' @param id module id
#' @importFrom shiny NS
#' @importFrom DT dataTableOutput
dot_table <- function(id){
  ns <- NS(id)
  DT::dataTableOutput(ns("dotplot_table"))

}


