#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#library(shiny)
#library(shinyFeedback)
#library(bslib)
#library(shinyWidgets)
#library(ggplot2)
#library(DT)
#library(officer)
#library(thematic)
#library(colourpicker)
#library(reshape2)
#library(RColorBrewer)
#library(gplots)
#library(gridExtra)
#library(cowplot)
#library(ggrepel)
#library(ggplotify)
#library(promises)
#library(future)



#my_theme <- bslib::bs_theme(bootswatch = "flatly")
#' The shiny app
#' @importFrom shiny fluidPage titlePanel tags radioButtons sidebarLayout sidebarPanel tabsetPanel tabPanel mainPanel observe session shinyApp
#' @importFrom bslib bs_theme bs_theme_update
#' @importFrom future plan multisession

#future::plan(future::multisession)
myApp <- function(...){
  my_theme <- bslib::bs_theme(bootswatch = "flatly")
  ui <- fluidPage(
    titlePanel("Metagenomics of Sorbitol Intolerant Mice"),
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"),
    theme = my_theme,
    radioButtons("current_theme", "App Theme:", c("Light" = "flatly", "Dark" = "darkly"), inline = TRUE),
    actionButton("about_data", label = "About this Dataset", icon = icon("question-circle", lib = "font-awesome"),style="color: #fff; background-color: #0694bf; border-color: #013747"),
    sidebarLayout(
      sidebarPanel(     tags$head(tags$style(type="text/css", "
             #loadmessage {
               position: fixed;
               top: 0px;
               left: 0px;
               width: 100%;
               padding: 5px 0px 5px 0px;
               text-align: center;
               font-weight: bold;
               font-size: 100%;
               color: #000000;
               background-color: #CCFF66;
               z-index: 105;
             }
          ")),conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                          tags$div("Loading...",id="loadmessage")),
        tabsetPanel(
          tabPanel("Plot Heatmaps", heatmap_ui("myModule1")),
          tabPanel("Plot Volcano", volcano_ui("myModule2"), color_output("myModule2")),
          tabPanel("Plot Dotplot", dotplot_ui("myModule3"), timepoint_color("myModule3")),
          tabPanel("Plot Gene Counts by Taxa", RA_ui("myModule4"), RA_timepoint_color("myModule4")))
         ),
      mainPanel(tabsetPanel(
          tabPanel("Heatmap", heatmapOutput("myModule1")),
          tabPanel("Volcano Plot", verticalLayout(volcanoOutput("myModule2"), vol_table("myModule2"))),
          tabPanel("Dotplot", verticalLayout(dotplotOutput("myModule3"), dot_table("myModule3"))),
          tabPanel("Gene Counts By Taxa", RAOutput("myModule4"))
        ), width = 6)))

  server <- function(input, output, session) {

    observe({
      # Make sure theme is kept current with desired
      session$setCurrentTheme(
        bslib::bs_theme_update(my_theme, bootswatch = input$current_theme)
      )
    })

    #help button pop up text box
    observeEvent(input$about_data,{
      showModal(modalDialog(
        title = "Help",
        HTML("Sample metagenomes were sequenced using an Illumina Novaseq 6000, paired end 150bp.<br>
        <br>
        Metagenome quality control, and binning were performed using the DOE JGI Metagenome workflow https://journals.asm.org/doi/10.1128/mSystems.00804-20. For assembly, metaspades v3.15.2 was
        used to perform a co-assembly to create higher quality bins. Following assembly, BBMAP was used to map each of the sequencing fastq files to the coassembly contigs. Contigs and coverage information
        were uploaded to IMG for annotation and binning.<br>
        <br>
        Count matrices of reads were generated using medium to high quality bins, and gene names were gathered from the KEGG API, along with brite hierarchy metadata.<br>
        <br>
        The experiment was a paired study with c57B6j mice, where mice were given a treatment of the broad spectrum antibiotic streptomycin and a diet switch to a high fat diet. Samples shown are before treatment
             and 4 weeks following treatment"
        ), easyClose = TRUE))
    })

    heatmap_server("myModule1")

    volcano_server("myModule2")

    dotplot_server("myModule3")

    RA_server("myModule4")

  }
  shinyApp(ui, server)

}

myApp()
