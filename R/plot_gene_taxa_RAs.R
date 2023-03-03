#'Plot Gene Counts by Taxonomic Groups
#' @param metagenome_data metagenomic data in long form with columns for genes counts, orthology groups,
#' sample data and taxonomic data
#' @param genes A string. The gene to plot
#' @param taxa_level A string, the taxonomic level to group by for plotting
#' @param font_size An integer. Plot font size
#' @param cols Colors for taxonomic groups, a character vector of length n, where n is the number of taxonomic levels
#' @param border An integer. Determines plot border thickness
#' @return A ggplot2 object
#' @importFrom stats aggregate as.formula
#' @importFrom ggplot2 ggplot geom_bar aes facet_grid scale_fill_manual ylab theme_bw theme element_text element_rect element_blank as_labeller
#' @importFrom spsComps shinyCatch


plot_gene_taxa_RAs <- function(metagenome_data, genes, taxa_level, font_size, cols, border){
  #create formula for aggregate
  formula <- stats::as.formula(paste0("Abundance", "~","Sample", "+", "Timepoint", "+", "NAME", "+", taxa_level))

  metagenome_data <- metagenome_data[!is.na(metagenome_data$Sample),]
  #sum counts by taxa
  if(length(metagenome_data[!is.na(metagenome_data[1]),])==0){

    spsComps::shinyCatch(stop("No taxonomic information for this gene"), blocking_level = "error")

  }

  if(eval(nrow(metagenome_data[metagenome_data[,taxa_level]=="None",]) == nrow(metagenome_data))==TRUE){

    spsComps::shinyCatch(stop("No taxonomic information for this gene at selected level"), blocking_level = "error")

  }


  metagenome_agg <-  stats::aggregate(formula, metagenome_data, sum)
  #create barplots
  ggplot2::ggplot() +
    ggplot2::geom_bar(data = metagenome_agg, mapping = ggplot2::aes(x = Sample, y = Abundance, fill = metagenome_agg[,taxa_level]),
             stat = "identity", position = "stack", color = "black", size = 1) +
    ggplot2::facet_grid(.~Timepoint, scales = "free", space = "free_x", labeller = ggplot2::as_labeller(c("Before_HF+Strep"="HF_before_Str", "After_HF+Strep" = "HF_4_weeks_after_Str"))) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::ylab(paste0(gsub(pattern = "_", replacement = " ", x = genes), "", "Abundance")) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),axis.text.x = ggplot2::element_text(angle = 90), axis.text = ggplot2::element_text(size = font_size),
          strip.text = ggplot2::element_text(size = font_size),legend.text = ggplot2::element_text(size = font_size), legend.position = "top",
          axis.title = ggplot2::element_text(size = font_size+1), legend.title = ggplot2::element_blank(), panel.border = ggplot2::element_rect(size = border))



}
