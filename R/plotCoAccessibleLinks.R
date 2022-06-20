#' plotCoAccessibleLinks
#' 
#' This function plots distribution of co-accessibility scores and percent accessible cells/aggregates.
#'
#' @param CoAccessibility GRanges object with co-accessible links from RWireX::addCoAccessibility and ArchR::getCoAccessibility.
#' @param Background List with background co-accessibility information from RWireX::getBackgroundCoAccessibility.
#' @param main Title for plot.
#' @return ggPlot.
#' @export


plotCoAccessibleLinks <- function(
    CoAccessibility = NULL,
    Background = NULL,
    main = NULL
    ){
    
    ArchR:::.validInput(input = CoAccessibility, name = "CoAccessibility", valid = c("GRanges"))
    ArchR:::.validInput(input = Background, name = "Background", valid = c("list"))
    ArchR:::.validInput(input = main, name = "main", valid = c("character"))
        
    coacc <- CoAccessibility@elementMetadata %>% as.data.frame(.)
    coacc$PercAccess <- apply(coacc[c("PercAccess1", "PercAccess2")], 1, mean)
    
    scatter_plot <- ggPoint(x = coacc$correlation, y = coacc$PercAccess, 
                            colorDensity = FALSE, pal = colorRampPalette(c("#f0f0f0", "black"))(100)[20:100], rastr = TRUE, legendSize = 8, baseSize = 16, labelSize = 8, 
                            xlabel = "Co-accessibility score", ylabel = "Percent accessible cells/aggregates", xlim = c(-1.1,1.1), ylim = c(-5,105)) +
                            xlim(c(-1.1,1.1)) + ylim(c(-0.1,100.1)) + ggtitle(main) +
                            geom_vline(xintercept= c(-Background$BackgroundCutoff, Background$BackgroundCutoff), linetype = 2, color = "red")
    
    return(ggExtra::ggMarginal(scatter_plot))
}