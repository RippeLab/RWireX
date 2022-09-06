#' plotCoAccessibleLinks
#' 
#' This function plots distribution of co-accessibility scores and percent accessible cells/aggregates.
#'
#' @param CoAccessibility GRanges object with co-accessible links or dataframe with element metadata from RWireX::addCoAccessibility and ArchR::getCoAccessibility.
#' @param Background List with background co-accessibility information from RWireX::getBackgroundCoAccessibility.
#' @param title Title for plot.
#' @param fast If TRUE, scattermore is used for fast plotting, but the detail is reduced.
#' @return ggPlot.
#' @export


plotCoAccessibleLinks <- function(
    CoAccessibility = NULL,
    Background = NULL,
    
    title = NULL,
    fast = TRUE
    ){
    
    ArchR:::.validInput(input = CoAccessibility, name = "CoAccessibility", valid = c("GRanges", "dataframe"))
    ArchR:::.validInput(input = Background, name = "Background", valid = c("list", "null"))
    ArchR:::.validInput(input = title, name = "title", valid = c("character"))
    ArchR:::.validInput(input = fast, name = "fast", valid = c("boolean"))
    
    if (is(CoAccessibility, "GRanges")){
      CoAccessibility <- CoAccessibility@elementMetadata %>% as.data.frame(.)
    }
    CoAccessibility$PercAccess <- apply(CoAccessibility[c("PercAccess1", "PercAccess2")], 1, mean)

    if (fast){
      stopifnot(requireNamespace("scattermore"))
      stopifnot(require(scattermore))
      #Have not figured out yet how to include aestethics into scattermost plot
      #For now use it as a fast analysis tool
      #ToDo: implement parameters for plotting
      scatter_plot = ggplot() + 
                            geom_scattermost(cbind(CoAccessibility$correlation,CoAccessibility$PercAccess), 
                            pixels=c(1000,1000), interpolate=TRUE) + 
                            xlab("Co-accessibility score")+ylab("Percent accessible cells/aggregates")+
                            xlim(c(-1.1,1.1)) + ylim(c(-0.1,100.1)) + ggtitle(title)
  
    } else{    

      scatter_plot = ggPoint(x = CoAccessibility$correlation, y = CoAccessibility$PercAccess, 
                            colorDensity = FALSE, pal = colorRampPalette(c("#f0f0f0", "black"))(100)[20:100], 
                            rastr = TRUE, legendSize = 8, baseSize = 16, labelSize = 8, 
                            xlabel = "Co-accessibility score", ylabel = "Percent accessible cells/aggregates", 
                            xlim = c(-1.1,1.1), ylim = c(-5,105)) +
                            xlim(c(-1.1,1.1)) + ylim(c(-0.1,100.1)) + ggtitle(title)
    }
    
    if (!is.null(Background)){
      scatter_plot = scatter_plot +
        geom_vline(xintercept= c(-Background$BackgroundCutoff, Background$BackgroundCutoff), linetype = 2, color = "red")
    }
    return(scatter_plot)
}