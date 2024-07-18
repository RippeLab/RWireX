#' Plot an ArchR Region Track
#' 
#' This function will plot the coverage at an input region in the style of a browser track. It allows for normalization of the signal
#' which enables direct comparison across samples.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param region A `GRanges` region that indicates the region to be plotted. If more than one region exists in the `GRanges` object,
#' all will be plotted. If no region is supplied, then the `geneSymbol` argument can be used to center the plot window at the
#' transcription start site of the supplied gene.
#' @param groupBy A string that indicates how cells should be grouped. This string corresponds to one of the standard or
#' user-supplied `cellColData` metadata columns (for example, "Clusters"). Cells with the same value annotated in this metadata
#' column will be grouped together and the average signal will be plotted.
#' @param useGroups A character vector that is used to select a subset of groups by name from the designated `groupBy` column in
#' `cellColData`. This limits the groups to be plotted.
#' @param plotSummary A character vector containing the features to be potted. Possible values include "bulkTrack" (the ATAC-seq signal),
#' "scTrack" (scATAC-seq signal), "featureTrack" (i.e. the peak regions), "geneTrack" (line diagrams of genes with introns and exons shown. 
#' Blue-colored genes are on the minus strand and red-colored genes are on the plus strand), and "loopTrack" (links between a peak and a gene).
#' @param sizes A numeric vector containing up to 3 values that indicate the sizes of the individual components passed to `plotSummary`.
#' The order must be the same as `plotSummary`.
#' @param features A `GRanges` object containing the "features" to be plotted via the "featureTrack". This should be thought of as a
#' bed track. i.e. the set of peaks obtained using `getPeakSet(ArchRProj))`. 
#' @param loops A `GRanges` object containing the "loops" to be plotted via the "loopTrack".
#' This `GRanges` object start represents the center position of one loop anchor and the end represents the center position of another loop anchor. 
#' A "loopTrack" draws an arc between two genomic regions that show some type of interaction. This type of track can be used 
#' to display chromosome conformation capture data or co-accessibility links obtained using `getCoAccessibility()`. 
#' @param geneSymbol If `region` is not supplied, plotting can be centered at the transcription start site corresponding to the gene symbol(s) passed here.
#' @param useMatrix If supplied geneSymbol, one can plot the corresponding GeneScores/GeneExpression within this matrix. I.E. "GeneScoreMatrix"
#' @param log2Norm If supplied geneSymbol, Log2 normalize the corresponding GeneScores/GeneExpression matrix before plotting.
#' @param upstream The number of basepairs upstream of the transcription start site of `geneSymbol` to extend the plotting window.
#' If `region` is supplied, this argument is ignored.
#' @param downstream The number of basepairs downstream of the transcription start site of `geneSymbol` to extend the plotting window.
#' If `region` is supplied, this argument is ignored.
#' @param tileSize The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
#' @param minCells The minimum number of cells contained within a cell group to allow for this cell group to be plotted. This argument can be
#' used to exclude pseudo-bulk replicates generated from low numbers of cells.
#' @param normMethod The name of the column in `cellColData` by which normalization should be performed. The recommended and default value
#' is "ReadsInTSS" which simultaneously normalizes tracks based on sequencing depth and sample data quality.
#' @param threads The number of threads to use for parallel execution.
#' @param ylim The numeric quantile y-axis limit to be used for for "bulkTrack" plotting. If not provided, the y-axis limit will be c(0, 0.999).
#' @param pal A custom palette (see `paletteDiscrete` or `ArchRPalettes`) used to override coloring for groups.
#' @param pal_loops A custom palette (see `paletteDiscrete` or `ArchRPalettes`) used to override coloring for loops.
#' @param scaleLim_loops A numeric vector of length 2 with lower and upper limit for co-accessibility score legend in loop track.
#' @param baseSize The numeric font size to be used in the plot. This applies to all plot labels.
#' @param scTileSize The width of the tiles in scTracks. Larger numbers may make cells overlap more. Default is 0.5 for about 100 cells.
#' @param scCellsMax The maximum number of cells for scTracks.
#' @param borderWidth The numeric line width to be used for plot borders.
#' @param tickWidth The numeric line width to be used for axis tick marks.
#' @param facetbaseSize The numeric font size to be used in the facets (gray boxes used to provide track labels) of the plot.
#' @param geneAnnotation The `geneAnnotation` object to be used for plotting the "geneTrack" object. See `createGeneAnnotation()` for more info.
#' @param title The title to add at the top of the plot next to the plot's genomic coordinates.
#' @param verbose A boolean value that determines whether standard output should be printed.
#' @param logFile The path to a file to be used for logging ArchR output.
#' @export

plotBrowserTrack <- function(
  ArchRProj = NULL, 
  region = NULL, 
  groupBy = "Clusters",
  useGroups = NULL, 
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 1.5, 3, 4),
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  geneSymbol = NULL,
  useMatrix = NULL,
  log2Norm = TRUE,
  upstream = 50000,
  downstream = 50000,
  tileSize = 250, 
  minCells = 25,
  normMethod = "ReadsInTSS",
  threads = getArchRThreads(), 
  ylim = NULL,
  pal = NULL,
  pal_loops = colorRampPalette(c("#f0f0f0", "black"))(100),
  scaleLim_loops = NULL,
  baseSize = 7,
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
  ){
  
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = "ArchRProj")
  ArchR:::.validInput(input = region, name = "region", valid = c("granges","null"))
  ArchR:::.validInput(input = groupBy, name = "groupBy", valid = "character")
  ArchR:::.validInput(input = useGroups, name = "useGroups", valid = c("character", "null"))
  ArchR:::.validInput(input = plotSummary, name = "plotSummary", valid = "character")
  ArchR:::.validInput(input = sizes, name = "sizes", valid = "numeric")
  ArchR:::.validInput(input = features, name = "features", valid = c("granges", "grangeslist", "null"))
  ArchR:::.validInput(input = loops, name = "loops", valid = c("granges", "grangeslist", "null"))
  ArchR:::.validInput(input = geneSymbol, name = "geneSymbol", valid = c("character", "null"))
  ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character", "null"))
  ArchR:::.validInput(input = log2Norm, name = "log2Norm", valid = c("boolean"))
  ArchR:::.validInput(input = upstream, name = "upstream", valid = c("integer"))
  ArchR:::.validInput(input = downstream, name = "downstream", valid = c("integer"))
  ArchR:::.validInput(input = tileSize, name = "tileSize", valid = c("integer"))
  ArchR:::.validInput(input = minCells, name = "minCells", valid = c("integer"))
  ArchR:::.validInput(input = normMethod, name = "normMethod", valid = c("character"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = ylim, name = "ylim", valid = c("numeric", "null"))
  ArchR:::.validInput(input = pal, name = "pal", valid = c("palette", "null"))
  ArchR:::.validInput(input = pal_loops, name = "pal_loops", valid = c("palette", "null"))
  ArchR:::.validInput(input = scaleLim_loops, name = "scaleLim_loops", valid = c("numeric", "null"))
  ArchR:::.validInput(input = baseSize, name = "baseSize", valid = "numeric")
  ArchR:::.validInput(input = scTileSize, name = "scTileSize", valid = "numeric")
  ArchR:::.validInput(input = scCellsMax, name = "scCellsMax", valid = "integer")
  ArchR:::.validInput(input = borderWidth, name = "borderWidth", valid = "numeric")
  ArchR:::.validInput(input = tickWidth, name = "tickWidth", valid = "numeric")
  ArchR:::.validInput(input = facetbaseSize, name = "facetbaseSize", valid = "numeric")
  geneAnnotation <- ArchR:::.validGeneAnnotation(geneAnnotation)
  ArchR:::.validInput(input = title, name = "title", valid = "character")

  tstart <- Sys.time()
  ArchR:::.startLogging(logFile=logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "plotBrowserTrack Input-Parameters", logFile = logFile)

  ##########################################################
  # Get Region Where Plot Will Occur (GenomicRanges)
  ##########################################################
  ArchR:::.logDiffTime("Validating Region", t1=tstart, verbose=verbose, logFile=logFile)
  if(is.null(region)){
    if(!is.null(geneSymbol)){
      region <- geneAnnotation$genes
      region <- region[which(tolower(mcols(region)$symbol) %in% tolower(geneSymbol))]
      region <- region[order(match(tolower(mcols(region)$symbol), tolower(geneSymbol)))]
      print(region)
      region <- resize(region, 1, "start")
      strand(region) <- "*"
      region <- extendGR(region, upstream = upstream, downstream = downstream)
    }
  }
  region <- ArchR:::.validGRanges(region)
  ArchR:::.logThis(region, "region", logFile = logFile)

  if(is.null(geneSymbol)){
    useMatrix <- NULL
  }

  if(!is.null(useMatrix)){
    featureMat <- ArchR:::.getMatrixValues(
      ArchRProj = ArchRProj,
      matrixName = useMatrix,
      name = mcols(region)$symbol
    )
    if(log2Norm){
      featureMat <- log2(featureMat + 1) 
    }
    featureMat <- data.frame(t(featureMat))
    featureMat$Group <- getCellColData(ArchRProj, groupBy, drop = FALSE)[rownames(featureMat), 1]
  }

  ggList <- lapply(seq_along(region), function(x){

    plotList <- list()

    ##########################################################
    # Bulk Tracks
    ##########################################################
    if("bulktrack" %in% tolower(plotSummary)){
      ArchR:::.logDiffTime(sprintf("Adding Bulk Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$bulktrack <- ArchR:::.bulkTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = minCells,
        pal = pal,
        ylim = ylim,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        facetbaseSize = facetbaseSize,
        normMethod = normMethod,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }
    
    ##########################################################
    # SC Tracks
    ##########################################################
    if("sctrack" %in% tolower(plotSummary)){
      ArchR:::.logDiffTime(sprintf("Adding SC Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$sctrack <- ArchR:::.scTracks(
        ArchRProj = ArchRProj, 
        region = region[x], 
        tileSize = tileSize, 
        groupBy = groupBy,
        threads = threads, 
        minCells = 5,
        maxCells = scCellsMax,
        pal = pal,
        baseSize = baseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth,
        scTileSize = scTileSize,
        facetbaseSize = facetbaseSize,
        geneAnnotation = geneAnnotation,
        title = title,
        useGroups = useGroups,
        tstart = tstart,
        logFile = logFile) + theme(plot.margin = unit(c(0.35, 0.75, 0.35, 0.75), "cm"))
    }

    ##########################################################
    # Feature Tracks
    ##########################################################
    if("featuretrack" %in% tolower(plotSummary)){
      if(!is.null(features)){
        ArchR:::.logDiffTime(sprintf("Adding Feature Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$featuretrack <- ArchR:::.featureTracks(
            features = features, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            title = "Peaks",
            logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }

    ##########################################################
    # Loop Tracks
    ##########################################################
    if("looptrack" %in% tolower(plotSummary)){
      if(!is.null(loops)){
        ArchR:::.logDiffTime(sprintf("Adding Loop Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
        plotList$looptrack <- .loopTracks(
            loops = loops, 
            region = region[x], 
            facetbaseSize = facetbaseSize,
            hideX = TRUE, 
            ### ADAPTED IL ###
            #hideY = TRUE,
            hideY = FALSE,
            pal = pal_loops,
            coaccess_lim = scaleLim_loops,
            ###
            title = "Loops",
            logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
      }
    }

    ##########################################################
    # Gene Tracks
    ##########################################################
    if("genetrack" %in% tolower(plotSummary)){
      ArchR:::.logDiffTime(sprintf("Adding Gene Tracks (%s of %s)",x,length(region)), t1=tstart, verbose=verbose, logFile=logFile)
      plotList$genetrack <- ArchR:::.geneTracks(
        geneAnnotation = geneAnnotation, 
        region = region[x], 
        facetbaseSize = facetbaseSize,
        title = "Genes",
        logFile = logFile) + theme(plot.margin = unit(c(0.1, 0.75, 0.1, 0.75), "cm"))
    }

    ##########################################################
    # Time to plot
    ##########################################################
    plotSummary <- tolower(plotSummary)
    names(sizes) <- plotSummary
    sizes <- sizes[order(plotSummary)]
    plotSummary <- plotSummary[order(plotSummary)]

    # nullSummary <- unlist(lapply(seq_along(plotSummary), function(x) is.null(eval(parse(text=paste0("plotList$", plotSummary[x]))))))
    # if(any(nullSummary)){
    #   sizes <- sizes[-which(nullSummary)]
    # }
    sizes <- sizes[tolower(names(plotList))]

    if(!is.null(useMatrix)){

      suppressWarnings(ArchR:::.combinedFeaturePlot(
        plotList = plotList,
        log2Norm = log2Norm,
        featureMat = featureMat,
        feature = region[x]$symbol[[1]],
        useMatrix = useMatrix,
        pal = pal,
        sizes = sizes,
        baseSize = baseSize,
        facetbaseSize = facetbaseSize,
        borderWidth = borderWidth,
        tickWidth = tickWidth
      ))

    }else{

      ArchR:::.logThis(names(plotList), sprintf("(%s of %s) names(plotList)",x,length(region)), logFile=logFile)
      ArchR:::.logThis(sizes, sprintf("(%s of %s) sizes",x,length(region)), logFile=logFile)
      #ArchR:::.logThis(nullSummary, sprintf("(%s of %s) nullSummary",x,length(region)), logFile=logFile)
      ArchR:::.logDiffTime("Plotting", t1=tstart, verbose=verbose, logFile=logFile)
      
      tryCatch({
        suppressWarnings(ggAlignPlots(plotList = plotList, sizes=sizes, draw = FALSE))
      }, error = function(e){
        ArchR:::.logMessage("Error with plotting, diagnosing each element", verbose = TRUE, logFile = logFile)
        for(i in seq_along(plotList)){
          tryCatch({
            print(plotList[[i]])
          }, error = function(f){
            ArchR:::.logError(f, fn = names(plotList)[i], info = "", errorList = NULL, logFile = logFile)
          })
        }
        ArchR:::.logError(e, fn = "ggAlignPlots", info = "", errorList = NULL, logFile = logFile)
      })

    }

  })

  if(!is.null(mcols(region)$symbol)){
    names(ggList) <- mcols(region)$symbol
  }else{
    if(length(ggList) == 1){
      ggList <- ggList[[1]]
    }
  }

  ArchR:::.endLogging(logFile=logFile)

  ggList

}
    



#######################################################
# Loop Tracks
#######################################################
.loopTracks <- function(
  loops = NULL, 
  region = NULL, 
  title = "LoopTrack", 
  pal = NULL,
  coaccess_lim = NULL,
  baseSize = 9, 
  facetbaseSize = 9,
  featureWidth = 2, 
  borderWidth = 0.4, 
  hideX = FALSE, 
  hideY = FALSE,
  logFile = NULL
  ){
  
  getArchDF <- function(lp, r = 1){
    angles <- seq(pi, 2*pi,length.out=100)
    rx <- (end(lp)-start(lp))/2
    rscale <- r * apply(data.frame(lp$PercAccess1, lp$PercAccess2), 1, mean)
    cx <- start(lp) + rx
    if(is.null(mcols(lp)$correlation)){
      mcols(lp)$correlation <- 1
    }
    df <- lapply(seq_along(cx), function(z){
      xz <- rx[z]*cos(angles)+cx[z]
      dfz <- DataFrame(x=xz, y=rscale[z]*sin(angles), id=Rle(paste0("l",z)), correlation = mcols(lp)$correlation[z])
    }) %>% Reduce("rbind",.)
    return(df)
  }

  if(!is.null(loops)){

    if(is(loops, "GRanges")){
      loops <- SimpleList(Loops = loops)
    }else if(ArchR:::.isGRList(loops)){
    }else{
      stop("Loops is not a GRanges or a list of GRanges! Please supply valid input!")
    }
        
    loopO <- lapply(seq_along(loops), function(x){
       subLoops <- subsetByOverlaps(loops[[x]], region, ignore.strand = TRUE, type = "within") 
       if(length(subLoops)>0){
         dfx <- getArchDF(subLoops)
         dfx$name <- Rle(paste0(names(loops)[x]))
         dfx
       }else{
         NULL
       }
    })
    loopO <- loopO[!vapply(loopO, is.null, logical(1))]
    loopO <- loopO %>% Reduce("rbind",.)
    ArchR:::.logThis(loopO, "loopO", logFile = logFile)
    
    if (is.null(coaccess_lim)){
        correlationMin <- min(loopO$correlation)
        correlationMax <- max(loopO$correlation)
    } else {
        correlationMin <- min(coaccess_lim)
        correlationMax <- max(coaccess_lim)
    }
    correlationMax_PercAccess <- min(loopO$y)
                                      
                                      
    testDim <- tryCatch({
      if(is.null(loopO)){
        FALSE
      }
      if(nrow(loopO) > 0){
        TRUE
      }else{
        FALSE
      }
    }, error = function(x){
      FALSE
    })

    if(testDim){

      loopO$facet <- title
      if(is.null(pal)){
        pal <- colorRampPalette(c("#E6E7E8","#3A97FF","#8816A7","black"))(100)
      }

      p <- ggplot(data = data.frame(loopO), aes(x = x, y = y, group = id, color = correlation)) +
        geom_line() +
        facet_grid(name ~ .) +
        ylab("Percent accessible cells/aggregates") + 
        coord_cartesian(ylim = c(correlationMax_PercAccess,0)) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        scale_color_gradientn(colors = pal, limits = c(correlationMin, correlationMax)) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth, legendPosition = "right") +
        theme(strip.text.y = element_text(size = facetbaseSize, angle = 0), strip.background = element_blank(),
          legend.box.background = element_rect(color = NA)) +
        guides(color= guide_colorbar(barwidth = 1, barheight = 4)) +
        theme(legend.text = element_text(size = 10), legend.title = element_text(size = facetbaseSize))
 
      ### From https://stackoverflow.com/questions/51784640/ggplot-how-to-retrieve-values-for-axis-labels (07/04/2022)
      if (utils::packageVersion("ggplot2") >= "3.3.0"){
          y_labs <- ggplot_build(p)$layout$panel_params[[1]]$y$get_labels()
          y_breaks <- ggplot_build(p)$layout$panel_params[[1]]$y$get_breaks()
      } else if (utils::packageVersion("ggplot2") >= "3.0.0"){
          y_labs <- ggplot_build(p)$layout$panel_params[[1]]$y.labels
          y_breaks <- ggplot_build(p)$layout$panel_params[[1]]$y.breaks
      } else {
          y_labs <- ggplot_build(p)$layout$panel_ranges[[1]]$y.labels
          y_breaks <- ggplot_build(p)$layout$panel_ranges[[1]]$y.breaks
      }
      
      p <- p + scale_y_continuous(name = "Percent accessible", breaks = y_breaks, labels = gsub("-", "", y_labs)) + 
                theme(axis.title.y=element_text(size=facetbaseSize), axis.text.y=element_text(size=10))
        
    }else{

      #create empty plot
      df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
      p <- ggplot(data = df, aes(start, end)) + 
        geom_point() +
        facet_grid(facet~.) +
        theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
        scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
        theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
        theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

    }

  }else{

    #create empty plot
    df <- data.frame(facet = "LoopTrack", start = 0, end = 0, strand = "*", symbol = "none")
    p <- ggplot(data = df, aes(start, end)) + 
      geom_point() +
      facet_grid(facet~.) +
      theme_ArchR(baseSize = baseSize, baseLineSize = borderWidth, baseRectSize = borderWidth) +
      scale_x_continuous(limits = c(start(region), end(region)), expand = c(0,0)) +
      theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
      theme(axis.title.y=element_blank(), axis.text.y=element_blank(),axis.ticks.y=element_blank())

  }

  if(hideX){
    p <- p + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  }

  if(hideY){
    p <- p + theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
  }

  if(!is.ggplot(p)){
    ArchR:::.logError("loopTracks is not a ggplot!", fn = "ArchR:::.loopTracks", info = "", errorList = NULL, logFile = logFile)
  }

  return(p)

}                                                            
                                  
### Copied from https://github.com/GreenleafLab/ArchR/blob/968e4421ce7187a8ac7ea1cf6077412126876d5f/R/ArchRBrowser.R on 05/04/2022
