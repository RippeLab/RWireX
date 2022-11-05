#' plotCoAccessibilityMap
#' 
#' This function returns co-accessibility map of a specified genomic region.
#'
#' @param ArchRProj An `ArchRProject` object.
#' @param CoAccessibility A GRanges object with co-accessible links from getCoAccessibility function.
#' @param region GRanges object of length 1 with region of interest for which co-accessibility should be visualized.
#' @param genome Genome and gene annotation to use.
#' @param useMatrix The name of the data matrix used for calculation of CoAccessibility. Options include "PeakMatrix" or non-binarized "TileMatrix" (binarize = FALSE in addTileMatrix).
#' @param main Title of co-accessibility map.
#' @param onlyPos Only visualize positive co-accessibility scores.
#' @param scaleLim Set limits of co-accessibility scale.
#' @param rescale If TRUE, genomic coordinate is rescaled to visualize peaks continously. 
#' @return Co-accessibility map
#' @keywords co-accessibility map coAccessibilityMap archr rwire
#' @examples Coming soon.
#' @export


plotCoAccessibilityMap <- function(
    ArchRProj, 
    CoAccessibility, 
    region, 
    genome = getArchRGenome(), 
    useMatrix = "PeakMatrix",
    main = NULL, 
    onlyPos = FALSE,
    scaleLim = NULL,
    rescale = FALSE){
    
    ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
    ArchR:::.validInput(input = CoAccessibility, name = "CoAccessibility", valid = c("GRanges"))
    ArchR:::.validInput(input = region, name = "region", valid = c("GRanges"))
    ArchR:::.validInput(input = genome, name = "genome", valid = c("character"))
    ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
    ArchR:::.validInput(input = main, name = "main", valid = c("character", "null"))
    ArchR:::.validInput(input = onlyPos, name = "onlyPos", valid = c("logical"))
    ArchR:::.validInput(input = scaleLim, name = "scaleLim", valid = c("numeric", "null"))
    ArchR:::.validInput(input = rescale, name = "rescale", valid = c("logical"))
    
    ### Get co-accessible links in region of interest
    dat <- CoAccessibility[findOverlaps(CoAccessibility, region, type = "within") %>% queryHits(.)] %>% as.data.frame(.)
    dat <- dat[, c("start", "end", "correlation")]

    ### Get Peak Set
    if (useMatrix == "PeakMatrix"){
        peaks <- getPeakSet(ArchRProj)
        seqlevels(peaks, pruning.mode="tidy") <- seqlevels(region)
    } else if (useMatrix == "TileMatrix"){
        tileSet <- getMatrixFromProject(ArchRProj, useMatrix = "TileMatrix")@elementMetadata
        peaks <- GRanges(seqnames = tileSet$seqnames, 
                           ranges = IRanges(start = tileSet$start, width = tileSet$start[2]-tileSet$start[1]))
        peaks$idx <- tileSet$idx
        peaks$id <- 1:length(peaks)
        seqlevels(peaks, pruning.mode="tidy") <- seqlevels(region)
        if (rescale){
          print("NOT allowed to RESCALE tile matrix. Setting RESCALE to FALSE!")
          rescale = FALSE
        }
    }
    else{
      stop("This matrix is not supported yet. Please use Peak or Tile Matrix instead.")
    }
    
    resolution <- peaks@ranges@width %>% unique(.)
    if (length(resolution) != 1){
      stop("Peaks must be the same width. Consider using another algorithm or resize your peaks.")
    }
    ### Get genome
    genome <- tolower(genome)
    
    ### If requested, rescale peaks to continous scale
    if (rescale){
        origTiles <- peaks[findOverlaps(region, peaks) %>% subjectHits(.)] %>% ranges(.) %>% as.data.frame(.)
        colnames(origTiles) <- c("PeakStart", "PeakEnd", "PeakWidth", "PeakOrigin")
        origTiles$PeakID <- rowMeans(origTiles[,c("PeakStart", "PeakEnd")])
        origTiles <- origTiles[order(origTiles$PeakID, decreasing = FALSE),]
    
        numTiles <- nrow(origTiles)
    
        newTiles <- tile(region, numTiles)[[1]] %>% as.data.frame(.)
        newTiles$ID <- rowMeans(newTiles[,c("start", "end")])
        newTiles$num <- 1:nrow(newTiles)
        newTiles <- cbind(newTiles, origTiles)
        
        dat$start <- newTiles$start[match(dat$start, newTiles$PeakID)]
        dat$end <- newTiles$start[match(dat$end, newTiles$PeakID)]
        
        resolution <- newTiles$width[1] 
    } else {
        ### Shift start and end points of co-accessible links by 1/2 peak size to the left (cannot use mid-points as link start)
        dat$start <- dat$start - resolution/2
        dat$end <- dat$end - resolution/2
    }
  
    if (is.null(main)){
        main <- as.character(region)
    }
    
    ## Create a plotgardener page
    if (useMatrix == "PeakMatrix" & !rescale){
        height <- 15
    } else if (useMatrix == "TileMatrix"){
        height <- 13
    } else if (useMatrix == "PeakMatrix" & rescale){
        height <- 10
    }
    pageCreate(width = 21, height = height, default.units = "inches", showGuides = FALSE)

    ## Set genomic and dimension parameters in a `params` object
    params_obj <- pgParams(chrom = region@seqnames@values, chromstart = region@ranges@start, chromend = region@ranges@start + region@ranges@width - 1, 
                          assembly = genome,
                          x = 0.5, width = 20, default.units = "inches")
    ## Main title
    plotText(label = main, fontcolor = "black", fontsize = 22,
                x = 1, y = 1, just = c("left","top"), default.units = "inches")

    ## Plot Hi-C triangle
    correlation_range <- c(-max(max(dat$correlation), abs(min(dat$correlation))), max(max(dat$correlation), abs(min(dat$correlation))))
    if (!is.null(scaleLim)){
        correlation_range <- scaleLim
    }
    
    pal <- colorRampPalette(colors = c("blue", "white", "red"))
    if (onlyPos){
        pal <- colorRampPalette(colors = c("white", "white", "red"))
    }
   
    hic_gm <- plotHicTriangle(data = dat, params = params_obj, bg ="white", palette = pal,
                                 zrange = correlation_range, resolution = resolution, x = 0.5, width = 20,
                                 y = 10, height = 10, just = c("left", "bottom"))

    ## Annotate Hi-C heatmap legend
    annoHeatmapLegend(plot = hic_gm, fontsize = 15, fontcolor = "black", digits = 2,
                      x = 17, y = 3, width = 3, height = 0.5, orientation = "h",
                      just = c("centre", "top"), default.units = "inches")
    plotText(label = "Co-accessibility score", fontcolor = "black", fontsize = 18,
                x = 17, y = 2.75, just = c("centre","center"), default.units = "inches")


    ## Plot genes
    if (useMatrix == "PeakMatrix" & !rescale){
        y_genes <- 12
        y_label <- 13.5
    } else if (useMatrix == "TileMatrix"){
        y_genes <- 10
        y_label <- 11.5
    }
    
    if (!rescale){
      plotGenes(params = params_obj, stroke = 1, fontsize = 15, fontcolor = c("red", "dodgerblue2"), fill = c("red", "dodgerblue2"), bg = "white",
                strandLabels = FALSE, y = y_genes, height = 3)
      plotText(label = "Genes", fontcolor = "black", fontsize = 18, rot = 90,
               x = 0.1, y = y_label, just = c("top","center"), default.units = "inches")
    }
   
    ## Plot peaks
    if (useMatrix == "PeakMatrix" & rescale){
        plotText(label = "Rescaled peaks", fontcolor = "black", fontsize = 18,
                x = 1, y = 2, just = c("left","top"), default.units = "inches")
    } else if (useMatrix == "PeakMatrix"  & !rescale){
        plotRanges(data = peaks, params = params_obj, y = 10, height = 2, bg = "white", fill = "grey")
        plotText(label = "Peaks", fontcolor = "black", fontsize = 18, rot = 90,
                x = 0.1, y = 11.5, just = c("top","left"), default.units = "inches")
    }  

    ## Annotate genome label
    plotGenomeLabel(params = params_obj, 
                           scale = "Mb", fontsize = 15, y = 10, length = 20)
}