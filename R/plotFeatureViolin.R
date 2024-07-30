#' Create a violin plot of cell values per group, e.g. samples.
#' 
#' This function returns a plot.
#' 
#' @param ArchRProj An `ArchRProject` object.
#' @param feature Peak to plot in the format chr:idx, for example "chr1:15" if 
#' you use PeakMatrix, gene name if you use GeneScoreMatrix, or one of the data columns names if you use cellColData.
#' @param useMatrix Type of Matrix to plot data from
#' @param binarize If you plot values from a binary matrix. MUST be specified if the matrix is binarized. 
#' @param matrix You can provide matrix directly for a faster run.
#' @param groupBy Which groups should be displayed in the plot.
#' @param normBy Specify if you want the values in the matrix normalized.
#' @param orderedSamples You can specify the order of the group in the plot here in a string vector.
#' @param useGroups Filter the groups.
#' @param pseudolog If you would like the values to be plotted in the pseudolog.
#' @param pal Color pallete
#' @export
plotFeatureViolin <- function(ArchRProj,
                              feature, # Gene name or peak id (chr:idx)
                              useMatrix = "PeakMatrix",
                              binarize = FALSE,
                              matrix=NULL,
                              groupBy = "Sample",
                              normBy="nFrags", # NULL for no normalization
                              orderedSamples=NULL, 
                              useGroups=NULL, 
                              pseudolog = FALSE,
                              pal = NULL
){
  if (binarize){
    message("Warning! While it is possible to create a violin plot for binarized data, please consider a bar plot.")
    if (!is.null(normBy)){
      stop("You can not normalize binary values of the matrix. Please set normBy=NULL")
    }
  }
  
  if (is.null(matrix)){
    if (useMatrix == "cellColData"){
      matrix <- ArchRProj@cellColData
    } else {
      matrix = ArchR:::getMatrixFromProject(ArchRProj, useMatrix = useMatrix, binarize = binarize)
    }
  }
  if (useMatrix == "PeakMatrix"){
    value_idx = which(paste(matrix@rowRanges@seqnames, matrix@rowRanges$idx, sep=":") == feature)
    values = matrix@assays@data[[1]][value_idx,]
  } else if (useMatrix == "GeneScoreMatrix"){
    value_idx = which(matrix@elementMetadata$name==feature)
    values = matrix@assays@data[[1]][value_idx,]
  } else if (useMatrix == "cellColData"){
    value_idx = which(colnames(matrix)==feature)
    values = matrix[, value_idx]
  }
  
  
  
  if (is.null(normBy)){
    final_values <- as.matrix(values)
  } else {
    norm_factors <- ArchR:::getCellColData(ArchRProj, normBy, drop=FALSE)
    norm_factors[,1] <- median(norm_factors[,1]) / norm_factors[,1]
    
    mColSums = Matrix::colSums(matrix@assays@data[[1]])
    scaleTo = 10^4
    norm_factors[,1] <- norm_factors[,1] * (scaleTo / median(norm_factors[names(mColSums), 1] * mColSums))
    final_values =  values * as.matrix(norm_factors)
  }

  samples = getSampleNames(ArchRProj)
  
  if (useMatrix == "cellColData"){
    df = data.frame(Group=ArchR:::getCellColData(ArchRProj, groupBy, drop=FALSE)[,1],
                    Value= final_values[,1])
    df$Group <- as.factor(df$Group)
    df$Value <- as.numeric(df$Value)
  }else{
    #this step is required because the names of the cells are not sotred 
    #in the same way as the are in metadata of the ArchRProject
    names(final_values) = names(values)
    df = data.frame(Group=str_extract(names(final_values), paste0(samples, collapse='|')),
                    Value= final_values[,1])
  }
  
  if (!is.null(useGroups)){
    df <- subset(df, subset = Group %in% useGroups)
  }
  
  if (useMatrix == "PeakMatrix"){
    range = matrix@rowRanges[paste(matrix@rowRanges@seqnames, matrix@rowRanges$idx, sep=":") == feature]
    main <- paste0("Peak coordinates: ",range@seqnames, ":", range@ranges@start, "-", as.character(range@ranges@start+range@ranges@width-1))
    if (is.null(normBy)){
      ylabel <- "Peak accessibility"
    } else {
      ylabel <- paste0("Peak accessibility (norm. by: ", normBy, ")")
    }
  } else if (useMatrix == "GeneScoreMatrix"){
    main <- feature
    if (is.null(normBy)){
      ylabel <- "Gene score"
    } else {
      ylabel <- paste0("Gene score (norm. by: ", normBy, ")")
    }
  } else if (useMatrix == "cellColData"){
    main <- feature
    if (is.null(normBy)){
      ylabel <- "Score"
    } else {
      ylabel <- paste0("Score (norm. by: ", normBy, ")")
    }
  }
  
  p <- ggplot(df, aes(x=Group, y=Value)) + geom_violin(alpha = 0.5, width = 1, aes(fill = Group)) + 
    geom_boxplot(width = 0.1) + 
    theme_minimal() +
    theme(legend.position = "none", plot.title = element_text(size=16), axis.title = element_text(size=16), axis.text = element_text(size=16), axis.text.x=element_text(angle=90, vjust=1, hjust=1), legend.title = element_text(size=18), legend.text = element_text(size=16)) + 
    xlab(groupBy) + 
    ylab(ylabel) +
    ggtitle(main)
  
  if (pseudolog){
    p = p + scale_y_continuous(trans="pseudo_log")
  }
  if (!is.null(orderedSamples)){
    p = p + scale_x_discrete(limits = orderedSamples)
  }
  if (!is.null(pal)){
    p <- p + scale_fill_manual(values = pal)
  }
  return(p)
}