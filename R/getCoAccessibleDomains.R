#' getCoAccessibleDomains
#' 
#' This function returns a GRanges object of co-accessibly enriched domains.
#'
#' @param CoAccessibility A GRanges object with co-accessible links from getCoAccessibility function.
#' @param enrich_cutoff Quantile cutoff to determine co-accessibility enrichment.
#' @return GRanges object
#' @keywords co-accessibility domain archr rwirex
#' @examples Coming soon.
#' @export


getCoAccessibleDomains <- function(
    CoAccessibility, 
    enrich_cutoff = 0.9, 
    logFile = createLogFile("getCoAccessibleDomains"), 
    verbose = TRUE){
    
    ArchR:::.validInput(input = CoAccessibility, name = "CoAccessibility", valid = c("list"))
    ArchR:::.validInput(input = enrich_cutoff, name = "genome", valid = c("numeric"))
    ArchR:::.validInput(input = verbose, name = "verbose", valid = c("logical"))
    
    if (class(CoAccessibility) == "DFrame"){
      CoAccessibility = .createLoopsSameChromosome(CoAccessibility)
    }
    
    tstart <- Sys.time()
    ArchR:::.startLogging(logFile = logFile)
    ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), "getCoAccessibleDomains Input-Parameters", logFile = logFile)

    ### Generate HiC-like interaction data frames
    ArchR:::.logDiffTime(main = "Generating HiC-like interaction data frames", t1 = tstart, verbose = verbose, logFile = logFile)
    chrs <- CoAccessibility$CoAccessibility@seqnames@values %>% unique(.)
    hic_list <- lapply(chrs, function(chr){hic_obj <- CoAccessibility$CoAccessibility[CoAccessibility$CoAccessibility@seqnames == chr] %>% as.data.frame(.);
                                           hic_obj <- hic_obj[, c("start", "end", "correlation")];
                                           hic_obj$correlation[hic_obj$correlation < 0] <- 0;
                                           return(hic_obj)})
    
    ### Call domains with SpectralTAD
    ArchR:::.logDiffTime(main = "Calling domains with SpectralTAD", t1 = tstart, verbose = verbose, logFile = logFile)
    features <- CoAccessibility@metadata$featureSet
    resolution <- features@ranges@width %>% unique(.)
    domains_small <- try(SpectralTAD::SpectralTAD_Par(cont_list = hic_list, grange = TRUE,
                                                  chr = chrs, resolution = resolution,
                                                  levels = 3, min_size = 2, window_size = 20,
                                                  qual_filter = FALSE, z_clust = FALSE), silent = TRUE)
    domains_large <- try(SpectralTAD::SpectralTAD_Par(cont_list = hic_list, grange = TRUE,
                                                  chr = chrs, resolution = resolution,
                                                  levels = 3, min_size = 20, window_size = 200,
                                                  qual_filter = FALSE, z_clust = FALSE), silent = TRUE)

    ### Merge hierachical levels and remove duplicates
    domains_small <- lapply(domains_small, function(x){x <- unlist(x);
                                                       x <- x[!duplicated(x)];
                                                       return(x)}) %>% do.call("c", .)
    domains_large <- lapply(domains_large, function(x){x <- unlist(x);
                                                       x <- x[!duplicated(x)];
                                                       return(x)}) %>% do.call("c", .)

    ### Extend domains from border midpoint
    domains_small <- resize(domains_small, fix = 'start', width = width(domains_small) + 5000) 
    domains_small <- resize(domains_small, fix = "end", width = width(domains_small) + 4999)
    domains_large <- resize(domains_large, fix = 'start', width = width(domains_large) + 5000) 
    domains_large <- resize(domains_large, fix = "end", width = width(domains_large) + 4999)
  
    ### Calculate domain co-accessibility scores
    ArchR:::.logDiffTime(main = "Calculating domain co-accessibility scores", t1 = tstart, verbose = verbose, logFile = logFile)
    domains_small$coaccessibility_score <- 0
    for (i in 1:length(domains_small)){
        lookup <- findOverlaps(CoAccessibility$CoAccessibility, domains_small[i], type = "within")
        cor_vals <- CoAccessibility$CoAccessibility$correlation[queryHits(lookup)]
        domains_small$coaccessibility_score[i] <- mean(cor_vals)
    }
    domains_large$coaccessibility_score <- 0
    for (i in 1:length(domains_large)){
        lookup <- findOverlaps(CoAccessibility$CoAccessibility, domains_large[i], type = "within")
        cor_vals <- CoAccessibility$CoAccessibility$correlation[queryHits(lookup)]
        domains_large$coaccessibility_score[i] <- mean(cor_vals)
    }
  
    ### Filter co-accessibility domains
    domains <- c(domains_small[domains_small$coaccessibility_score >= quantile(domains_small$coaccessibility_score, enrich_cutoff)],
                 domains_large[domains_large$coaccessibility_score >= quantile(domains_large$coaccessibility_score, enrich_cutoff)])
    ArchR:::.endLogging(logFile = logFile)
    return(domains)
}
