# RWireX
RWireX is an extension to the ArchR software package. 
It performs co-accessibility analysis of scATAC-seq data acquired by droplet-based single cell sequencing. 

RWireX employs a local background model, computes single-cell link activity scores, and provides co-accessibility maps of large genomic regions to identify interacting loci without pre-selecting specific viewpoints.

## Functions

All methods return the co-accessibility object and not ArchR Object, as it is in the ArchR package.
### GetCoAccessibility
An extended version of the addCoAccessibility by ArchR. Allows unique aggregation and accessibility analyis without aggregation - on the single cell level.
Additionally provides scaled/normalizred matrix data and feature set used to calculate co-accessibility ad metadata. 

### GetCoAccessibilityChromosomeWise
Similar to GetCoAccessibility, but the co-accessibility is calculated for all peak-pairs in the chromosome and not for peak overlaps within maxDist, as done in the original method.

### GetBackgroundCoAccessibility
This function calculates background co-accessibility scores for peaks of a given ArchRProject. The background co-accessibility is determined by shuffling of features and cells.

## In development: 
1. Background Chromosomewise
2. Peaks co-accessibility between chromosomes
