# RWireX
RWireX is a R software package and contains various functions to analyze co-accessibility in scATAC-seq data acquired by plate- or droplet-based single cell sequencing. It is implemented as an extension to the ArchR software package building on its existing funcitonalities. 

RWireX provides 
- separated workflows to compute cis and trans co-accessibility, 
- employs a local background model for co-accessibility link evaluation, 
- computes single-cell link activity scores that enable correlation to gene expression output, 
- and provides co-accessibility maps of large genomic regions to identify interacting loci without pre-selecting specific viewpoints.

## Installation

Install required packages:
BiocManager::install("plotgardener") # only needed for plotCoAccessibilityMap functionality
devtools::install_github("GreenleafLab/ArchR", ref="dev", repos = BiocManager::repositories())

Create and save your personal token for installation at https://github.com/settings/tokens 
personal_token = ""
devtools::install_github("https://github.com/RippeLab/RWireX", auth_token = personal_token)


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
