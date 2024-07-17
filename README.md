# RWireX
RWireX is a R software package that provides various functions to analyze co-accessibility in scATAC-seq data. It is implemented as an extension to the ArchR software package (Granja et al. 2021) building on its existing funcitonalities. It computes Pearson correlation coefficients across different cell populations and at varying levels of resolution to identify both autonomous links of co-accessibility (ACs) and domains of contiguous co-accessibility (DCs). <br />

%% scheme from Fig. 3A

The “single-cell co-accessibility” workflow identifies ACs from stochastic accessibility changes in ATAC peaks with a homogeneous population of single cells as input. Pearson correlation coefficients between two peaks are assessed against a local background model. Background co-accessibility is determined from the 99th percentile of co-accessibility from accessibility matrices per chromosome shuffled over cells and peaks. The stability of ACs is assessed from their prevalence in the single cell population by computing the average percent accessible cells of the linked peaks.  <br />
<br />
The “metacell co-accessibility” workflow identifies DCs from perturbed accessibility changes in 10 kb genomic tiles and aggregated metacell profiles of cells with similar chromatin accessibility profiles. It requires cell populations that are heterogeneous in respect to the perturbation. This allows to identify broader genomic patterns of depleted or enriched co-accessibility along the genomic coordinate.

## Installation

Install required packages: <br />
BiocManager::install("plotgardener") # only needed for plotCoAccessibilityMap functionality <br />
devtools::install_github("GreenleafLab/ArchR", ref="dev", repos = BiocManager::repositories()) <br />
devtools::install_github("dozmorovlab/SpectralTAD") <br />
<br />
Install RWireX: <br />
devtools::install_github("RippeLab/RWireX") <br />

## Vignettes

see Wiki <br />

## How to cite RWireX

%% doi of biorxiv

## References
Granja, J.M. et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nat Genet 53, 403-411 (2021).

