# RWireX
RWireX is a R software package that provides various functions to analyze co-accessibility in scATAC-seq data. It is implemented as an extension to the ArchR software package (Granja et al. 2021) building on its existing funcitonalities. It computes Pearson correlation coefficients across different cell populations and at varying levels of resolution to identify both autonomous links of co-accessibility (ACs) and domains of contiguous co-accessibility (DCs). <br />

![RWireX_scheme](/figures/RWireX_scheme.png)

The “single-cell co-accessibility” workflow identifies ACs from stochastic accessibility changes in ATAC peaks with a homogeneous population of single cells as input. Pearson correlation coefficients between two peaks are assessed against a local background model. Background co-accessibility is determined from the 99th percentile of co-accessibility from accessibility matrices per chromosome shuffled over cells and peaks. The stability of ACs is assessed from their prevalence in the single cell population by computing the average percent accessible cells of the linked peaks.  <br />
<br />
The “metacell co-accessibility” workflow identifies DCs from perturbed accessibility changes in genomic tiles and aggregated metacell profiles of cells with similar chromatin accessibility profiles. It requires cell populations that are heterogeneous in respect to a perturbation or cell type/state. This allows to identify broader genomic patterns of depleted or enriched co-accessibility along the genomic coordinate.

## Installation

Install required packages: <br />
```ruby
BiocManager::install("plotgardener")
```
```ruby
devtools::install_github("GreenleafLab/ArchR", repos = BiocManager::repositories())
```
```ruby
devtools::install_github("dozmorovlab/SpectralTAD")
```
<br />
Install RWireX: <br />

```ruby
devtools::install_github("RippeLab/RWireX")
```

## Vignettes
We provide test data and have vignettes for single-cell as well as metacell co-accessibility. Check out our [wiki](https://github.com/RippeLab/RWireX/wiki).

## How to cite RWireX

%% doi of biorxiv

## References
Granja, J.M. et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nat Genet 53, 403-411 (2021).

