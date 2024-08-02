# RWireX: Co-accessibility analysis for scATAC-seq data
RWireX is a R software package that provides various functions to analyze co-accessibility in scATAC-seq data. It is implemented as an extension to the ArchR software package (Granja et al. 2021) building on its existing funcitonalities. RWireX offers two distinct workflows to uncover different aspects of chromatin organization: <br />
 - Single cell co-accessibility workflow: <br />
   Identifies autonomous links of co-accessibility (ACs) by analyzing stochastic accessibility changes in individual cells.
 - Metacell co-accessibility workflow: <br />
   Detects domains of contiguous co-accessibility (DCs) using aggregated profiles (metacells) from cells with similar chromatin accessibility.
<br />

![RWireX_scheme](/figures/RWireX_scheme_short.png)

These workflows enable researchers to identify both long-range interactions (ACs) and local domains of contiguous co-accessibility (DCs), providing a comprehensive view of chromatin organization and its relationship to gene regulation.
Key features:
- Compatible with ArchR-processed scATAC-seq data
- Flexible analysis options for different experimental designs
- Visualization tools for co-accessibility matrices
- Statistical methods to assess the significance of identified features
<br />
RWireX has been validated on data from various cellular systems, including human endothelial cells, mouse embryonic stem cells, and leukemia models. It is particularly effective with deep coverage scATAC-seq data (≥20,000 unique fragments/cell). A test data set of 200 cells is available at https://www.doi.org/10.5281/zenodo.13142236. For further details on RWireX see our bioRxiv preprint by Seufert et al.
By uncovering ACs and DCs, RWireX provides valuable insights into the chromatin-level mechanisms coordinating gene expression programs, applicable to a wide range of biological questions and model systems. <br />
<br />

## Installation

Install required packages: <br />

```r
install.packages("devtools")
install.packages("BiocManager")
```
```r
BiocManager::install("plotgardener")
```
```r
devtools::install_github("GreenleafLab/ArchR", repos = BiocManager::repositories())
ArchR::installExtraPackages()
```
```r
devtools::install_github("dozmorovlab/SpectralTAD")
```
<br />
Install RWireX: <br />

```r
devtools::install_github("RippeLab/RWireX")
```
<br />

## Download test data
We provide test data from untreated and 30 min/240 min TNFa-treated HUVECs to run the vignette. The test data is approximately 0.6 GB in size.

```r
library(RWireX)
data_paths <- downloadVignetteData()
```
<br />

## How to run RWireX
We provide vignettes as well as extended explanations on how to run single cell as well as metacell co-accessibility workflows in our [wiki](https://github.com/RippeLab/RWireX/wiki).
<br /><br />

## How to cite RWireX
Seufert, I. et al.
<br /><br />

## References
Granja, J.M. et al. ArchR is a scalable software package for integrative single-cell chromatin accessibility analysis. Nat Genet 53, 403-411 (2021).

