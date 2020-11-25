# TCA: Testis Cell Atlas

## What is TCA?
Testis Cell Atlas (TCA) is a scRNA-seq database focusing on spermatogenesis. TCA provides detailed cell-type annotation at the single-cell level, enabling the exploration of testis across different species.

## Link
http://tca.xielab.tech/

## Requirements
NOTE: Works with R > 3.6.0 or greater.
Please install the required packages.
```
pkgs<-c("Seurat", "dplyr", "R.utils", "data.table", "monocle", "GSVA")
Install.packages(pkgs)
```


## What is TCA?
TCA provides 5 major modules which allow users to explore cellular diversity of spermatogenesis, including gene expression visualization, cell type automatically annotation, pathway enrichment, literature mining, and gene correlation analysis. TCA is an all-in-one resource of spermatogenesis and will help reproductive biologists to test their hypotheses.

## Usage
To repeat our analysis, you can clone this repository and run

```
git clone github.com/XGCLab/TCA
Rscript 01_pretreat_data.R
```


## Contact
Prof. Gangcai Xie, PhD

gangcai (at) ntu.edu.cn

If you have technical problems, please contact Yingcheng Wu (wuyc at mail.com).

## Citation
NA


