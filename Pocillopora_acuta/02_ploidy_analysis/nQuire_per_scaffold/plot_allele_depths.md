---
title: "Plot allele depths for each *P. acuta* scaffold"
author: "Timothy Stephens"
date: "24/09/2022"
output: 
  html_document:
    keep_md: yes
---



# Setup

Setup R env. Load packages and set default image export formats, size and resolution.


```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.height = 12, 
                      fig.width = 8, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(tibble)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
options(scipen = 999) #Prevent scientific notation
```





# List files to plot

List files with allelic proportions that we want to plot. Also get titles that we want to use for the plots.

```r
awk 'NR>1' Pocillopora_acuta_KBHIv2.assembly.fasta.1Mbp_scaffold_features.bed.lrdmodel.nQuire_results.tsv \
  | sort -k19,19 -k1,1 \
  | awk '{split($1, a, "-"); print $1".denoise.bin.coverage.sitesProp.gz\t"a[2]" ("$19")"}' \
  | sed -e 's/.coordsorted-/.coordsorted.lowerProp-/' \
  > samples_list.txt 
```





# Functions

Function for plotting AP


```r
plot_AP = function(x.lim.min, x.lim.max)
{
  files2plot <- read.table("samples_list.txt", header = F, sep='\t')
  nsamples <- nrow(files2plot)
  pdf(paste("plot_xLim_", x.lim.min, "-", x.lim.max, "_Pacuta_scaffolds_AP_combined.pdf", sep=''))
  par(mar=c(5, 3, 3, 2), cex=1.5, mfrow=c(4,2)) # number of plots per page
  for (i in 1:nsamples) {
    AP <- read.table(files2plot[i,1], header = F)
    d <- density(AP[,1], from=x.lim.min, to=x.lim.max, bw=0.01, na.rm =T)
    plot(d, xlim = c(x.lim.min,x.lim.max), main=as.character(files2plot[i,2]), col="blue", xlab = (dim(AP)[1])/2, lwd=2)
  }
  dev.off()
}
```





# Load datasets into R

Load allele depth files and plot distributions.


```r
plot_AP(0.0, 1.0)
```

```
## pdf 
##   2
```




```r
plot_AP(0.1, 0.9)
```

```
## pdf 
##   2
```




```r
plot_AP(0.2, 0.8)
```

```
## pdf 
##   2
```





# Session Info


```r
sessionInfo()
```

```
## R version 4.1.2 (2021-11-01)
## Platform: x86_64-apple-darwin18.7.0 (64-bit)
## Running under: macOS Mojave 10.14.6
## 
## Matrix products: default
## BLAS:   /usr/local/Cellar/openblas/0.3.18/lib/libopenblasp-r0.3.18.dylib
## LAPACK: /usr/local/Cellar/r/4.1.2/lib/R/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] RColorBrewer_1.1-3 ggplot2_3.3.6      reshape2_1.4.4     tibble_3.1.8      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.9       pillar_1.8.1     bslib_0.4.0      compiler_4.1.2  
##  [5] jquerylib_0.1.4  plyr_1.8.7       tools_4.1.2      digest_0.6.29   
##  [9] jsonlite_1.8.0   evaluate_0.16    lifecycle_1.0.2  gtable_0.3.1    
## [13] pkgconfig_2.0.3  rlang_1.0.5      DBI_1.1.3        cli_3.4.0       
## [17] rstudioapi_0.14  yaml_2.3.5       xfun_0.33        fastmap_1.1.0   
## [21] withr_2.5.0      stringr_1.4.1    dplyr_1.0.10     knitr_1.40      
## [25] generics_0.1.3   vctrs_0.4.1      sass_0.4.2       tidyselect_1.1.2
## [29] grid_4.1.2       glue_1.6.2       R6_2.5.1         fansi_1.0.3     
## [33] rmarkdown_2.16   purrr_0.3.4      magrittr_2.0.3   scales_1.2.1    
## [37] htmltools_0.5.3  assertthat_0.2.1 colorspace_2.0-3 utf8_1.2.2      
## [41] stringi_1.7.8    munsell_0.5.0    cachem_1.0.6
```
