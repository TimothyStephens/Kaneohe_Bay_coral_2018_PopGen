---
title: "Plot nQuire results"
author: "Timothy Stephens"
date: "13/09/2022"
output: 
  html_document:
    keep_md: yes
---

Takes nQuire delta log-liklihood values for diploid/triploid/tetraploid models and plots values in a single dot plot.

# Setup

Setup R env. Load packages and set default image export formats, size and resolution.


```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.width = 8, 
                      fig.height = 12, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(tibble)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(patchwork)
library(cowplot)
```

```
## 
## Attaching package: 'cowplot'
```

```
## The following object is masked from 'package:patchwork':
## 
##     align_plots
```

```r
options(scipen = 999) #Prevent scientific notation
```




# Load results and plot samples

## *M. capitata* RNA-Seq samples from this study


```r
xlab <- "M. capitata RNA-seq samples"

file.name.1 <- "../../../Montipora_capitata/02_ploidy_analysis/nQuire/nQuire_results.tsv"
data.1 <- read.table(file.name.1, header = T, sep = '\t')
file.name.2 <- "../../../samples_from_SRA/Montipora/02_ploidy_analysis/nQuire/nQuire_results.tsv"
data.2 <- read.table(file.name.2, header = T, sep = '\t')

data <- rbind(data.1, data.2)

data.melted <- data %>% 
  select(c("ID", "Ploidy", "Diploid", "Triploid", "Tetraploid")) %>% 
  melt(id.vars = c("ID", "Ploidy")) %>%
  mutate(ID = factor(ID, data$ID)) %>%
  mutate(variable = factor(variable, c("Diploid", "Triploid", "Tetraploid")))

data.denoised.melted <- data %>% 
  select(c("ID", "Ploidy", "Diploid_denoised", "Triploid_denoised", "Tetraploid_denoised")) %>% 
  melt(id.vars = c("ID", "Ploidy")) %>%
  mutate(ID = factor(ID, data$ID)) %>%
  mutate(variable = factor(variable, c("Diploid_denoised", "Triploid_denoised", "Tetraploid_denoised")))

# Plot raw sites
MC.raw <- ggplot(data.melted, aes(x=ID, y=value, group=variable)) +
  geom_point(aes(color=variable)) +
  scale_color_brewer(palette="Dark2") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab(xlab) +
  ylab("delta Log-Likelihood") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 3))

# Plot denoised sites
MC.denoised <- ggplot(data.denoised.melted, aes(x=ID, y=value, group=variable)) +
  geom_point(aes(color=variable)) +
  scale_color_brewer(palette="Dark2") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab(xlab) +
  ylab("delta Log-Likelihood") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 3))

# Display plots
(MC.raw + theme(legend.position = "none")) + cowplot::get_legend(MC.raw) + plot_layout(ncol = 2, widths = c(8,2))
```

![](plot_nQuire_delta_loglikelihood_files/figure-html/Plot_Mcapitata-1.png)<!-- -->

```r
(MC.denoised + theme(legend.position = "none")) + cowplot::get_legend(MC.denoised) + plot_layout(ncol = 2, widths = c(8,2)) 
```

![](plot_nQuire_delta_loglikelihood_files/figure-html/Plot_Mcapitata-2.png)<!-- -->





## *P. acuta* RNA-Seq samples from this study


```r
xlab <- "P. acuta RNA-seq samples"

file.name.1 <- "../../../Pocillopora_acuta/02_ploidy_analysis/nQuire/nQuire_results.tsv"
data.1 <- read.table(file.name.1, header = T, sep = '\t')
file.name.2 <- "../../../samples_from_SRA/Pocillopora/02_ploidy_analysis/nQuire/nQuire_results.tsv"
data.2 <- read.table(file.name.2, header = T, sep = '\t')

data <- rbind(data.1, data.2)

data.melted <- data %>% 
  select(c("ID", "Ploidy", "Diploid", "Triploid", "Tetraploid")) %>% 
  melt(id.vars = c("ID", "Ploidy")) %>%
  mutate(ID = factor(ID, data$ID)) %>%
  mutate(variable = factor(variable, c("Diploid", "Triploid", "Tetraploid")))

data.denoised.melted <- data %>% 
  select(c("ID", "Ploidy", "Diploid_denoised", "Triploid_denoised", "Tetraploid_denoised")) %>% 
  melt(id.vars = c("ID", "Ploidy")) %>%
  mutate(ID = factor(ID, data$ID)) %>%
  mutate(variable = factor(variable, c("Diploid_denoised", "Triploid_denoised", "Tetraploid_denoised")))

# Plot raw sites
PA.raw <- ggplot(data.melted, aes(x=ID, y=value, group=variable)) +
  geom_point(aes(color=variable)) +
  scale_color_brewer(palette="Dark2") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab(xlab) +
  ylab("delta Log-Likelihood") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 3))

# Plot denoised sites
PA.denoised <- ggplot(data.denoised.melted, aes(x=ID, y=value, group=variable)) +
  geom_point(aes(color=variable)) +
  scale_color_brewer(palette="Dark2") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab(xlab) +
  ylab("delta Log-Likelihood") +
  theme(panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size = 3))

# Display plots
(PA.raw + theme(legend.position = "none")) + cowplot::get_legend(PA.raw) + plot_layout(ncol = 2, widths = c(8,2))
```

![](plot_nQuire_delta_loglikelihood_files/figure-html/Plot_Pacuta-1.png)<!-- -->

```r
(PA.denoised + theme(legend.position = "none")) + cowplot::get_legend(PA.denoised) + plot_layout(ncol = 2, widths = c(8,2)) 
```

![](plot_nQuire_delta_loglikelihood_files/figure-html/Plot_Pacuta-2.png)<!-- -->





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
## [1] cowplot_1.1.1      patchwork_1.1.2    RColorBrewer_1.1-3 ggplot2_3.3.6     
## [5] reshape2_1.4.4     dplyr_1.0.10       tibble_3.1.8      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.9       highr_0.9        pillar_1.8.1     bslib_0.4.0     
##  [5] compiler_4.1.2   jquerylib_0.1.4  plyr_1.8.7       tools_4.1.2     
##  [9] digest_0.6.29    gtable_0.3.1     jsonlite_1.8.0   evaluate_0.16   
## [13] lifecycle_1.0.1  pkgconfig_2.0.3  rlang_1.0.5      cli_3.3.0       
## [17] DBI_1.1.3        rstudioapi_0.14  yaml_2.3.5       xfun_0.32       
## [21] fastmap_1.1.0    withr_2.5.0      stringr_1.4.1    knitr_1.40      
## [25] generics_0.1.3   vctrs_0.4.1      sass_0.4.2       grid_4.1.2      
## [29] tidyselect_1.1.2 glue_1.6.2       R6_2.5.1         fansi_1.0.3     
## [33] rmarkdown_2.16   farver_2.1.1     purrr_0.3.4      magrittr_2.0.3  
## [37] ellipsis_0.3.2   scales_1.2.1     htmltools_0.5.3  assertthat_0.2.1
## [41] colorspace_2.0-3 labeling_0.4.2   utf8_1.2.2       stringi_1.7.8   
## [45] munsell_0.5.0    cachem_1.0.6
```
