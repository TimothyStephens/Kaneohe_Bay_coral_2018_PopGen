---
title: "Plot `vcftools --relatedness2` results for *M. capitata* RNA-seq samples from this study and SRA"
author: "Timothy Stephens"
date: "26/09/2022"
output: 
  html_document:
    keep_md: yes
---



# Setup

Setup R env. Load packages and set default image export formats, size and resolution.

```r
knitr::opts_chunk$set(echo = TRUE,
                      fig.height = 12, 
                      fig.width = 12, 
                      dev = c("png", "pdf"),
                      dpi = 1000)
library(readxl)
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
library(ggplot2)
library(gplots)
```

```
## 
## Attaching package: 'gplots'
```

```
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
library(reshape2)
library(RcppCNPy)
library(tibble)
library(ggdendro)
library(cowplot)
library(RColorBrewer)
library(viridis)
```

```
## Loading required package: viridisLite
```

```r
library(knitr)
library(plotly)
```

```
## 
## Attaching package: 'plotly'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     last_plot
```

```
## The following object is masked from 'package:stats':
## 
##     filter
```

```
## The following object is masked from 'package:graphics':
## 
##     layout
```

```r
library(phylogram)
options(scipen = 999) #Prevent scientific notation
cexSize <- 0.3
```





# Metadata 

Load file with annotation for each sample. 

```r
samples.info <- read.table("../../../samples_Mcapitata_ALL.annotations.txt", header=T, comment.char='')
rownames(samples.info) <- samples.info$sample
samples.info
```

```
##                                            sample   species
## Mcapitata_ATAC_TP11_2302 Mcapitata_ATAC_TP11_2302 Mcapitata
## Mcapitata_HTAC_TP11_2380 Mcapitata_HTAC_TP11_2380 Mcapitata
## Mcapitata_ATAC_TP6_2402   Mcapitata_ATAC_TP6_2402 Mcapitata
## Mcapitata_ATAC_TP12_2403 Mcapitata_ATAC_TP12_2403 Mcapitata
## Mcapitata_ATAC_TP1_1037   Mcapitata_ATAC_TP1_1037 Mcapitata
## Mcapitata_ATAC_TP1_1600   Mcapitata_ATAC_TP1_1600 Mcapitata
## Mcapitata_ATAC_TP1_1652   Mcapitata_ATAC_TP1_1652 Mcapitata
## Mcapitata_ATAC_TP3_1101   Mcapitata_ATAC_TP3_1101 Mcapitata
## Mcapitata_ATAC_TP3_1548   Mcapitata_ATAC_TP3_1548 Mcapitata
## Mcapitata_ATAC_TP3_1628   Mcapitata_ATAC_TP3_1628 Mcapitata
## Mcapitata_ATAC_TP4_1108   Mcapitata_ATAC_TP4_1108 Mcapitata
## Mcapitata_ATAC_TP4_1609   Mcapitata_ATAC_TP4_1609 Mcapitata
## Mcapitata_ATAC_TP4_1651   Mcapitata_ATAC_TP4_1651 Mcapitata
## Mcapitata_ATAC_TP5_1196   Mcapitata_ATAC_TP5_1196 Mcapitata
## Mcapitata_ATAC_TP5_1610   Mcapitata_ATAC_TP5_1610 Mcapitata
## Mcapitata_ATAC_TP5_1776   Mcapitata_ATAC_TP5_1776 Mcapitata
## Mcapitata_ATAC_TP6_1114   Mcapitata_ATAC_TP6_1114 Mcapitata
## Mcapitata_ATAC_TP6_1611   Mcapitata_ATAC_TP6_1611 Mcapitata
## Mcapitata_ATAC_TP7_1058   Mcapitata_ATAC_TP7_1058 Mcapitata
## Mcapitata_ATAC_TP7_1455   Mcapitata_ATAC_TP7_1455 Mcapitata
## Mcapitata_ATAC_TP7_1499   Mcapitata_ATAC_TP7_1499 Mcapitata
## Mcapitata_ATAC_TP8_1083   Mcapitata_ATAC_TP8_1083 Mcapitata
## Mcapitata_ATAC_TP8_1436   Mcapitata_ATAC_TP8_1436 Mcapitata
## Mcapitata_ATAC_TP8_1779   Mcapitata_ATAC_TP8_1779 Mcapitata
## Mcapitata_ATAC_TP9_1121   Mcapitata_ATAC_TP9_1121 Mcapitata
## Mcapitata_ATAC_TP9_1420   Mcapitata_ATAC_TP9_1420 Mcapitata
## Mcapitata_ATAC_TP9_1580   Mcapitata_ATAC_TP9_1580 Mcapitata
## Mcapitata_ATAC_TP10_1095 Mcapitata_ATAC_TP10_1095 Mcapitata
## Mcapitata_ATAC_TP10_1561 Mcapitata_ATAC_TP10_1561 Mcapitata
## Mcapitata_ATAC_TP10_1631 Mcapitata_ATAC_TP10_1631 Mcapitata
## Mcapitata_ATAC_TP11_1076 Mcapitata_ATAC_TP11_1076 Mcapitata
## Mcapitata_ATAC_TP11_1644 Mcapitata_ATAC_TP11_1644 Mcapitata
## Mcapitata_ATAC_TP12_1120 Mcapitata_ATAC_TP12_1120 Mcapitata
## Mcapitata_ATAC_TP12_1452 Mcapitata_ATAC_TP12_1452 Mcapitata
## Mcapitata_ATHC_TP1_1218   Mcapitata_ATHC_TP1_1218 Mcapitata
## Mcapitata_ATHC_TP1_1826   Mcapitata_ATHC_TP1_1826 Mcapitata
## Mcapitata_ATHC_TP1_2068   Mcapitata_ATHC_TP1_2068 Mcapitata
## Mcapitata_ATHC_TP3_1544   Mcapitata_ATHC_TP3_1544 Mcapitata
## Mcapitata_ATHC_TP3_2731   Mcapitata_ATHC_TP3_2731 Mcapitata
## Mcapitata_ATHC_TP3_2866   Mcapitata_ATHC_TP3_2866 Mcapitata
## Mcapitata_ATHC_TP4_1221   Mcapitata_ATHC_TP4_1221 Mcapitata
## Mcapitata_ATHC_TP4_2561   Mcapitata_ATHC_TP4_2561 Mcapitata
## Mcapitata_ATHC_TP4_2734   Mcapitata_ATHC_TP4_2734 Mcapitata
## Mcapitata_ATHC_TP5_1229   Mcapitata_ATHC_TP5_1229 Mcapitata
## Mcapitata_ATHC_TP5_1706   Mcapitata_ATHC_TP5_1706 Mcapitata
## Mcapitata_ATHC_TP5_2986   Mcapitata_ATHC_TP5_2986 Mcapitata
## Mcapitata_ATHC_TP6_1212   Mcapitata_ATHC_TP6_1212 Mcapitata
## Mcapitata_ATHC_TP6_2016   Mcapitata_ATHC_TP6_2016 Mcapitata
## Mcapitata_ATHC_TP6_2555   Mcapitata_ATHC_TP6_2555 Mcapitata
## Mcapitata_ATHC_TP7_1223   Mcapitata_ATHC_TP7_1223 Mcapitata
## Mcapitata_ATHC_TP7_2860   Mcapitata_ATHC_TP7_2860 Mcapitata
## Mcapitata_ATHC_TP7_2875   Mcapitata_ATHC_TP7_2875 Mcapitata
## Mcapitata_ATHC_TP8_1260   Mcapitata_ATHC_TP8_1260 Mcapitata
## Mcapitata_ATHC_TP8_2735   Mcapitata_ATHC_TP8_2735 Mcapitata
## Mcapitata_ATHC_TP8_2753   Mcapitata_ATHC_TP8_2753 Mcapitata
## Mcapitata_ATHC_TP9_1148   Mcapitata_ATHC_TP9_1148 Mcapitata
## Mcapitata_ATHC_TP9_2862   Mcapitata_ATHC_TP9_2862 Mcapitata
## Mcapitata_ATHC_TP9_2995   Mcapitata_ATHC_TP9_2995 Mcapitata
## Mcapitata_ATHC_TP10_1204 Mcapitata_ATHC_TP10_1204 Mcapitata
## Mcapitata_ATHC_TP10_2554 Mcapitata_ATHC_TP10_2554 Mcapitata
## Mcapitata_ATHC_TP10_2737 Mcapitata_ATHC_TP10_2737 Mcapitata
## Mcapitata_ATHC_TP11_1237 Mcapitata_ATHC_TP11_1237 Mcapitata
## Mcapitata_ATHC_TP11_2188 Mcapitata_ATHC_TP11_2188 Mcapitata
## Mcapitata_ATHC_TP11_2756 Mcapitata_ATHC_TP11_2756 Mcapitata
## Mcapitata_ATHC_TP12_1154 Mcapitata_ATHC_TP12_1154 Mcapitata
## Mcapitata_ATHC_TP12_2736 Mcapitata_ATHC_TP12_2736 Mcapitata
## Mcapitata_ATHC_TP12_2990 Mcapitata_ATHC_TP12_2990 Mcapitata
## Mcapitata_HTAC_TP1_1579   Mcapitata_HTAC_TP1_1579 Mcapitata
## Mcapitata_HTAC_TP1_2153   Mcapitata_HTAC_TP1_2153 Mcapitata
## Mcapitata_HTAC_TP1_2183   Mcapitata_HTAC_TP1_2183 Mcapitata
## Mcapitata_HTAC_TP3_1289   Mcapitata_HTAC_TP3_1289 Mcapitata
## Mcapitata_HTAC_TP3_1751   Mcapitata_HTAC_TP3_1751 Mcapitata
## Mcapitata_HTAC_TP3_2021   Mcapitata_HTAC_TP3_2021 Mcapitata
## Mcapitata_HTAC_TP4_1269   Mcapitata_HTAC_TP4_1269 Mcapitata
## Mcapitata_HTAC_TP4_1481   Mcapitata_HTAC_TP4_1481 Mcapitata
## Mcapitata_HTAC_TP4_2000   Mcapitata_HTAC_TP4_2000 Mcapitata
## Mcapitata_HTAC_TP5_1321   Mcapitata_HTAC_TP5_1321 Mcapitata
## Mcapitata_HTAC_TP5_1583   Mcapitata_HTAC_TP5_1583 Mcapitata
## Mcapitata_HTAC_TP5_1997   Mcapitata_HTAC_TP5_1997 Mcapitata
## Mcapitata_HTAC_TP6_1496   Mcapitata_HTAC_TP6_1496 Mcapitata
## Mcapitata_HTAC_TP6_1588   Mcapitata_HTAC_TP6_1588 Mcapitata
## Mcapitata_HTAC_TP6_1705   Mcapitata_HTAC_TP6_1705 Mcapitata
## Mcapitata_HTAC_TP7_1278   Mcapitata_HTAC_TP7_1278 Mcapitata
## Mcapitata_HTAC_TP7_1645   Mcapitata_HTAC_TP7_1645 Mcapitata
## Mcapitata_HTAC_TP7_1722   Mcapitata_HTAC_TP7_1722 Mcapitata
## Mcapitata_HTAC_TP8_1235   Mcapitata_HTAC_TP8_1235 Mcapitata
## Mcapitata_HTAC_TP8_2386   Mcapitata_HTAC_TP8_2386 Mcapitata
## Mcapitata_HTAC_TP8_2410   Mcapitata_HTAC_TP8_2410 Mcapitata
## Mcapitata_HTAC_TP9_1306   Mcapitata_HTAC_TP9_1306 Mcapitata
## Mcapitata_HTAC_TP9_1467   Mcapitata_HTAC_TP9_1467 Mcapitata
## Mcapitata_HTAC_TP9_2412   Mcapitata_HTAC_TP9_2412 Mcapitata
## Mcapitata_HTAC_TP10_1315 Mcapitata_HTAC_TP10_1315 Mcapitata
## Mcapitata_HTAC_TP10_1478 Mcapitata_HTAC_TP10_1478 Mcapitata
## Mcapitata_HTAC_TP10_1754 Mcapitata_HTAC_TP10_1754 Mcapitata
## Mcapitata_HTAC_TP11_1248 Mcapitata_HTAC_TP11_1248 Mcapitata
## Mcapitata_HTAC_TP11_1562 Mcapitata_HTAC_TP11_1562 Mcapitata
## Mcapitata_HTAC_TP12_1632 Mcapitata_HTAC_TP12_1632 Mcapitata
## Mcapitata_HTAC_TP12_1729 Mcapitata_HTAC_TP12_1729 Mcapitata
## Mcapitata_HTAC_TP12_2007 Mcapitata_HTAC_TP12_2007 Mcapitata
## Mcapitata_HTHC_TP1_1145   Mcapitata_HTHC_TP1_1145 Mcapitata
## Mcapitata_HTHC_TP1_1323   Mcapitata_HTHC_TP1_1323 Mcapitata
## Mcapitata_HTHC_TP1_2081   Mcapitata_HTHC_TP1_2081 Mcapitata
## Mcapitata_HTHC_TP3_1128   Mcapitata_HTHC_TP3_1128 Mcapitata
## Mcapitata_HTHC_TP3_1277   Mcapitata_HTHC_TP3_1277 Mcapitata
## Mcapitata_HTHC_TP3_2518   Mcapitata_HTHC_TP3_2518 Mcapitata
## Mcapitata_HTHC_TP4_1124   Mcapitata_HTHC_TP4_1124 Mcapitata
## Mcapitata_HTHC_TP4_1328   Mcapitata_HTHC_TP4_1328 Mcapitata
## Mcapitata_HTHC_TP4_2204   Mcapitata_HTHC_TP4_2204 Mcapitata
## Mcapitata_HTHC_TP5_1345   Mcapitata_HTHC_TP5_1345 Mcapitata
## Mcapitata_HTHC_TP5_1449   Mcapitata_HTHC_TP5_1449 Mcapitata
## Mcapitata_HTHC_TP5_1694   Mcapitata_HTHC_TP5_1694 Mcapitata
## Mcapitata_HTHC_TP6_1164   Mcapitata_HTHC_TP6_1164 Mcapitata
## Mcapitata_HTHC_TP6_1317   Mcapitata_HTHC_TP6_1317 Mcapitata
## Mcapitata_HTHC_TP6_1604   Mcapitata_HTHC_TP6_1604 Mcapitata
## Mcapitata_HTHC_TP7_1126   Mcapitata_HTHC_TP7_1126 Mcapitata
## Mcapitata_HTHC_TP7_1250   Mcapitata_HTHC_TP7_1250 Mcapitata
## Mcapitata_HTHC_TP7_2419   Mcapitata_HTHC_TP7_2419 Mcapitata
## Mcapitata_HTHC_TP8_1082   Mcapitata_HTHC_TP8_1082 Mcapitata
## Mcapitata_HTHC_TP8_1246   Mcapitata_HTHC_TP8_1246 Mcapitata
## Mcapitata_HTHC_TP8_2067   Mcapitata_HTHC_TP8_2067 Mcapitata
## Mcapitata_HTHC_TP9_1078   Mcapitata_HTHC_TP9_1078 Mcapitata
## Mcapitata_HTHC_TP9_1331   Mcapitata_HTHC_TP9_1331 Mcapitata
## Mcapitata_HTHC_TP9_2009   Mcapitata_HTHC_TP9_2009 Mcapitata
## Mcapitata_HTHC_TP10_1074 Mcapitata_HTHC_TP10_1074 Mcapitata
## Mcapitata_HTHC_TP10_1332 Mcapitata_HTHC_TP10_1332 Mcapitata
## Mcapitata_HTHC_TP10_1689 Mcapitata_HTHC_TP10_1689 Mcapitata
## Mcapitata_HTHC_TP11_1178 Mcapitata_HTHC_TP11_1178 Mcapitata
## Mcapitata_HTHC_TP11_1270 Mcapitata_HTHC_TP11_1270 Mcapitata
## Mcapitata_HTHC_TP11_2511 Mcapitata_HTHC_TP11_2511 Mcapitata
## Mcapitata_HTHC_TP12_1140 Mcapitata_HTHC_TP12_1140 Mcapitata
## Mcapitata_HTHC_TP12_1274 Mcapitata_HTHC_TP12_1274 Mcapitata
## Mcapitata_HTHC_TP12_2190 Mcapitata_HTHC_TP12_2190 Mcapitata
## SRR5453739                             SRR5453739 Mcapitata
## SRR5453740                             SRR5453740 Mcapitata
## SRR5453741                             SRR5453741 Mcapitata
## SRR5453742                             SRR5453742 Mcapitata
## SRR5453743                             SRR5453743 Mcapitata
## SRR5453744                             SRR5453744 Mcapitata
## SRR5453745                             SRR5453745 Mcapitata
## SRR5453746                             SRR5453746 Mcapitata
## SRR5453747                             SRR5453747 Mcapitata
## SRR5453748                             SRR5453748 Mcapitata
## SRR5453749                             SRR5453749 Mcapitata
## SRR5453750                             SRR5453750 Mcapitata
## SRR5453751                             SRR5453751 Mcapitata
## SRR5453752                             SRR5453752 Mcapitata
## SRR5453753                             SRR5453753 Mcapitata
## SRR5453754                             SRR5453754 Mcapitata
## SRR5453755                             SRR5453755 Mcapitata
## SRR5453756                             SRR5453756 Mcapitata
## SRR5453757                             SRR5453757 Mcapitata
## SRR5453758                             SRR5453758 Mcapitata
## SRR5453759                             SRR5453759 Mcapitata
## SRR5453760                             SRR5453760 Mcapitata
## SRR5453761                             SRR5453761 Mcapitata
## SRR5453762                             SRR5453762 Mcapitata
## SRR5453763                             SRR5453763 Mcapitata
## SRR5453764                             SRR5453764 Mcapitata
## SRR5453765                             SRR5453765 Mcapitata
##                                                            treatment timepoint
## Mcapitata_ATAC_TP11_2302                                        ATAC      TP11
## Mcapitata_HTAC_TP11_2380                                        HTAC      TP11
## Mcapitata_ATAC_TP6_2402                                         ATAC       TP6
## Mcapitata_ATAC_TP12_2403                                        ATAC      TP12
## Mcapitata_ATAC_TP1_1037                                         ATAC       TP1
## Mcapitata_ATAC_TP1_1600                                         ATAC       TP1
## Mcapitata_ATAC_TP1_1652                                         ATAC       TP1
## Mcapitata_ATAC_TP3_1101                                         ATAC       TP3
## Mcapitata_ATAC_TP3_1548                                         ATAC       TP3
## Mcapitata_ATAC_TP3_1628                                         ATAC       TP3
## Mcapitata_ATAC_TP4_1108                                         ATAC       TP4
## Mcapitata_ATAC_TP4_1609                                         ATAC       TP4
## Mcapitata_ATAC_TP4_1651                                         ATAC       TP4
## Mcapitata_ATAC_TP5_1196                                         ATAC       TP5
## Mcapitata_ATAC_TP5_1610                                         ATAC       TP5
## Mcapitata_ATAC_TP5_1776                                         ATAC       TP5
## Mcapitata_ATAC_TP6_1114                                         ATAC       TP6
## Mcapitata_ATAC_TP6_1611                                         ATAC       TP6
## Mcapitata_ATAC_TP7_1058                                         ATAC       TP7
## Mcapitata_ATAC_TP7_1455                                         ATAC       TP7
## Mcapitata_ATAC_TP7_1499                                         ATAC       TP7
## Mcapitata_ATAC_TP8_1083                                         ATAC       TP8
## Mcapitata_ATAC_TP8_1436                                         ATAC       TP8
## Mcapitata_ATAC_TP8_1779                                         ATAC       TP8
## Mcapitata_ATAC_TP9_1121                                         ATAC       TP9
## Mcapitata_ATAC_TP9_1420                                         ATAC       TP9
## Mcapitata_ATAC_TP9_1580                                         ATAC       TP9
## Mcapitata_ATAC_TP10_1095                                        ATAC      TP10
## Mcapitata_ATAC_TP10_1561                                        ATAC      TP10
## Mcapitata_ATAC_TP10_1631                                        ATAC      TP10
## Mcapitata_ATAC_TP11_1076                                        ATAC      TP11
## Mcapitata_ATAC_TP11_1644                                        ATAC      TP11
## Mcapitata_ATAC_TP12_1120                                        ATAC      TP12
## Mcapitata_ATAC_TP12_1452                                        ATAC      TP12
## Mcapitata_ATHC_TP1_1218                                         ATHC       TP1
## Mcapitata_ATHC_TP1_1826                                         ATHC       TP1
## Mcapitata_ATHC_TP1_2068                                         ATHC       TP1
## Mcapitata_ATHC_TP3_1544                                         ATHC       TP3
## Mcapitata_ATHC_TP3_2731                                         ATHC       TP3
## Mcapitata_ATHC_TP3_2866                                         ATHC       TP3
## Mcapitata_ATHC_TP4_1221                                         ATHC       TP4
## Mcapitata_ATHC_TP4_2561                                         ATHC       TP4
## Mcapitata_ATHC_TP4_2734                                         ATHC       TP4
## Mcapitata_ATHC_TP5_1229                                         ATHC       TP5
## Mcapitata_ATHC_TP5_1706                                         ATHC       TP5
## Mcapitata_ATHC_TP5_2986                                         ATHC       TP5
## Mcapitata_ATHC_TP6_1212                                         ATHC       TP6
## Mcapitata_ATHC_TP6_2016                                         ATHC       TP6
## Mcapitata_ATHC_TP6_2555                                         ATHC       TP6
## Mcapitata_ATHC_TP7_1223                                         ATHC       TP7
## Mcapitata_ATHC_TP7_2860                                         ATHC       TP7
## Mcapitata_ATHC_TP7_2875                                         ATHC       TP7
## Mcapitata_ATHC_TP8_1260                                         ATHC       TP8
## Mcapitata_ATHC_TP8_2735                                         ATHC       TP8
## Mcapitata_ATHC_TP8_2753                                         ATHC       TP8
## Mcapitata_ATHC_TP9_1148                                         ATHC       TP9
## Mcapitata_ATHC_TP9_2862                                         ATHC       TP9
## Mcapitata_ATHC_TP9_2995                                         ATHC       TP9
## Mcapitata_ATHC_TP10_1204                                        ATHC      TP10
## Mcapitata_ATHC_TP10_2554                                        ATHC      TP10
## Mcapitata_ATHC_TP10_2737                                        ATHC      TP10
## Mcapitata_ATHC_TP11_1237                                        ATHC      TP11
## Mcapitata_ATHC_TP11_2188                                        ATHC      TP11
## Mcapitata_ATHC_TP11_2756                                        ATHC      TP11
## Mcapitata_ATHC_TP12_1154                                        ATHC      TP12
## Mcapitata_ATHC_TP12_2736                                        ATHC      TP12
## Mcapitata_ATHC_TP12_2990                                        ATHC      TP12
## Mcapitata_HTAC_TP1_1579                                         HTAC       TP1
## Mcapitata_HTAC_TP1_2153                                         HTAC       TP1
## Mcapitata_HTAC_TP1_2183                                         HTAC       TP1
## Mcapitata_HTAC_TP3_1289                                         HTAC       TP3
## Mcapitata_HTAC_TP3_1751                                         HTAC       TP3
## Mcapitata_HTAC_TP3_2021                                         HTAC       TP3
## Mcapitata_HTAC_TP4_1269                                         HTAC       TP4
## Mcapitata_HTAC_TP4_1481                                         HTAC       TP4
## Mcapitata_HTAC_TP4_2000                                         HTAC       TP4
## Mcapitata_HTAC_TP5_1321                                         HTAC       TP5
## Mcapitata_HTAC_TP5_1583                                         HTAC       TP5
## Mcapitata_HTAC_TP5_1997                                         HTAC       TP5
## Mcapitata_HTAC_TP6_1496                                         HTAC       TP6
## Mcapitata_HTAC_TP6_1588                                         HTAC       TP6
## Mcapitata_HTAC_TP6_1705                                         HTAC       TP6
## Mcapitata_HTAC_TP7_1278                                         HTAC       TP7
## Mcapitata_HTAC_TP7_1645                                         HTAC       TP7
## Mcapitata_HTAC_TP7_1722                                         HTAC       TP7
## Mcapitata_HTAC_TP8_1235                                         HTAC       TP8
## Mcapitata_HTAC_TP8_2386                                         HTAC       TP8
## Mcapitata_HTAC_TP8_2410                                         HTAC       TP8
## Mcapitata_HTAC_TP9_1306                                         HTAC       TP9
## Mcapitata_HTAC_TP9_1467                                         HTAC       TP9
## Mcapitata_HTAC_TP9_2412                                         HTAC       TP9
## Mcapitata_HTAC_TP10_1315                                        HTAC      TP10
## Mcapitata_HTAC_TP10_1478                                        HTAC      TP10
## Mcapitata_HTAC_TP10_1754                                        HTAC      TP10
## Mcapitata_HTAC_TP11_1248                                        HTAC      TP11
## Mcapitata_HTAC_TP11_1562                                        HTAC      TP11
## Mcapitata_HTAC_TP12_1632                                        HTAC      TP12
## Mcapitata_HTAC_TP12_1729                                        HTAC      TP12
## Mcapitata_HTAC_TP12_2007                                        HTAC      TP12
## Mcapitata_HTHC_TP1_1145                                         HTHC       TP1
## Mcapitata_HTHC_TP1_1323                                         HTHC       TP1
## Mcapitata_HTHC_TP1_2081                                         HTHC       TP1
## Mcapitata_HTHC_TP3_1128                                         HTHC       TP3
## Mcapitata_HTHC_TP3_1277                                         HTHC       TP3
## Mcapitata_HTHC_TP3_2518                                         HTHC       TP3
## Mcapitata_HTHC_TP4_1124                                         HTHC       TP4
## Mcapitata_HTHC_TP4_1328                                         HTHC       TP4
## Mcapitata_HTHC_TP4_2204                                         HTHC       TP4
## Mcapitata_HTHC_TP5_1345                                         HTHC       TP5
## Mcapitata_HTHC_TP5_1449                                         HTHC       TP5
## Mcapitata_HTHC_TP5_1694                                         HTHC       TP5
## Mcapitata_HTHC_TP6_1164                                         HTHC       TP6
## Mcapitata_HTHC_TP6_1317                                         HTHC       TP6
## Mcapitata_HTHC_TP6_1604                                         HTHC       TP6
## Mcapitata_HTHC_TP7_1126                                         HTHC       TP7
## Mcapitata_HTHC_TP7_1250                                         HTHC       TP7
## Mcapitata_HTHC_TP7_2419                                         HTHC       TP7
## Mcapitata_HTHC_TP8_1082                                         HTHC       TP8
## Mcapitata_HTHC_TP8_1246                                         HTHC       TP8
## Mcapitata_HTHC_TP8_2067                                         HTHC       TP8
## Mcapitata_HTHC_TP9_1078                                         HTHC       TP9
## Mcapitata_HTHC_TP9_1331                                         HTHC       TP9
## Mcapitata_HTHC_TP9_2009                                         HTHC       TP9
## Mcapitata_HTHC_TP10_1074                                        HTHC      TP10
## Mcapitata_HTHC_TP10_1332                                        HTHC      TP10
## Mcapitata_HTHC_TP10_1689                                        HTHC      TP10
## Mcapitata_HTHC_TP11_1178                                        HTHC      TP11
## Mcapitata_HTHC_TP11_1270                                        HTHC      TP11
## Mcapitata_HTHC_TP11_2511                                        HTHC      TP11
## Mcapitata_HTHC_TP12_1140                                        HTHC      TP12
## Mcapitata_HTHC_TP12_1274                                        HTHC      TP12
## Mcapitata_HTHC_TP12_2190                                        HTHC      TP12
## SRR5453739               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453740               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453741               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453742                             healthy_colony_healthy_tissue      <NA>
## SRR5453743                             healthy_colony_healthy_tissue      <NA>
## SRR5453744                             healthy_colony_healthy_tissue      <NA>
## SRR5453745                      growth_anomaly_colony_healthy_tissue      <NA>
## SRR5453746                      growth_anomaly_colony_healthy_tissue      <NA>
## SRR5453747                      growth_anomaly_colony_healthy_tissue      <NA>
## SRR5453748               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453749               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453750               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453751               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453752               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453753               growth_anomaly_colony_growth_anomaly_tissue      <NA>
## SRR5453754                             healthy_colony_healthy_tissue      <NA>
## SRR5453755                             healthy_colony_healthy_tissue      <NA>
## SRR5453756                             healthy_colony_healthy_tissue      <NA>
## SRR5453757                             healthy_colony_healthy_tissue      <NA>
## SRR5453758                             healthy_colony_healthy_tissue      <NA>
## SRR5453759                             healthy_colony_healthy_tissue      <NA>
## SRR5453760                      growth_anomaly_colony_healthy_tissue      <NA>
## SRR5453761                      growth_anomaly_colony_healthy_tissue      <NA>
## SRR5453762                      growth_anomaly_colony_healthy_tissue      <NA>
## SRR5453763                      growth_anomaly_colony_healthy_tissue      <NA>
## SRR5453764                      growth_anomaly_colony_healthy_tissue      <NA>
## SRR5453765                      growth_anomaly_colony_healthy_tissue      <NA>
##                             plugid            reef reef_color ploidy
## Mcapitata_ATAC_TP11_2302      2302         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP11_2380      2380         Reef.18    #d8daeb      2
## Mcapitata_ATAC_TP6_2402       2402         Reef.18    #d8daeb      2
## Mcapitata_ATAC_TP12_2403      2403         Reef.18    #d8daeb      2
## Mcapitata_ATAC_TP1_1037       1037      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP1_1600       1600      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP1_1652       1652            HIMB    #b35806      2
## Mcapitata_ATAC_TP3_1101       1101            HIMB    #b35806      2
## Mcapitata_ATAC_TP3_1548       1548      Reef.11.13    #fee0b6      2
## Mcapitata_ATAC_TP3_1628       1628            HIMB    #b35806      2
## Mcapitata_ATAC_TP4_1108       1108      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP4_1609       1609            HIMB    #b35806      2
## Mcapitata_ATAC_TP4_1651       1651            HIMB    #b35806      2
## Mcapitata_ATAC_TP5_1196       1196         Reef.18    #d8daeb      2
## Mcapitata_ATAC_TP5_1610       1610      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP5_1776       1776      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP6_1114       1114      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP6_1611       1611      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP7_1058       1058 Lilipuna.Fringe    #f1a340      2
## Mcapitata_ATAC_TP7_1455       1455      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP7_1499       1499      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP8_1083       1083      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP8_1436       1436      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP8_1779       1779      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP9_1121       1121            HIMB    #b35806      2
## Mcapitata_ATAC_TP9_1420       1420      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP9_1580       1580      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP10_1095      1095      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP10_1561      1561      Reef.11.13    #fee0b6      2
## Mcapitata_ATAC_TP10_1631      1631            HIMB    #b35806      2
## Mcapitata_ATAC_TP11_1076      1076      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP11_1644      1644            HIMB    #b35806      2
## Mcapitata_ATAC_TP12_1120      1120      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP12_1452      1452      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP1_1218       1218      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP1_1826       1826      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP1_2068       2068      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP3_1544       1544         Reef.18    #d8daeb      2
## Mcapitata_ATHC_TP3_2731       2731      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP3_2866       2866            HIMB    #b35806      2
## Mcapitata_ATHC_TP4_1221       1221      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP4_2561       2561         Reef.18    #d8daeb      2
## Mcapitata_ATHC_TP4_2734       2734      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP5_1229       1229      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP5_1706       1706      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP5_2986       2986      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP6_1212       1212      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP6_2016       2016      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP6_2555       2555         Reef.18    #d8daeb      2
## Mcapitata_ATHC_TP7_1223       1223      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP7_2860       2860            HIMB    #b35806      2
## Mcapitata_ATHC_TP7_2875       2875            HIMB    #b35806      2
## Mcapitata_ATHC_TP8_1260       1260      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP8_2735       2735      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP8_2753       2753      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP9_1148       1148      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP9_2862       2862            HIMB    #b35806      2
## Mcapitata_ATHC_TP9_2995       2995      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP10_1204      1204      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP10_2554      2554         Reef.18    #d8daeb      2
## Mcapitata_ATHC_TP10_2737      2737      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP11_1237      1237      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP11_2188      2188      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP11_2756      2756      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP12_1154      1154            HIMB    #b35806      2
## Mcapitata_ATHC_TP12_2736      2736      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP12_2990      2990         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP1_1579       1579 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTAC_TP1_2153       2153         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP1_2183       2183      Reef.42.43    #542788      2
## Mcapitata_HTAC_TP3_1289       1289      Reef.35.36    #998ec3      2
## Mcapitata_HTAC_TP3_1751       1751      Reef.42.43    #542788      2
## Mcapitata_HTAC_TP3_2021       2021      Reef.42.43    #542788      2
## Mcapitata_HTAC_TP4_1269       1269      Reef.35.36    #998ec3      2
## Mcapitata_HTAC_TP4_1481       1481      Reef.35.36    #998ec3      2
## Mcapitata_HTAC_TP4_2000       2000            HIMB    #b35806      2
## Mcapitata_HTAC_TP5_1321       1321      Reef.11.13    #fee0b6      2
## Mcapitata_HTAC_TP5_1583       1583 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTAC_TP5_1997       1997            HIMB    #b35806      2
## Mcapitata_HTAC_TP6_1496       1496      Reef.35.36    #998ec3      2
## Mcapitata_HTAC_TP6_1588       1588         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP6_1705       1705      Reef.35.36    #998ec3      2
## Mcapitata_HTAC_TP7_1278       1278      Reef.11.13    #fee0b6      2
## Mcapitata_HTAC_TP7_1645       1645      Reef.42.43    #542788      2
## Mcapitata_HTAC_TP7_1722       1722      Reef.35.36    #998ec3      2
## Mcapitata_HTAC_TP8_1235       1235      Reef.42.43    #542788      2
## Mcapitata_HTAC_TP8_2386       2386         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP8_2410       2410         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP9_1306       1306 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTAC_TP9_1467       1467      Reef.11.13    #fee0b6      2
## Mcapitata_HTAC_TP9_2412       2412         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP10_1315      1315 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTAC_TP10_1478      1478      Reef.11.13    #fee0b6      2
## Mcapitata_HTAC_TP10_1754      1754            HIMB    #b35806      2
## Mcapitata_HTAC_TP11_1248      1248            HIMB    #b35806      2
## Mcapitata_HTAC_TP11_1562      1562      Reef.11.13    #fee0b6      2
## Mcapitata_HTAC_TP12_1632      1632            HIMB    #b35806      4
## Mcapitata_HTAC_TP12_1729      1729         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP12_2007      2007      Reef.42.43    #542788      2
## Mcapitata_HTHC_TP1_1145       1145            HIMB    #b35806      2
## Mcapitata_HTHC_TP1_1323       1323      Reef.35.36    #998ec3      2
## Mcapitata_HTHC_TP1_2081       2081            HIMB    #b35806      2
## Mcapitata_HTHC_TP3_1128       1128            HIMB    #b35806      2
## Mcapitata_HTHC_TP3_1277       1277            HIMB    #b35806      2
## Mcapitata_HTHC_TP3_2518       2518         Reef.18    #d8daeb      2
## Mcapitata_HTHC_TP4_1124       1124      Reef.42.43    #542788      2
## Mcapitata_HTHC_TP4_1328       1328      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP4_2204       2204      Reef.42.43    #542788      2
## Mcapitata_HTHC_TP5_1345       1345         Reef.18    #d8daeb      2
## Mcapitata_HTHC_TP5_1449       1449 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTHC_TP5_1694       1694      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP6_1164       1164      Reef.35.36    #998ec3      2
## Mcapitata_HTHC_TP6_1317       1317      Reef.35.36    #998ec3      2
## Mcapitata_HTHC_TP6_1604       1604      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP7_1126       1126      Reef.42.43    #542788      2
## Mcapitata_HTHC_TP7_1250       1250      Reef.42.43    #542788      2
## Mcapitata_HTHC_TP7_2419       2419         Reef.18    #d8daeb      2
## Mcapitata_HTHC_TP8_1082       1082      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP8_1246       1246      Reef.42.43    #542788      2
## Mcapitata_HTHC_TP8_2067       2067            HIMB    #b35806      2
## Mcapitata_HTHC_TP9_1078       1078      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP9_1331       1331      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP9_2009       2009 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTHC_TP10_1074      1074      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP10_1332      1332      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP10_1689      1689      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP11_1178      1178 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTHC_TP11_1270      1270            HIMB    #b35806      2
## Mcapitata_HTHC_TP11_2511      2511         Reef.18    #d8daeb      2
## Mcapitata_HTHC_TP12_1140      1140            HIMB    #b35806      2
## Mcapitata_HTHC_TP12_1274      1274            HIMB    #b35806      2
## Mcapitata_HTHC_TP12_2190      2190      Reef.42.43    #542788      2
## SRR5453739                Mcap_KA2          Kiholo    #33a02c      2
## SRR5453740                Mcap_KA3          Kiholo    #33a02c      2
## SRR5453741                Mcap_KA5          Kiholo    #33a02c      2
## SRR5453742                Mcap_KH1          Kiholo    #33a02c      2
## SRR5453743                Mcap_KH2          Kiholo    #33a02c      2
## SRR5453744                Mcap_KH8          Kiholo    #33a02c      2
## SRR5453745                Mcap_KU2          Kiholo    #33a02c      2
## SRR5453746                Mcap_KU3          Kiholo    #33a02c      2
## SRR5453747                Mcap_KU5          Kiholo    #33a02c      2
## SRR5453748                Mcap_WA8         Waiopae    #1f78b4      2
## SRR5453749                Mcap_WA1         Waiopae    #1f78b4      2
## SRR5453750                Mcap_WA2         Waiopae    #1f78b4      2
## SRR5453751               Mcap_WA12         Waiopae    #1f78b4      2
## SRR5453752                Mcap_WA3         Waiopae    #1f78b4      2
## SRR5453753                Mcap_WA5         Waiopae    #1f78b4      2
## SRR5453754                Mcap_WH3         Waiopae    #1f78b4      2
## SRR5453755                Mcap_WH8         Waiopae    #1f78b4      2
## SRR5453756                Mcap_WH1         Waiopae    #1f78b4      2
## SRR5453757                Mcap_WH5         Waiopae    #1f78b4      2
## SRR5453758                Mcap_WH2         Waiopae    #1f78b4      2
## SRR5453759               Mcap_WH12         Waiopae    #1f78b4      2
## SRR5453760                Mcap_WU8         Waiopae    #1f78b4      2
## SRR5453761                Mcap_WU1         Waiopae    #1f78b4      2
## SRR5453762                Mcap_WU2         Waiopae    #1f78b4      2
## SRR5453763               Mcap_WU12         Waiopae    #1f78b4      2
## SRR5453764                Mcap_WU3         Waiopae    #1f78b4      2
## SRR5453765                Mcap_WU5         Waiopae    #1f78b4      2
##                          ploidy_color       group group_color
## Mcapitata_ATAC_TP11_2302      #01665e      Group1     #1f78b4
## Mcapitata_HTAC_TP11_2380      #01665e      Group1     #1f78b4
## Mcapitata_ATAC_TP6_2402       #01665e      Group2     #33a02c
## Mcapitata_ATAC_TP12_2403      #01665e      Group2     #33a02c
## Mcapitata_ATAC_TP1_1037       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP1_1600       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP1_1652       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP3_1101       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP3_1548       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP3_1628       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP4_1108       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP4_1609       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP4_1651       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP5_1196       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP5_1610       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP5_1776       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP6_1114       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP6_1611       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP7_1058       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP7_1455       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP7_1499       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP8_1083       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP8_1436       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP8_1779       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP9_1121       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP9_1420       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP9_1580       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP10_1095      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP10_1561      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP10_1631      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP11_1076      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP11_1644      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP12_1120      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP12_1452      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP1_1218       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP1_1826       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP1_2068       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP3_1544       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP3_2731       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP3_2866       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP4_1221       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP4_2561       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP4_2734       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP5_1229       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP5_1706       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP5_2986       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP6_1212       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP6_2016       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP6_2555       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP7_1223       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP7_2860       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP7_2875       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP8_1260       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP8_2735       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP8_2753       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP9_1148       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP9_2862       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP9_2995       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP10_1204      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP10_2554      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP10_2737      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP11_1237      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP11_2188      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP11_2756      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP12_1154      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP12_2736      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP12_2990      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP1_1579       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP1_2153       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP1_2183       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP3_1289       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP3_1751       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP3_2021       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP4_1269       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP4_1481       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP4_2000       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP5_1321       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP5_1583       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP5_1997       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP6_1496       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP6_1588       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP6_1705       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP7_1278       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP7_1645       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP7_1722       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP8_1235       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP8_2386       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP8_2410       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP9_1306       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP9_1467       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP9_2412       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP10_1315      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP10_1478      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP10_1754      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP11_1248      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP11_1562      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP12_1632      #7570b3   Ungrouped     #808080
## Mcapitata_HTAC_TP12_1729      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP12_2007      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP1_1145       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP1_1323       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP1_2081       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP3_1128       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP3_1277       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP3_2518       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP4_1124       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP4_1328       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP4_2204       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP5_1345       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP5_1449       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP5_1694       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP6_1164       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP6_1317       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP6_1604       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP7_1126       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP7_1250       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP7_2419       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP8_1082       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP8_1246       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP8_2067       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP9_1078       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP9_1331       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP9_2009       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP10_1074      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP10_1332      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP10_1689      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP11_1178      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP11_1270      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP11_2511      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP12_1140      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP12_1274      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP12_2190      #01665e   Ungrouped     #808080
## SRR5453739                    #01665e  SRA-Kiholo     #543005
## SRR5453740                    #01665e  SRA-Kiholo     #543005
## SRR5453741                    #01665e  SRA-Kiholo     #543005
## SRR5453742                    #01665e  SRA-Kiholo     #543005
## SRR5453743                    #01665e  SRA-Kiholo     #543005
## SRR5453744                    #01665e  SRA-Kiholo     #543005
## SRR5453745                    #01665e  SRA-Kiholo     #543005
## SRR5453746                    #01665e  SRA-Kiholo     #543005
## SRR5453747                    #01665e  SRA-Kiholo     #543005
## SRR5453748                    #01665e SRA-Waiopae     #dfc27d
## SRR5453749                    #01665e SRA-Waiopae     #dfc27d
## SRR5453750                    #01665e SRA-Waiopae     #dfc27d
## SRR5453751                    #01665e SRA-Waiopae     #dfc27d
## SRR5453752                    #01665e SRA-Waiopae     #dfc27d
## SRR5453753                    #01665e SRA-Waiopae     #dfc27d
## SRR5453754                    #01665e SRA-Waiopae     #dfc27d
## SRR5453755                    #01665e SRA-Waiopae     #dfc27d
## SRR5453756                    #01665e SRA-Waiopae     #dfc27d
## SRR5453757                    #01665e SRA-Waiopae     #dfc27d
## SRR5453758                    #01665e SRA-Waiopae     #dfc27d
## SRR5453759                    #01665e SRA-Waiopae     #dfc27d
## SRR5453760                    #01665e SRA-Waiopae     #dfc27d
## SRR5453761                    #01665e SRA-Waiopae     #dfc27d
## SRR5453762                    #01665e SRA-Waiopae     #dfc27d
## SRR5453763                    #01665e SRA-Waiopae     #dfc27d
## SRR5453764                    #01665e SRA-Waiopae     #dfc27d
## SRR5453765                    #01665e SRA-Waiopae     #dfc27d
```





# Load similarity dendrogram

Load similarity dendrogram from `vcf_clone_detect.py`. Have to recompute the clusters because a dendro written and read from file does not seem to be the same as one generated directly from hclust (heatmap.2 wont work with a dendro loaded from file, only a fresh one generated in this session). 

```r
pairwise_percent_shared <- read.table("../vcf_clone_detect/GVCFall.filtered.recode.vcf.gz.allelic_similarity.full.tsv", sep='\t', header=T)
pairwise_percent_shared.matrix <- xtabs(match_perc ~ ind1 + ind2, data=pairwise_percent_shared)
pairwise_percent_shared.matrix.dist   <- dist(pairwise_percent_shared.matrix, method="euclidean")
pairwise_percent_shared.matrix.hclust <- hclust(pairwise_percent_shared.matrix.dist, method="complete")
simScore.dendro <- as.dendrogram(pairwise_percent_shared.matrix.hclust)
write.dendrogram(simScore.dendro, "cluster_dendrogram.tre", edges=TRUE)
plot(pairwise_percent_shared.matrix.hclust)
```

![](plot_vcftools_relatedness2_results_files/figure-html/load_simScore_dendrogram-1.png)<!-- -->





# Load pairwise relatedness scores (relScore) and cluster samples

Load pairwise relatedness scores and convert if from one line per pair format into a symetric matrix uisng the `xtabs` function.

```r
pairwise_relatedness <- read.table("GVCFall.filtered.recode.relatedness2", sep='\t', header=T)
pairwise_relatedness.matrix <- xtabs(RELATEDNESS_PHI ~ INDV1 + INDV2, data=pairwise_relatedness)
```



Process the pairwise relatedness scores data to get it ready for plotting.

```r
# Make a copy of input matrix
tmp.mat <- pairwise_relatedness.matrix

# Set negative relatedness values to 0
tmp.mat[tmp.mat < 0] <- 0

# Generate clusters from the relatedness scores matrix
tmp.mat.dist   <- dist(tmp.mat, method="euclidean")
tmp.mat.hclust <- hclust(tmp.mat.dist, method="complete")

# Get sample colors for sidebar in the same order as the matrix
Rcols <- samples.info[rownames(tmp.mat),]$group_color
Ccols <- samples.info[colnames(tmp.mat),]$group_color

# Get min and max values of matrix
tmp.mat.min <- min(tmp.mat); tmp.mat.min
```

```
## [1] 0
```

```r
tmp.mat.max <- max(tmp.mat); tmp.mat.max
```

```
## [1] 0.5
```



Write ordered matrix to file incase we want it later.

```r
tmp.mat.order <- dendro_data(simScore.dendro, type="rectangle")[["labels"]][["label"]]
tmp.mat.ordered <- tmp.mat[
  order(match(rownames(tmp.mat), rev(tmp.mat.order))),
  order(match(colnames(tmp.mat), tmp.mat.order))
  ]
write.table(as.matrix(tmp.mat.ordered), 
            "GVCFall.filtered.recode.relatedness2.matrix.tsv", 
            sep='\t', quote=FALSE, 
            row.names=TRUE, col.names=TRUE)
```





# Plot pairwise relatedness scores

Plot relScores heatmap.

```r
par(cex.main=0.7, cex.lab=0.9, cex.axis=0.9) # Change the font size of the legend title,labels and axis
heatmap.2(x=tmp.mat,
          Rowv = simScore.dendro,
          Colv = simScore.dendro,
          col=cividis(100),
          trace="none", # Dont draw "trace" line
          key=TRUE, # Show color-key
          margins=c(8,8), # Numeric vector of length 2 containing the margins for column and row names, respectively.
          ColSideColors=Ccols,
          RowSideColors=Rcols,
          offsetRow=0,
          offsetCol=0,
          cexRow=cexSize, # Positive numbers, used as cex.axis in for the row or column axis labeling.
          cexCol=cexSize,
          cellnote=formatC(round(tmp.mat, 3), big.mark=","),
          notecex=0.1,
          key.xtickfun=function() {
                 breaks <- parent.frame()$breaks
                 return(list(
                      at=parent.frame()$scale01(seq(0, 0.5, 0.05)),
                      labels=      as.character(seq(0, 0.5, 0.05))
                      ))
          },
          )
```

![](plot_vcftools_relatedness2_results_files/figure-html/plot_relScore_heatmap-1.png)<!-- -->



Plot relScores heatmap without cell notes/values.

```r
par(cex.main=0.7, cex.lab=0.9, cex.axis=0.9) # Change the font size of the legend title,labels and axis
heatmap.2(x=tmp.mat,
          Rowv = simScore.dendro,
          Colv = simScore.dendro,
          col=cividis(100),
          trace="none", # Don't draw "trace" line
          key=TRUE, # Show color-key
          margins=c(8,8), # Numeric vector of length 2 containing the margins for column and row names, respectively.
          ColSideColors=Ccols,
          RowSideColors=Rcols,
          offsetRow=0,
          offsetCol=0,
          cexRow=cexSize, # Positive numbers, used as cex.axis in for the row or column axis labeling.
          cexCol=cexSize,
          notecex=0.1,
          key.xtickfun=function() {
                 breaks <- parent.frame()$breaks
                 return(list(
                      at=parent.frame()$scale01(seq(0, 0.5, 0.05)),
                      labels=      as.character(seq(0, 0.5, 0.05))
                      ))
          },
          )
```

![](plot_vcftools_relatedness2_results_files/figure-html/plot_relScore_heatmap_noCellNotes-1.png)<!-- -->



Plot relScores heatmap with specific colors applied to selected ranges of values (i.e., value ranges picked from the histogram on the top left).

```r
par(cex.main=0.7, cex.lab=0.9, cex.axis=0.9) # Change the font size of the legend title,labels and axis
heatmap.2(x=tmp.mat,
          Rowv = simScore.dendro,
          Colv = simScore.dendro,
          col=c(rep(cividis(5)[1], length(seq(0.01, 0.18, 0.01))),
                rep(cividis(5)[2], length(seq(0.19, 0.26, 0.01))),
                rep(cividis(5)[3], length(seq(0.27, 0.31, 0.01))),
                rep(cividis(5)[4], length(seq(0.32, 0.42, 0.01))),
                rep(cividis(5)[5], length(seq(0.43, 0.50, 0.01)))
                ),
          breaks=seq(0.00, 0.50, 0.01),
          trace="none", # Dont draw "trace" line
          key=TRUE, # Show color-key
          margins=c(8,8), # Numeric vector of length 2 containing the margins for column and row names, respectively.
          ColSideColors=Ccols,
          RowSideColors=Rcols,
          offsetRow=0,
          offsetCol=0,
          cexRow=cexSize, # Positive numbers, used as cex.axis in for the row or column axis labeling.
          cexCol=cexSize,
          cellnote=formatC(round(tmp.mat, 3), big.mark=","),
          notecex=0.1,
          key.xtickfun=function() {
                 breaks <- parent.frame()$breaks
                 return(list(
                      at=parent.frame()$scale01(seq(0, 0.5, 0.05)),
                      labels=      as.character(seq(0, 0.5, 0.05))
                      ))
          },
          )
```

![](plot_vcftools_relatedness2_results_files/figure-html/plot_relScore_heatmap_withCutoff-1.png)<!-- -->



Plot relScores heatmap with specific colors applied to selected ranges of values (i.e., value ranges picked from the histogram on the top left) without cell notes/values.

```r
par(cex.main=0.7, cex.lab=0.9, cex.axis=0.9) # Change the font size of the legend title,labels and axis
heatmap.2(x=tmp.mat,
          Rowv = simScore.dendro,
          Colv = simScore.dendro,
          col=c(rep(cividis(5)[1], length(seq(0.01, 0.18, 0.01))),
                rep(cividis(5)[2], length(seq(0.19, 0.26, 0.01))),
                rep(cividis(5)[3], length(seq(0.27, 0.31, 0.01))),
                rep(cividis(5)[4], length(seq(0.32, 0.42, 0.01))),
                rep(cividis(5)[5], length(seq(0.43, 0.50, 0.01)))
                ),
          breaks=seq(0.00, 0.50, 0.01),
          trace="none", # Dont draw "trace" line
          key=TRUE, # Show color-key
          margins=c(8,8), # Numeric vector of length 2 containing the margins for column and row names, respectively.
          ColSideColors=Ccols,
          RowSideColors=Rcols,
          offsetRow=0,
          offsetCol=0,
          cexRow=cexSize, # Positive numbers, used as cex.axis in for the row or column axis labeling.
          cexCol=cexSize,
          notecex=0.1,
          key.xtickfun=function() {
                 breaks <- parent.frame()$breaks
                 return(list(
                      at=parent.frame()$scale01(seq(0, 0.5, 0.05)),
                      labels=      as.character(seq(0, 0.5, 0.05))
                      ))
          },
          )
```

![](plot_vcftools_relatedness2_results_files/figure-html/plot_relScore_heatmap_withCutoff_noCellNotes-1.png)<!-- -->





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
##  [1] phylogram_2.1.0    plotly_4.10.0      knitr_1.40         viridis_0.6.2     
##  [5] viridisLite_0.4.1  RColorBrewer_1.1-3 cowplot_1.1.1      ggdendro_0.1.23   
##  [9] tibble_3.1.8       RcppCNPy_0.2.11    reshape2_1.4.4     gplots_3.1.3      
## [13] ggplot2_3.3.6      dplyr_1.0.10       readxl_1.4.1      
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_1.0.9         lattice_0.20-45    ape_5.6-2          tidyr_1.2.1       
##  [5] gtools_3.9.3       assertthat_0.2.1   digest_0.6.29      utf8_1.2.2        
##  [9] R6_2.5.1           cellranger_1.1.0   plyr_1.8.7         evaluate_0.16     
## [13] highr_0.9          httr_1.4.4         pillar_1.8.1       rlang_1.0.5       
## [17] lazyeval_0.2.2     rstudioapi_0.14    data.table_1.14.2  jquerylib_0.1.4   
## [21] rmarkdown_2.16     stringr_1.4.1      htmlwidgets_1.5.4  munsell_0.5.0     
## [25] compiler_4.1.2     xfun_0.33          pkgconfig_2.0.3    htmltools_0.5.3   
## [29] tidyselect_1.1.2   gridExtra_2.3      fansi_1.0.3        withr_2.5.0       
## [33] MASS_7.3-58.1      bitops_1.0-7       grid_4.1.2         nlme_3.1-159      
## [37] jsonlite_1.8.0     gtable_0.3.1       lifecycle_1.0.2    DBI_1.1.3         
## [41] magrittr_2.0.3     scales_1.2.1       KernSmooth_2.23-20 cli_3.4.0         
## [45] stringi_1.7.8      cachem_1.0.6       bslib_0.4.0        generics_0.1.3    
## [49] vctrs_0.4.1        tools_4.1.2        glue_1.6.2         purrr_0.3.4       
## [53] parallel_4.1.2     fastmap_1.1.0      yaml_2.3.5         colorspace_2.0-3  
## [57] caTools_1.18.2     sass_0.4.2
```
