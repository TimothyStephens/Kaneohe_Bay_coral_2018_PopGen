---
title: "Plot `ANGSD` results for *M. capitata* RNA-seq samples from this study and SRA"
author: "Timothy Stephens"
date: "22/09/2022"
output: 
  html_document:
    keep_md: yes
---



## Setup

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
library(knitr)
library(phylogram)
options(scipen = 999) #Prevent scientific notation
```





# Load Metadata 

Load file with annotation for each sample and order it using the order of the samples in the results files. 

```r
samples.info <- read.table("../../../samples_Mcapitata_ALL.annotations.txt", header=T, comment.char='')
rownames(samples.info) <- samples.info$sample

samples.order <- read.table("bam.filelist.labels", header=F)
samples.info <- samples.info[samples.order$V1,]
samples.info
```

```
##                                            sample   species
## Mcapitata_ATAC_TP10_1095 Mcapitata_ATAC_TP10_1095 Mcapitata
## Mcapitata_ATAC_TP10_1561 Mcapitata_ATAC_TP10_1561 Mcapitata
## Mcapitata_ATAC_TP10_1631 Mcapitata_ATAC_TP10_1631 Mcapitata
## Mcapitata_ATAC_TP1_1037   Mcapitata_ATAC_TP1_1037 Mcapitata
## Mcapitata_ATAC_TP11_1076 Mcapitata_ATAC_TP11_1076 Mcapitata
## Mcapitata_ATAC_TP11_1644 Mcapitata_ATAC_TP11_1644 Mcapitata
## Mcapitata_ATAC_TP11_2302 Mcapitata_ATAC_TP11_2302 Mcapitata
## Mcapitata_ATAC_TP1_1600   Mcapitata_ATAC_TP1_1600 Mcapitata
## Mcapitata_ATAC_TP1_1652   Mcapitata_ATAC_TP1_1652 Mcapitata
## Mcapitata_ATAC_TP12_1120 Mcapitata_ATAC_TP12_1120 Mcapitata
## Mcapitata_ATAC_TP12_1452 Mcapitata_ATAC_TP12_1452 Mcapitata
## Mcapitata_ATAC_TP12_2403 Mcapitata_ATAC_TP12_2403 Mcapitata
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
## Mcapitata_ATAC_TP6_2402   Mcapitata_ATAC_TP6_2402 Mcapitata
## Mcapitata_ATAC_TP7_1058   Mcapitata_ATAC_TP7_1058 Mcapitata
## Mcapitata_ATAC_TP7_1455   Mcapitata_ATAC_TP7_1455 Mcapitata
## Mcapitata_ATAC_TP7_1499   Mcapitata_ATAC_TP7_1499 Mcapitata
## Mcapitata_ATAC_TP8_1083   Mcapitata_ATAC_TP8_1083 Mcapitata
## Mcapitata_ATAC_TP8_1436   Mcapitata_ATAC_TP8_1436 Mcapitata
## Mcapitata_ATAC_TP8_1779   Mcapitata_ATAC_TP8_1779 Mcapitata
## Mcapitata_ATAC_TP9_1121   Mcapitata_ATAC_TP9_1121 Mcapitata
## Mcapitata_ATAC_TP9_1420   Mcapitata_ATAC_TP9_1420 Mcapitata
## Mcapitata_ATAC_TP9_1580   Mcapitata_ATAC_TP9_1580 Mcapitata
## Mcapitata_ATHC_TP10_1204 Mcapitata_ATHC_TP10_1204 Mcapitata
## Mcapitata_ATHC_TP10_2554 Mcapitata_ATHC_TP10_2554 Mcapitata
## Mcapitata_ATHC_TP10_2737 Mcapitata_ATHC_TP10_2737 Mcapitata
## Mcapitata_ATHC_TP11_1237 Mcapitata_ATHC_TP11_1237 Mcapitata
## Mcapitata_ATHC_TP11_2188 Mcapitata_ATHC_TP11_2188 Mcapitata
## Mcapitata_ATHC_TP1_1218   Mcapitata_ATHC_TP1_1218 Mcapitata
## Mcapitata_ATHC_TP11_2756 Mcapitata_ATHC_TP11_2756 Mcapitata
## Mcapitata_ATHC_TP1_1826   Mcapitata_ATHC_TP1_1826 Mcapitata
## Mcapitata_ATHC_TP1_2068   Mcapitata_ATHC_TP1_2068 Mcapitata
## Mcapitata_ATHC_TP12_1154 Mcapitata_ATHC_TP12_1154 Mcapitata
## Mcapitata_ATHC_TP12_2736 Mcapitata_ATHC_TP12_2736 Mcapitata
## Mcapitata_ATHC_TP12_2990 Mcapitata_ATHC_TP12_2990 Mcapitata
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
## Mcapitata_HTAC_TP10_1315 Mcapitata_HTAC_TP10_1315 Mcapitata
## Mcapitata_HTAC_TP10_1478 Mcapitata_HTAC_TP10_1478 Mcapitata
## Mcapitata_HTAC_TP10_1754 Mcapitata_HTAC_TP10_1754 Mcapitata
## Mcapitata_HTAC_TP11_1248 Mcapitata_HTAC_TP11_1248 Mcapitata
## Mcapitata_HTAC_TP11_1562 Mcapitata_HTAC_TP11_1562 Mcapitata
## Mcapitata_HTAC_TP11_2380 Mcapitata_HTAC_TP11_2380 Mcapitata
## Mcapitata_HTAC_TP1_1579   Mcapitata_HTAC_TP1_1579 Mcapitata
## Mcapitata_HTAC_TP1_2153   Mcapitata_HTAC_TP1_2153 Mcapitata
## Mcapitata_HTAC_TP12_1632 Mcapitata_HTAC_TP12_1632 Mcapitata
## Mcapitata_HTAC_TP12_1729 Mcapitata_HTAC_TP12_1729 Mcapitata
## Mcapitata_HTAC_TP1_2183   Mcapitata_HTAC_TP1_2183 Mcapitata
## Mcapitata_HTAC_TP12_2007 Mcapitata_HTAC_TP12_2007 Mcapitata
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
## Mcapitata_HTHC_TP10_1074 Mcapitata_HTHC_TP10_1074 Mcapitata
## Mcapitata_HTHC_TP10_1332 Mcapitata_HTHC_TP10_1332 Mcapitata
## Mcapitata_HTHC_TP10_1689 Mcapitata_HTHC_TP10_1689 Mcapitata
## Mcapitata_HTHC_TP11_1178 Mcapitata_HTHC_TP11_1178 Mcapitata
## Mcapitata_HTHC_TP11_1270 Mcapitata_HTHC_TP11_1270 Mcapitata
## Mcapitata_HTHC_TP1_1145   Mcapitata_HTHC_TP1_1145 Mcapitata
## Mcapitata_HTHC_TP11_2511 Mcapitata_HTHC_TP11_2511 Mcapitata
## Mcapitata_HTHC_TP1_1323   Mcapitata_HTHC_TP1_1323 Mcapitata
## Mcapitata_HTHC_TP1_2081   Mcapitata_HTHC_TP1_2081 Mcapitata
## Mcapitata_HTHC_TP12_1140 Mcapitata_HTHC_TP12_1140 Mcapitata
## Mcapitata_HTHC_TP12_1274 Mcapitata_HTHC_TP12_1274 Mcapitata
## Mcapitata_HTHC_TP12_2190 Mcapitata_HTHC_TP12_2190 Mcapitata
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
## Mcapitata_ATAC_TP10_1095                                        ATAC      TP10
## Mcapitata_ATAC_TP10_1561                                        ATAC      TP10
## Mcapitata_ATAC_TP10_1631                                        ATAC      TP10
## Mcapitata_ATAC_TP1_1037                                         ATAC       TP1
## Mcapitata_ATAC_TP11_1076                                        ATAC      TP11
## Mcapitata_ATAC_TP11_1644                                        ATAC      TP11
## Mcapitata_ATAC_TP11_2302                                        ATAC      TP11
## Mcapitata_ATAC_TP1_1600                                         ATAC       TP1
## Mcapitata_ATAC_TP1_1652                                         ATAC       TP1
## Mcapitata_ATAC_TP12_1120                                        ATAC      TP12
## Mcapitata_ATAC_TP12_1452                                        ATAC      TP12
## Mcapitata_ATAC_TP12_2403                                        ATAC      TP12
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
## Mcapitata_ATAC_TP6_2402                                         ATAC       TP6
## Mcapitata_ATAC_TP7_1058                                         ATAC       TP7
## Mcapitata_ATAC_TP7_1455                                         ATAC       TP7
## Mcapitata_ATAC_TP7_1499                                         ATAC       TP7
## Mcapitata_ATAC_TP8_1083                                         ATAC       TP8
## Mcapitata_ATAC_TP8_1436                                         ATAC       TP8
## Mcapitata_ATAC_TP8_1779                                         ATAC       TP8
## Mcapitata_ATAC_TP9_1121                                         ATAC       TP9
## Mcapitata_ATAC_TP9_1420                                         ATAC       TP9
## Mcapitata_ATAC_TP9_1580                                         ATAC       TP9
## Mcapitata_ATHC_TP10_1204                                        ATHC      TP10
## Mcapitata_ATHC_TP10_2554                                        ATHC      TP10
## Mcapitata_ATHC_TP10_2737                                        ATHC      TP10
## Mcapitata_ATHC_TP11_1237                                        ATHC      TP11
## Mcapitata_ATHC_TP11_2188                                        ATHC      TP11
## Mcapitata_ATHC_TP1_1218                                         ATHC       TP1
## Mcapitata_ATHC_TP11_2756                                        ATHC      TP11
## Mcapitata_ATHC_TP1_1826                                         ATHC       TP1
## Mcapitata_ATHC_TP1_2068                                         ATHC       TP1
## Mcapitata_ATHC_TP12_1154                                        ATHC      TP12
## Mcapitata_ATHC_TP12_2736                                        ATHC      TP12
## Mcapitata_ATHC_TP12_2990                                        ATHC      TP12
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
## Mcapitata_HTAC_TP10_1315                                        HTAC      TP10
## Mcapitata_HTAC_TP10_1478                                        HTAC      TP10
## Mcapitata_HTAC_TP10_1754                                        HTAC      TP10
## Mcapitata_HTAC_TP11_1248                                        HTAC      TP11
## Mcapitata_HTAC_TP11_1562                                        HTAC      TP11
## Mcapitata_HTAC_TP11_2380                                        HTAC      TP11
## Mcapitata_HTAC_TP1_1579                                         HTAC       TP1
## Mcapitata_HTAC_TP1_2153                                         HTAC       TP1
## Mcapitata_HTAC_TP12_1632                                        HTAC      TP12
## Mcapitata_HTAC_TP12_1729                                        HTAC      TP12
## Mcapitata_HTAC_TP1_2183                                         HTAC       TP1
## Mcapitata_HTAC_TP12_2007                                        HTAC      TP12
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
## Mcapitata_HTHC_TP10_1074                                        HTHC      TP10
## Mcapitata_HTHC_TP10_1332                                        HTHC      TP10
## Mcapitata_HTHC_TP10_1689                                        HTHC      TP10
## Mcapitata_HTHC_TP11_1178                                        HTHC      TP11
## Mcapitata_HTHC_TP11_1270                                        HTHC      TP11
## Mcapitata_HTHC_TP1_1145                                         HTHC       TP1
## Mcapitata_HTHC_TP11_2511                                        HTHC      TP11
## Mcapitata_HTHC_TP1_1323                                         HTHC       TP1
## Mcapitata_HTHC_TP1_2081                                         HTHC       TP1
## Mcapitata_HTHC_TP12_1140                                        HTHC      TP12
## Mcapitata_HTHC_TP12_1274                                        HTHC      TP12
## Mcapitata_HTHC_TP12_2190                                        HTHC      TP12
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
## Mcapitata_ATAC_TP10_1095      1095      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP10_1561      1561      Reef.11.13    #fee0b6      2
## Mcapitata_ATAC_TP10_1631      1631            HIMB    #b35806      2
## Mcapitata_ATAC_TP1_1037       1037      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP11_1076      1076      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP11_1644      1644            HIMB    #b35806      2
## Mcapitata_ATAC_TP11_2302      2302         Reef.18    #d8daeb      2
## Mcapitata_ATAC_TP1_1600       1600      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP1_1652       1652            HIMB    #b35806      2
## Mcapitata_ATAC_TP12_1120      1120      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP12_1452      1452      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP12_2403      2403         Reef.18    #d8daeb      2
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
## Mcapitata_ATAC_TP6_2402       2402         Reef.18    #d8daeb      2
## Mcapitata_ATAC_TP7_1058       1058 Lilipuna.Fringe    #f1a340      2
## Mcapitata_ATAC_TP7_1455       1455      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP7_1499       1499      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP8_1083       1083      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP8_1436       1436      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP8_1779       1779      Reef.42.43    #542788      2
## Mcapitata_ATAC_TP9_1121       1121            HIMB    #b35806      2
## Mcapitata_ATAC_TP9_1420       1420      Reef.35.36    #998ec3      2
## Mcapitata_ATAC_TP9_1580       1580      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP10_1204      1204      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP10_2554      2554         Reef.18    #d8daeb      2
## Mcapitata_ATHC_TP10_2737      2737      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP11_1237      1237      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP11_2188      2188      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP1_1218       1218      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP11_2756      2756      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP1_1826       1826      Reef.11.13    #fee0b6      2
## Mcapitata_ATHC_TP1_2068       2068      Reef.35.36    #998ec3      2
## Mcapitata_ATHC_TP12_1154      1154            HIMB    #b35806      2
## Mcapitata_ATHC_TP12_2736      2736      Reef.42.43    #542788      2
## Mcapitata_ATHC_TP12_2990      2990         Reef.18    #d8daeb      2
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
## Mcapitata_HTAC_TP10_1315      1315 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTAC_TP10_1478      1478      Reef.11.13    #fee0b6      2
## Mcapitata_HTAC_TP10_1754      1754            HIMB    #b35806      2
## Mcapitata_HTAC_TP11_1248      1248            HIMB    #b35806      2
## Mcapitata_HTAC_TP11_1562      1562      Reef.11.13    #fee0b6      2
## Mcapitata_HTAC_TP11_2380      2380         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP1_1579       1579 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTAC_TP1_2153       2153         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP12_1632      1632            HIMB    #b35806      4
## Mcapitata_HTAC_TP12_1729      1729         Reef.18    #d8daeb      2
## Mcapitata_HTAC_TP1_2183       2183      Reef.42.43    #542788      2
## Mcapitata_HTAC_TP12_2007      2007      Reef.42.43    #542788      2
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
## Mcapitata_HTHC_TP10_1074      1074      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP10_1332      1332      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP10_1689      1689      Reef.11.13    #fee0b6      2
## Mcapitata_HTHC_TP11_1178      1178 Lilipuna.Fringe    #f1a340      2
## Mcapitata_HTHC_TP11_1270      1270            HIMB    #b35806      2
## Mcapitata_HTHC_TP1_1145       1145            HIMB    #b35806      2
## Mcapitata_HTHC_TP11_2511      2511         Reef.18    #d8daeb      2
## Mcapitata_HTHC_TP1_1323       1323      Reef.35.36    #998ec3      2
## Mcapitata_HTHC_TP1_2081       2081            HIMB    #b35806      2
## Mcapitata_HTHC_TP12_1140      1140            HIMB    #b35806      2
## Mcapitata_HTHC_TP12_1274      1274            HIMB    #b35806      2
## Mcapitata_HTHC_TP12_2190      2190      Reef.42.43    #542788      2
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
## Mcapitata_ATAC_TP10_1095      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP10_1561      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP10_1631      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP1_1037       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP11_1076      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP11_1644      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP11_2302      #01665e      Group1     #1f78b4
## Mcapitata_ATAC_TP1_1600       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP1_1652       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP12_1120      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP12_1452      #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP12_2403      #01665e      Group2     #33a02c
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
## Mcapitata_ATAC_TP6_2402       #01665e      Group2     #33a02c
## Mcapitata_ATAC_TP7_1058       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP7_1455       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP7_1499       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP8_1083       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP8_1436       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP8_1779       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP9_1121       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP9_1420       #01665e   Ungrouped     #808080
## Mcapitata_ATAC_TP9_1580       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP10_1204      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP10_2554      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP10_2737      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP11_1237      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP11_2188      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP1_1218       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP11_2756      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP1_1826       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP1_2068       #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP12_1154      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP12_2736      #01665e   Ungrouped     #808080
## Mcapitata_ATHC_TP12_2990      #01665e   Ungrouped     #808080
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
## Mcapitata_HTAC_TP10_1315      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP10_1478      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP10_1754      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP11_1248      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP11_1562      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP11_2380      #01665e      Group1     #1f78b4
## Mcapitata_HTAC_TP1_1579       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP1_2153       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP12_1632      #7570b3   Ungrouped     #808080
## Mcapitata_HTAC_TP12_1729      #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP1_2183       #01665e   Ungrouped     #808080
## Mcapitata_HTAC_TP12_2007      #01665e   Ungrouped     #808080
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
## Mcapitata_HTHC_TP10_1074      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP10_1332      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP10_1689      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP11_1178      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP11_1270      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP1_1145       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP11_2511      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP1_1323       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP1_2081       #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP12_1140      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP12_1274      #01665e   Ungrouped     #808080
## Mcapitata_HTHC_TP12_2190      #01665e   Ungrouped     #808080
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

![](plot_ANGSD_results_files/figure-html/load_simScore_dendrogram-1.png)<!-- -->





# Plot Admixture based on two PC

A function that wraps up all the functions required to parse and plot the input data.
Will use the supplied hclust object to create dendrogram and order bars or will perform clustering itself.

```r
cluster_plot <- function(data, samples, selected_columns, 
                         plot.dendro=NULL, colors=NULL,
                         dist.method="euclidean", hclust.method="complete",
                         x.label.font.size=5, legend.size.scale=0.5,
                         legend.title="Groups",
                         rel_heights=c(0.3, 1.0, 0.3)) {
  # Code modified from: https://stackoverflow.com/questions/44646488/stacked-bar-plot-with-hierarchical-clustering-dendrogram
  
  #########################
  #### Cluster samples ####
  #########################
  if(is.null(plot.dendro)) {
    d <- dist(as.matrix(data), method=dist.method)
    hc <- hclust(d, method=hclust.method)
    plot.dendro <- as.dendrogram(hc)
  }
  
  ##############
  #### Plot ####
  ##############
  # Convert dendrogram to segment data
  dend_data <- dendro_data(plot.dendro, type="rectangle")
  segment_data <- dend_data[["segments"]]
  dend_labels <- gsub('\'', '', dend_data$labels$label)

  # Sample positions df
  sample_pos_table <- with(dend_data$labels, 
                           data.frame(x_center = x, sample = as.character(label), width = 0.9))
  # X-axis limits
  axis_limits <- with(sample_pos_table, 
                      c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))) + 
                      0.1 * c(-1, 1) # extra spacing: 0.1
  
  # Plot clustering dendrogram
  p.dendr <- ggplot(segment_data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_y_continuous(expand = c(0, 0, 0, 0.1)) + 
    scale_x_continuous(breaks = sample_pos_table$x_center, 
                       labels = NULL, 
                       limits = axis_limits, 
                       expand = c(0, 0)) + 
    labs(x = "", y = "Distance", colour = "", size = "") +
    theme_bw() + 
    theme(panel.grid.minor=element_blank(),
          panel.grid.major=element_blank(),
          plot.margin=unit(c(0,0,0,0), "mm"))
  
  # Plot bar chart
  p.barplot <- data %>%
    rownames_to_column() %>%
    arrange(factor(rowname, levels = dend_labels), dend_labels) %>%
    column_to_rownames("rowname") %>%
    as.matrix() %>%
    melt() %>%
    rename(Sample=Var1, ITS2type=Var2, Abundance=value) %>%
    # Plot bar chart
    ggplot(aes(x=Sample, y=Abundance, fill=ITS2type)) +
      geom_bar(position="stack", stat="identity") +
      scale_y_continuous(expand = c(0, 0)) +
      scale_fill_manual(values = colorRampPalette(brewer.pal(12,"Paired"))(12)) +
      labs(x = "", y = "Frequency", fill=legend.title) +
      theme_bw() +
      theme(axis.text.x=element_blank(),
            legend.key.height= unit(2*legend.size.scale, 'mm'),
            legend.key.width= unit(2*legend.size.scale, 'mm'),
            legend.title = element_text(size=4*legend.size.scale),
            legend.text = element_text(size=4*legend.size.scale),
            plot.margin=unit(c(0,0,0,0), "mm")) +
      guides(fill=guide_legend(ncol=1))
  
  # Plot Meta Information about samples
  if(is.null(colors)) {
    colors <- rep(brewer.pal(brewer.pal.info["Set3", "maxcolors"], "Set3"), ceiling(nrow(samples)/12))
  }
  p.sampleInfo <- samples %>%
    select(all_of(selected_columns)) %>%
    rownames_to_column() %>% 
    arrange(factor(rowname, levels = dend_labels), dend_labels) %>% 
    column_to_rownames("rowname") %>% 
    as.matrix() %>% 
    melt() %>% 
    rename(Sample=Var1, MetaInfo=Var2, Value=value) %>% 
    mutate(Count = 1) %>% 
    filter(MetaInfo!="PlugID") %>%
    ggplot(aes(x=Sample, y=Count, fill=Value)) + 
      geom_bar(position="stack", stat="identity") + 
      scale_y_continuous(expand = c(0, 0)) + 
      scale_fill_manual(values=colors) +
      labs(x = "Samples", y = "Meta Information", fill="") +
      theme_bw() +
      theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1, size=x.label.font.size),
            axis.ticks.y=element_blank(),
            axis.text.y=element_blank(),
            legend.key.height= unit(2*legend.size.scale, 'mm'),
            legend.key.width= unit(2*legend.size.scale, 'mm'),
            legend.title = element_text(size=4*legend.size.scale),
            legend.text = element_text(size=4*legend.size.scale),
            plot.margin=unit(c(0,0,0,0), "mm")) +
      guides(fill=guide_legend(ncol=1))
  
  # Combine plots and print to screen
  comb <- plot_grid(p.dendr, p.barplot, p.sampleInfo, align = 'v', ncol = 1, axis = "lr", rel_heights = rel_heights)
  print(comb)
}
```


Load Admixture results into `R`.

```r
q  <- read.table("PCAngsd.angsd.beagle.gz.Admixture.admix.2.Q")
rownames(q) <- samples.order$V1
```



Setup the color pallet that we want to use for the sample metadata.

```r
t.reef <- samples.info %>% select(all_of(c("reef", "reef_color"))) %>% unique()
colnames(t.reef) <- c("label", "color")
t.ploidy<- samples.info %>% select(all_of(c("ploidy", "ploidy_color"))) %>% unique()
colnames(t.ploidy) <- c("label", "color")
t.group <- samples.info %>% select(all_of(c("group", "group_color"))) %>% unique()
colnames(t.group) <- c("label", "color")
t <- rbind(t.reef, t.ploidy, t.group)
label.colors <- as.vector(t[,2])
names(label.colors) <- as.vector(t[,1])
label.colors
```

```
##      Reef.35.36      Reef.11.13            HIMB         Reef.18      Reef.42.43 
##       "#998ec3"       "#fee0b6"       "#b35806"       "#d8daeb"       "#542788" 
## Lilipuna.Fringe          Kiholo         Waiopae               2               4 
##       "#f1a340"       "#33a02c"       "#1f78b4"       "#01665e"       "#7570b3" 
##       Ungrouped          Group1          Group2      SRA-Kiholo     SRA-Waiopae 
##       "#808080"       "#1f78b4"       "#33a02c"       "#543005"       "#dfc27d"
```



Plot Admixture results. Use the admixture matrix to generate the degdrogram and sample order.

```r
cluster_plot(q, samples.info, c("reef", "ploidy", "group"), colors=label.colors, x.label.font.size=2)
```

![](plot_ANGSD_results_files/figure-html/plot_Admixture-1.png)<!-- -->



Plot Admixture results using a supplied dendrogram object built from the pairwise similairty scores produced by `vcf_clone_detect.py`.

```r
cluster_plot(q, samples.info, c("reef", "ploidy", "group"), plot.dendro=simScore.dendro, colors=label.colors, x.label.font.size=2)
```

![](plot_ANGSD_results_files/figure-html/plot_Admixture_withPSSclust-1.png)<!-- -->





# Estimating Individual Allele Frequencies

Load PCAngsd results with estimated individual allele frequencies.

```r
C <- as.matrix(read.table("PCAngsd.angsd.beagle.gz.IndAlleleFreq.cov"))
e <- eigen(C)
pop <- samples.info
```



Plot PC1 and PC2 of PCAngsd results with estimated individual allele frequencies.

```r
A <- 1
B <- 2

plot(e$vectors[,c(A,B)], col=pop$group_color,  pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Group")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC2-1.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$reef_color,   pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Reef")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC2-2.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$ploidy_color, pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Ploidy")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC2-3.png)<!-- -->



Plot PC1 and PC2 of PCAngsd results with estimated individual allele frequencies using the `plotly` package.

```r
w<-800
h<-800

df <- as.data.frame(e$vectors)
df <- cbind(df, samples.info)

# Group
p <- plot_ly(df, x=~V1, y=~V2, text=rownames(df), type = "scatter",
             mode="markers", color=~group, colors=label.colors, marker=list(size=11),
             width=w, height=h)
p <- layout(p,title="Without individual allele frequency - colored by Group",
            xaxis=list(title="PC1"),
            yaxis=list(title="PC2"))
p
```

```{=html}
<div id="htmlwidget-cdcd891d93cde0f526e8" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-cdcd891d93cde0f526e8">{"x":{"visdat":{"148844be0bd0":["function () ","plotlyVisDat"]},"cur_data":"148844be0bd0","attrs":{"148844be0bd0":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009","SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747","SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#33a02c","#1f78b4","#01665e","#7570b3","#808080","#1f78b4","#33a02c","#543005","#dfc27d"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Group","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0301690858687305,-0.0308583956418033],"y":[0.0384051391045254,0.0379059230536485],"text":["Mcapitata_ATAC_TP11_2302","Mcapitata_HTAC_TP11_2380"],"mode":"markers","marker":{"color":"rgba(31,120,180,1)","size":11,"line":{"color":"rgba(31,120,180,1)"}},"type":"scatter","name":"Group1","textfont":{"color":"rgba(31,120,180,1)"},"error_y":{"color":"rgba(31,120,180,1)"},"error_x":{"color":"rgba(31,120,180,1)"},"line":{"color":"rgba(31,120,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0444844359594012,-0.0432579249161303],"y":[-0.0624821710739555,-0.0635391652669186],"text":["Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP6_2402"],"mode":"markers","marker":{"color":"rgba(51,160,44,1)","size":11,"line":{"color":"rgba(51,160,44,1)"}},"type":"scatter","name":"Group2","textfont":{"color":"rgba(51,160,44,1)"},"error_y":{"color":"rgba(51,160,44,1)"},"error_x":{"color":"rgba(51,160,44,1)"},"line":{"color":"rgba(51,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.145997445819095,0.174955367864438,0.174902797925439,0.151102413283901,0.159919495340554,0.157942533764846,0.146420407313794,0.173851484550974,0.176377286930581],"y":[0.192883367283657,0.158032196681908,0.158940044226475,0.121758939195166,0.124271170183639,0.12000084654977,0.191061022008456,0.152828046521998,0.158656760250303],"text":["SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747"],"mode":"markers","marker":{"color":"rgba(84,48,5,1)","size":11,"line":{"color":"rgba(84,48,5,1)"}},"type":"scatter","name":"SRA-Kiholo","textfont":{"color":"rgba(84,48,5,1)"},"error_y":{"color":"rgba(84,48,5,1)"},"error_x":{"color":"rgba(84,48,5,1)"},"line":{"color":"rgba(84,48,5,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.192963783550781,0.187604556302613,0.195375556575114,0.173015234268865,0.174672358864291,0.184831419673122,0.169261625046705,0.161388905503109,0.160889770031219,0.15102065362944,0.178247124058197,0.179315378194888,0.194924584218308,0.188394179356329,0.195885141133129,0.142341074734009,0.155785419073372,0.154281672744728],"y":[-0.285872345487089,-0.169758811039003,-0.365188395712793,0.0809955243749636,0.122361962615188,-0.1038389960395,0.0525769728321265,0.0742335665983742,0.0191980288362958,0.0754682482966374,0.0819136645421538,0.0788716541902661,-0.343242651831009,-0.174828341971528,-0.368863867281678,0.105571228632626,0.146926092689027,-0.131702627315251],"text":["SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"color":"rgba(223,194,125,1)","size":11,"line":{"color":"rgba(223,194,125,1)"}},"type":"scatter","name":"SRA-Waiopae","textfont":{"color":"rgba(223,194,125,1)"},"error_y":{"color":"rgba(223,194,125,1)"},"error_x":{"color":"rgba(223,194,125,1)"},"line":{"color":"rgba(223,194,125,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.038161275659333,-0.0421627007388314,-0.0409021428609706,-0.0424372605131672,-0.0251071471310451,-0.0300776310961218,-0.0424218404440991,-0.0341432085099119,-0.0409297339482273,-0.0384193674029452,-0.0383481495015512,-0.0415466130007104,-0.0353744846455771,-0.041322081264036,-0.0429862960628803,-0.0367741472908528,-0.0403910130036353,-0.0386611368416218,-0.0418941372149584,-0.0352231684185208,-0.0276144204171559,-0.0367865998346106,-0.0365942390142216,-0.0386956320907914,-0.0343099305214256,-0.0375158481923677,-0.0418184503600527,-0.043916895065903,-0.0434559971792787,-0.036522296581305,-0.0417863937530024,-0.0432284688412846,-0.0449171792230922,-0.0352531096814879,-0.0438554443732597,-0.0420559244963384,-0.0401321382073382,-0.038528708498348,-0.0415858740136894,-0.0417239999356786,-0.0431043284945177,-0.041346197491548,-0.034432892544808,-0.0340751409212606,-0.0450929742370307,-0.0357963191135727,-0.0453580542141781,-0.0393053327125924,-0.0250162508503767,-0.0428634836693019,-0.0425228690952905,-0.0416939069927504,-0.0400800767087653,-0.0369076493608476,-0.0376277151567743,-0.0407211595190554,-0.0424375789248122,-0.0401615667926335,-0.035951315750608,-0.0443338310494654,-0.0355626390344675,-0.0436186966860185,-0.0455104213432407,-0.0432175243805484,-0.0432947119881588,-0.0432104890859358,-0.0358863681204628,-0.0374592113216433,-0.0426174938801476,-0.0446746194853807,-0.00231477606165097,-0.0428917045999489,-0.0388310508411041,-0.0354156723664713,-0.038587576871584,-0.0357062549057946,-0.0359140516147747,-0.0392197393072526,-0.0368889701151265,-0.0409243721165635,-0.039205879312509,-0.0344529346260811,-0.0396878904783603,-0.0256453096775335,-0.0423511768820486,-0.0447719697873943,-0.0458657554297095,-0.0406350378847866,-0.036100942456952,-0.0449296633226352,-0.0405628464312889,-0.0443299939550568,-0.0430541204611418,-0.0378267520761935,-0.0416304687336589,-0.0398485229715056,-0.0377304521387802,-0.0352134871000872,-0.0370078118823983,-0.0444659413906534,-0.0428268516859975,-0.0377798107420932,-0.0407892481300821,-0.0425881822761257,-0.0443725124776152,-0.0369920012954237,-0.0372284258548125,-0.0426555115585991,-0.0405934763100242,-0.0434383848345785,-0.0433075978132937,-0.0446282273063618,-0.0338749257133586,-0.0453351167541222,-0.0328893322908196,-0.0422816511677065,-0.0447968642647755,-0.0445188732653981,-0.0466995429605546,-0.0408171231071773,-0.0399142985324884,-0.0414071110492921,-0.0404617095420034,-0.0373038099664481,-0.0433588871152374,-0.0402921191283145,-0.0345428457388193,-0.04124667119586],"y":[-0.0022913811750488,-0.00678110608660428,-0.00687869629892149,0.024244138625765,0.0826564897910184,0.00482710622628474,-0.00612622178007891,-0.0252662741331315,-0.00781666346971058,0.0757826118811694,-0.0214309746074533,-0.0271150982216078,-0.0473365810683661,0.0109009844663345,-0.028398580943622,-0.0326832586934853,-0.0035200777212032,-0.00641878299224775,-0.0179116612011385,0.0817548313386855,-0.00195642536407966,0.00283049398660965,-0.0281643604952341,0.0146026697072063,-0.00109268211434981,0.0189343615171429,-0.00656677068500975,-0.0342598772655464,0.00729501011753349,0.0247496919757275,-0.00990877202504682,-0.00195412190658987,-0.0209824093155836,0.00379201048536207,0.00170100214407012,-0.0312411951256165,-0.011208869210558,0.0122471713317324,0.0337539030173617,-0.0133783283641473,-0.0216265511169747,-0.000995647117874521,0.0299537743236343,0.0145662442367171,0.011521391375259,0.0247567315878766,0.00306269848545635,-0.0283761421804631,0.0624719562036495,-0.0206773481723813,-0.0131269207335542,-0.0346490846641029,-0.0233908383197439,-0.0110427744666152,-0.0132236687351996,-0.00395511034693245,-0.00675931049257799,0.052471448846228,0.125344348093315,-0.00435486605905008,0.138002731106726,0.000607598946701553,-0.0331690892956602,-0.0226392150346187,0.0191828600685639,-0.00932777563117196,0.0240259896272834,-0.0132669295118956,-0.012932239009143,0.0141527032527715,0.0392652129071469,-0.00498876774516954,-0.000257319698432706,-0.0195093567319892,-0.0104540395483091,0.0243190276975909,-0.0392351236687136,-0.0357169259428478,-0.0250775765432533,-0.0171669516151884,-0.0266245751869755,-0.00345365020706523,0.000260266330813794,0.0403583344268908,-0.0151992846569537,0.0121968363247761,0.00528248039466447,-0.0285279957879809,0.0752680546166768,-0.00657393842168716,0.00844834564133509,-0.0128363865004791,0.00548697871730987,-0.0120996105986602,0.00441463446591582,-0.0257682366174799,0.0199282456204836,0.0363875380148792,-0.0144209738076071,-0.0201341838153781,-0.0347576893897962,-0.0299065027285982,0.00145468191265555,0.00484380643802521,-0.033330013686318,-0.0102560114626397,-0.045106781938338,0.00828410066804319,0.0120109565251325,-0.00040800865991926,-0.0432961069461484,-0.0103219734883704,0.0737414542200591,0.00561337252112178,-0.0245740791993678,-0.0358206474670995,-0.012206736219983,0.0123054908142295,-0.0338957648231648,-0.0154775189485483,0.00126982220656942,-0.0125241666892226,-0.020246061825422,0.0935189412626638,0.0122716008701112,-0.0319591806468576,0.0575143284509431,-0.0335521895465036],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"color":"rgba(128,128,128,1)","size":11,"line":{"color":"rgba(128,128,128,1)"}},"type":"scatter","name":"Ungrouped","textfont":{"color":"rgba(128,128,128,1)"},"error_y":{"color":"rgba(128,128,128,1)"},"error_x":{"color":"rgba(128,128,128,1)"},"line":{"color":"rgba(128,128,128,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
# Reef
p <- plot_ly(df, x=~V1, y=~V2, text=rownames(df), type = "scatter",
             mode="markers", color=~reef, colors=label.colors, marker=list(size=11),
             width=w, height=h)
p <- layout(p,title="Without individual allele frequency - colored by Reef",
            xaxis=list(title="PC1"),
            yaxis=list(title="PC2"))
p
```

```{=html}
<div id="htmlwidget-d6785aab2ccde65d240c" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-d6785aab2ccde65d240c">{"x":{"visdat":{"1488458f9871e":["function () ","plotlyVisDat"]},"cur_data":"1488458f9871e","attrs":{"1488458f9871e":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009","SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747","SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#33a02c","#1f78b4","#01665e","#7570b3","#808080","#1f78b4","#33a02c","#543005","#dfc27d"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Reef","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0409021428609706,-0.0300776310961218,-0.0341432085099119,-0.0383481495015512,-0.0353744846455771,-0.0429862960628803,-0.0367741472908528,-0.043916895065903,-0.0417239999356786,-0.0450929742370307,-0.0407211595190554,-0.0424375789248122,-0.0436186966860185,-0.0432104890859358,-0.0358863681204628,-0.00231477606165097,-0.0409243721165635,-0.0396878904783603,-0.0444659413906534,-0.0428268516859975,-0.0425881822761257,-0.0443725124776152,-0.0369920012954237,-0.0426555115585991,-0.0405934763100242,-0.0433588871152374],"y":[-0.00687869629892149,0.00482710622628474,-0.0252662741331315,-0.0214309746074533,-0.0473365810683661,-0.028398580943622,-0.0326832586934853,-0.0342598772655464,-0.0133783283641473,0.011521391375259,-0.00395511034693245,-0.00675931049257799,0.000607598946701553,-0.00932777563117196,0.0240259896272834,0.0392652129071469,-0.0171669516151884,0.000260266330813794,-0.0201341838153781,-0.0347576893897962,0.00484380643802521,-0.033330013686318,-0.0102560114626397,0.00828410066804319,0.0120109565251325,0.0122716008701112],"text":["Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP9_1121","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP9_2862","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1997","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP8_2067"],"mode":"markers","marker":{"color":"rgba(179,88,6,1)","size":11,"line":{"color":"rgba(179,88,6,1)"}},"type":"scatter","name":"HIMB","textfont":{"color":"rgba(179,88,6,1)"},"error_y":{"color":"rgba(179,88,6,1)"},"error_x":{"color":"rgba(179,88,6,1)"},"line":{"color":"rgba(179,88,6,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.145997445819095,0.174955367864438,0.174902797925439,0.151102413283901,0.159919495340554,0.157942533764846,0.146420407313794,0.173851484550974,0.176377286930581],"y":[0.192883367283657,0.158032196681908,0.158940044226475,0.121758939195166,0.124271170183639,0.12000084654977,0.191061022008456,0.152828046521998,0.158656760250303],"text":["SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747"],"mode":"markers","marker":{"color":"rgba(51,160,44,1)","size":11,"line":{"color":"rgba(51,160,44,1)"}},"type":"scatter","name":"Kiholo","textfont":{"color":"rgba(51,160,44,1)"},"error_y":{"color":"rgba(51,160,44,1)"},"error_x":{"color":"rgba(51,160,44,1)"},"line":{"color":"rgba(51,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0367865998346106,-0.0432175243805484,-0.0426174938801476,-0.0344529346260811,-0.0430541204611418,-0.0370078118823983,-0.0328893322908196,-0.04124667119586],"y":[0.00283049398660965,-0.0226392150346187,-0.012932239009143,-0.00345365020706523,0.00548697871730987,-0.0144209738076071,-0.0245740791993678,-0.0335521895465036],"text":["Mcapitata_ATAC_TP7_1058","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP9_1306","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"color":"rgba(241,163,64,1)","size":11,"line":{"color":"rgba(241,163,64,1)"}},"type":"scatter","name":"Lilipuna.Fringe","textfont":{"color":"rgba(241,163,64,1)"},"error_y":{"color":"rgba(241,163,64,1)"},"error_x":{"color":"rgba(241,163,64,1)"},"line":{"color":"rgba(241,163,64,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0421627007388314,-0.0415466130007104,-0.0417863937530024,-0.0420559244963384,-0.038528708498348,-0.0428634836693019,-0.0416939069927504,-0.0455104213432407,-0.0432947119881588,-0.0374592113216433,-0.039205879312509,-0.0458657554297095,-0.0378267520761935,-0.0398485229715056,-0.0377304521387802,-0.0352134871000872,-0.0446282273063618,-0.0422816511677065,-0.0466995429605546,-0.0404617095420034,-0.0402921191283145,-0.0345428457388193],"y":[-0.00678110608660428,-0.0271150982216078,-0.00990877202504682,-0.0312411951256165,0.0122471713317324,-0.0206773481723813,-0.0346490846641029,-0.0331690892956602,0.0191828600685639,-0.0132669295118956,-0.0266245751869755,0.00528248039466447,-0.0120996105986602,-0.0257682366174799,0.0199282456204836,0.0363875380148792,-0.0103219734883704,-0.0358206474670995,-0.0338957648231648,-0.020246061825422,-0.0319591806468576,0.0575143284509431],"text":["Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP3_1548","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP9_1467","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331"],"mode":"markers","marker":{"color":"rgba(254,224,182,1)","size":11,"line":{"color":"rgba(254,224,182,1)"}},"type":"scatter","name":"Reef.11.13","textfont":{"color":"rgba(254,224,182,1)"},"error_y":{"color":"rgba(254,224,182,1)"},"error_x":{"color":"rgba(254,224,182,1)"},"line":{"color":"rgba(254,224,182,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0301690858687305,-0.0444844359594012,-0.0403910130036353,-0.0432579249161303,-0.0432284688412846,-0.041346197491548,-0.034432892544808,-0.0453580542141781,-0.0369076493608476,-0.0308583956418033,-0.0446746194853807,-0.0428917045999489,-0.0423511768820486,-0.0405628464312889,-0.0443299939550568,-0.0416304687336589,-0.0377798107420932,-0.0434383848345785,-0.0453351167541222,-0.0414071110492921],"y":[0.0384051391045254,-0.0624821710739555,-0.0035200777212032,-0.0635391652669186,-0.00195412190658987,-0.000995647117874521,0.0299537743236343,0.00306269848545635,-0.0110427744666152,0.0379059230536485,0.0141527032527715,-0.00498876774516954,-0.0151992846569537,0.00844834564133509,-0.0128363865004791,0.00441463446591582,-0.0299065027285982,-0.00040800865991926,0.00561337252112178,-0.0125241666892226],"text":["Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP6_2402","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP6_2555","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP7_2419"],"mode":"markers","marker":{"color":"rgba(216,218,235,1)","size":11,"line":{"color":"rgba(216,218,235,1)"}},"type":"scatter","name":"Reef.18","textfont":{"color":"rgba(216,218,235,1)"},"error_y":{"color":"rgba(216,218,235,1)"},"error_x":{"color":"rgba(216,218,235,1)"},"line":{"color":"rgba(216,218,235,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.038161275659333,-0.0424372605131672,-0.0251071471310451,-0.0384193674029452,-0.041322081264036,-0.0365942390142216,-0.0386956320907914,-0.0343099305214256,-0.0375158481923677,-0.0434559971792787,-0.036522296581305,-0.0449171792230922,-0.0352531096814879,-0.0438554443732597,-0.0415858740136894,-0.0357963191135727,-0.0425228690952905,-0.0400800767087653,-0.0401615667926335,-0.038587576871584,-0.0392197393072526,-0.0368889701151265,-0.0256453096775335,-0.0447719697873943,-0.036100942456952,-0.0407892481300821,-0.0447968642647755,-0.0445188732653981],"y":[-0.0022913811750488,0.024244138625765,0.0826564897910184,0.0757826118811694,0.0109009844663345,-0.0281643604952341,0.0146026697072063,-0.00109268211434981,0.0189343615171429,0.00729501011753349,0.0247496919757275,-0.0209824093155836,0.00379201048536207,0.00170100214407012,0.0337539030173617,0.0247567315878766,-0.0131269207335542,-0.0233908383197439,0.052471448846228,-0.0104540395483091,-0.0357169259428478,-0.0250775765432533,0.0403583344268908,0.0121968363247761,0.0752680546166768,0.00145468191265555,-0.012206736219983,0.0123054908142295],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP8_1260","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1722","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317"],"mode":"markers","marker":{"color":"rgba(153,142,195,1)","size":11,"line":{"color":"rgba(153,142,195,1)"}},"type":"scatter","name":"Reef.35.36","textfont":{"color":"rgba(153,142,195,1)"},"error_y":{"color":"rgba(153,142,195,1)"},"error_x":{"color":"rgba(153,142,195,1)"},"line":{"color":"rgba(153,142,195,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0424218404440991,-0.0409297339482273,-0.0386611368416218,-0.0418941372149584,-0.0352231684185208,-0.0276144204171559,-0.0418184503600527,-0.0401321382073382,-0.0431043284945177,-0.0340751409212606,-0.0393053327125924,-0.0250162508503767,-0.0376277151567743,-0.035951315750608,-0.0443338310494654,-0.0355626390344675,-0.0388310508411041,-0.0354156723664713,-0.0357062549057946,-0.0359140516147747,-0.0406350378847866,-0.0449296633226352,-0.0372284258548125,-0.0433075978132937,-0.0338749257133586,-0.0408171231071773,-0.0399142985324884,-0.0373038099664481],"y":[-0.00612622178007891,-0.00781666346971058,-0.00641878299224775,-0.0179116612011385,0.0817548313386855,-0.00195642536407966,-0.00656677068500975,-0.011208869210558,-0.0216265511169747,0.0145662442367171,-0.0283761421804631,0.0624719562036495,-0.0132236687351996,0.125344348093315,-0.00435486605905008,0.138002731106726,-0.000257319698432706,-0.0195093567319892,0.0243190276975909,-0.0392351236687136,-0.0285279957879809,-0.00657393842168716,-0.045106781938338,-0.0432961069461484,0.0737414542200591,-0.0154775189485483,0.00126982220656942,0.0935189412626638],"text":["Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP8_1779","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP8_1235","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP8_1246"],"mode":"markers","marker":{"color":"rgba(84,39,136,1)","size":11,"line":{"color":"rgba(84,39,136,1)"}},"type":"scatter","name":"Reef.42.43","textfont":{"color":"rgba(84,39,136,1)"},"error_y":{"color":"rgba(84,39,136,1)"},"error_x":{"color":"rgba(84,39,136,1)"},"line":{"color":"rgba(84,39,136,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.192963783550781,0.187604556302613,0.195375556575114,0.173015234268865,0.174672358864291,0.184831419673122,0.169261625046705,0.161388905503109,0.160889770031219,0.15102065362944,0.178247124058197,0.179315378194888,0.194924584218308,0.188394179356329,0.195885141133129,0.142341074734009,0.155785419073372,0.154281672744728],"y":[-0.285872345487089,-0.169758811039003,-0.365188395712793,0.0809955243749636,0.122361962615188,-0.1038389960395,0.0525769728321265,0.0742335665983742,0.0191980288362958,0.0754682482966374,0.0819136645421538,0.0788716541902661,-0.343242651831009,-0.174828341971528,-0.368863867281678,0.105571228632626,0.146926092689027,-0.131702627315251],"text":["SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"color":"rgba(31,120,180,1)","size":11,"line":{"color":"rgba(31,120,180,1)"}},"type":"scatter","name":"Waiopae","textfont":{"color":"rgba(31,120,180,1)"},"error_y":{"color":"rgba(31,120,180,1)"},"error_x":{"color":"rgba(31,120,180,1)"},"line":{"color":"rgba(31,120,180,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
# Ploidy
p <- plot_ly(df, x=~V1, y=~V2, text=rownames(df), type = "scatter",
             mode="markers", color=~ploidy, colors=label.colors, marker=list(size=11),
             width=w, height=h)
p <- layout(p,title="Without individual allele frequency - colored by Ploidy",
            xaxis=list(title="PC1"),
            yaxis=list(title="PC2"))
p
```

```{=html}
<div id="htmlwidget-9277f511b41f371ecdb8" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-9277f511b41f371ecdb8">{"x":{"visdat":{"148846615f034":["function () ","plotlyVisDat"]},"cur_data":"148846615f034","attrs":{"148846615f034":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009","SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747","SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#33a02c","#1f78b4","#01665e","#7570b3","#808080","#1f78b4","#33a02c","#543005","#dfc27d"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Ploidy","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.038161275659333,-0.0421627007388314,-0.0409021428609706,-0.0424372605131672,-0.0251071471310451,-0.0300776310961218,-0.0301690858687305,-0.0424218404440991,-0.0341432085099119,-0.0409297339482273,-0.0384193674029452,-0.0444844359594012,-0.0383481495015512,-0.0415466130007104,-0.0353744846455771,-0.041322081264036,-0.0429862960628803,-0.0367741472908528,-0.0403910130036353,-0.0386611368416218,-0.0418941372149584,-0.0352231684185208,-0.0276144204171559,-0.0432579249161303,-0.0367865998346106,-0.0365942390142216,-0.0386956320907914,-0.0343099305214256,-0.0375158481923677,-0.0418184503600527,-0.043916895065903,-0.0434559971792787,-0.036522296581305,-0.0417863937530024,-0.0432284688412846,-0.0449171792230922,-0.0352531096814879,-0.0438554443732597,-0.0420559244963384,-0.0401321382073382,-0.038528708498348,-0.0415858740136894,-0.0417239999356786,-0.0431043284945177,-0.041346197491548,-0.034432892544808,-0.0340751409212606,-0.0450929742370307,-0.0357963191135727,-0.0453580542141781,-0.0393053327125924,-0.0250162508503767,-0.0428634836693019,-0.0425228690952905,-0.0416939069927504,-0.0400800767087653,-0.0369076493608476,-0.0376277151567743,-0.0407211595190554,-0.0424375789248122,-0.0401615667926335,-0.035951315750608,-0.0443338310494654,-0.0355626390344675,-0.0436186966860185,-0.0455104213432407,-0.0432175243805484,-0.0432947119881588,-0.0432104890859358,-0.0358863681204628,-0.0374592113216433,-0.0308583956418033,-0.0426174938801476,-0.0446746194853807,-0.00231477606165097,-0.0428917045999489,-0.0388310508411041,-0.0354156723664713,-0.038587576871584,-0.0357062549057946,-0.0359140516147747,-0.0392197393072526,-0.0368889701151265,-0.0409243721165635,-0.039205879312509,-0.0344529346260811,-0.0396878904783603,-0.0256453096775335,-0.0423511768820486,-0.0447719697873943,-0.0458657554297095,-0.0406350378847866,-0.036100942456952,-0.0449296633226352,-0.0405628464312889,-0.0443299939550568,-0.0430541204611418,-0.0378267520761935,-0.0416304687336589,-0.0398485229715056,-0.0377304521387802,-0.0352134871000872,-0.0370078118823983,-0.0444659413906534,-0.0428268516859975,-0.0377798107420932,-0.0407892481300821,-0.0425881822761257,-0.0443725124776152,-0.0369920012954237,-0.0372284258548125,-0.0426555115585991,-0.0405934763100242,-0.0434383848345785,-0.0433075978132937,-0.0446282273063618,-0.0338749257133586,-0.0453351167541222,-0.0328893322908196,-0.0422816511677065,-0.0447968642647755,-0.0445188732653981,-0.0466995429605546,-0.0408171231071773,-0.0399142985324884,-0.0414071110492921,-0.0404617095420034,-0.0373038099664481,-0.0433588871152374,-0.0402921191283145,-0.0345428457388193,-0.04124667119586,0.145997445819095,0.174955367864438,0.174902797925439,0.151102413283901,0.159919495340554,0.157942533764846,0.146420407313794,0.173851484550974,0.176377286930581,0.192963783550781,0.187604556302613,0.195375556575114,0.173015234268865,0.174672358864291,0.184831419673122,0.169261625046705,0.161388905503109,0.160889770031219,0.15102065362944,0.178247124058197,0.179315378194888,0.194924584218308,0.188394179356329,0.195885141133129,0.142341074734009,0.155785419073372,0.154281672744728],"y":[-0.0022913811750488,-0.00678110608660428,-0.00687869629892149,0.024244138625765,0.0826564897910184,0.00482710622628474,0.0384051391045254,-0.00612622178007891,-0.0252662741331315,-0.00781666346971058,0.0757826118811694,-0.0624821710739555,-0.0214309746074533,-0.0271150982216078,-0.0473365810683661,0.0109009844663345,-0.028398580943622,-0.0326832586934853,-0.0035200777212032,-0.00641878299224775,-0.0179116612011385,0.0817548313386855,-0.00195642536407966,-0.0635391652669186,0.00283049398660965,-0.0281643604952341,0.0146026697072063,-0.00109268211434981,0.0189343615171429,-0.00656677068500975,-0.0342598772655464,0.00729501011753349,0.0247496919757275,-0.00990877202504682,-0.00195412190658987,-0.0209824093155836,0.00379201048536207,0.00170100214407012,-0.0312411951256165,-0.011208869210558,0.0122471713317324,0.0337539030173617,-0.0133783283641473,-0.0216265511169747,-0.000995647117874521,0.0299537743236343,0.0145662442367171,0.011521391375259,0.0247567315878766,0.00306269848545635,-0.0283761421804631,0.0624719562036495,-0.0206773481723813,-0.0131269207335542,-0.0346490846641029,-0.0233908383197439,-0.0110427744666152,-0.0132236687351996,-0.00395511034693245,-0.00675931049257799,0.052471448846228,0.125344348093315,-0.00435486605905008,0.138002731106726,0.000607598946701553,-0.0331690892956602,-0.0226392150346187,0.0191828600685639,-0.00932777563117196,0.0240259896272834,-0.0132669295118956,0.0379059230536485,-0.012932239009143,0.0141527032527715,0.0392652129071469,-0.00498876774516954,-0.000257319698432706,-0.0195093567319892,-0.0104540395483091,0.0243190276975909,-0.0392351236687136,-0.0357169259428478,-0.0250775765432533,-0.0171669516151884,-0.0266245751869755,-0.00345365020706523,0.000260266330813794,0.0403583344268908,-0.0151992846569537,0.0121968363247761,0.00528248039466447,-0.0285279957879809,0.0752680546166768,-0.00657393842168716,0.00844834564133509,-0.0128363865004791,0.00548697871730987,-0.0120996105986602,0.00441463446591582,-0.0257682366174799,0.0199282456204836,0.0363875380148792,-0.0144209738076071,-0.0201341838153781,-0.0347576893897962,-0.0299065027285982,0.00145468191265555,0.00484380643802521,-0.033330013686318,-0.0102560114626397,-0.045106781938338,0.00828410066804319,0.0120109565251325,-0.00040800865991926,-0.0432961069461484,-0.0103219734883704,0.0737414542200591,0.00561337252112178,-0.0245740791993678,-0.0358206474670995,-0.012206736219983,0.0123054908142295,-0.0338957648231648,-0.0154775189485483,0.00126982220656942,-0.0125241666892226,-0.020246061825422,0.0935189412626638,0.0122716008701112,-0.0319591806468576,0.0575143284509431,-0.0335521895465036,0.192883367283657,0.158032196681908,0.158940044226475,0.121758939195166,0.124271170183639,0.12000084654977,0.191061022008456,0.152828046521998,0.158656760250303,-0.285872345487089,-0.169758811039003,-0.365188395712793,0.0809955243749636,0.122361962615188,-0.1038389960395,0.0525769728321265,0.0742335665983742,0.0191980288362958,0.0754682482966374,0.0819136645421538,0.0788716541902661,-0.343242651831009,-0.174828341971528,-0.368863867281678,0.105571228632626,0.146926092689027,-0.131702627315251],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009","SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747","SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"colorbar":{"title":"ploidy","ticklen":2},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(214,189,189,1)"],["0.0833333333333333","rgba(244,201,153,1)"],["0.125","rgba(201,121,56,1)"],["0.166666666666667","rgba(201,130,88,1)"],["0.208333333333333","rgba(217,207,216,1)"],["0.25","rgba(152,126,185,1)"],["0.291666666666667","rgba(102,48,132,1)"],["0.333333333333333","rgba(195,119,97,1)"],["0.375","rgba(202,165,59,1)"],["0.416666666666667","rgba(97,162,47,1)"],["0.458333333333333","rgba(68,143,107,1)"],["0.5","rgba(31,120,180,1)"],["0.541666666666667","rgba(23,110,129,1)"],["0.583333333333333","rgba(38,104,108,1)"],["0.625","rgba(96,111,157,1)"],["0.666666666666667","rgba(122,117,162,1)"],["0.708333333333333","rgba(128,127,132,1)"],["0.75","rgba(97,124,154,1)"],["0.791666666666667","rgba(45,123,170,1)"],["0.833333333333333","rgba(67,146,96,1)"],["0.875","rgba(74,132,33,1)"],["0.916666666666667","rgba(86,67,11,1)"],["0.958333333333333","rgba(140,105,53,1)"],["1","rgba(223,194,125,1)"]],"showscale":false,"color":[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"size":11,"line":{"colorbar":{"title":"","ticklen":2},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(214,189,189,1)"],["0.0833333333333333","rgba(244,201,153,1)"],["0.125","rgba(201,121,56,1)"],["0.166666666666667","rgba(201,130,88,1)"],["0.208333333333333","rgba(217,207,216,1)"],["0.25","rgba(152,126,185,1)"],["0.291666666666667","rgba(102,48,132,1)"],["0.333333333333333","rgba(195,119,97,1)"],["0.375","rgba(202,165,59,1)"],["0.416666666666667","rgba(97,162,47,1)"],["0.458333333333333","rgba(68,143,107,1)"],["0.5","rgba(31,120,180,1)"],["0.541666666666667","rgba(23,110,129,1)"],["0.583333333333333","rgba(38,104,108,1)"],["0.625","rgba(96,111,157,1)"],["0.666666666666667","rgba(122,117,162,1)"],["0.708333333333333","rgba(128,127,132,1)"],["0.75","rgba(97,124,154,1)"],["0.791666666666667","rgba(45,123,170,1)"],["0.833333333333333","rgba(67,146,96,1)"],["0.875","rgba(74,132,33,1)"],["0.916666666666667","rgba(86,67,11,1)"],["0.958333333333333","rgba(140,105,53,1)"],["1","rgba(223,194,125,1)"]],"showscale":false,"color":[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}},"type":"scatter","xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0466995429605546,0.195885141133129],"y":[-0.368863867281678,0.192883367283657],"type":"scatter","mode":"markers","opacity":0,"hoverinfo":"none","showlegend":false,"marker":{"colorbar":{"title":"ploidy","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(214,189,189,1)"],["0.0833333333333333","rgba(244,201,153,1)"],["0.125","rgba(201,121,56,1)"],["0.166666666666667","rgba(201,130,88,1)"],["0.208333333333333","rgba(217,207,216,1)"],["0.25","rgba(152,126,185,1)"],["0.291666666666667","rgba(102,48,132,1)"],["0.333333333333333","rgba(195,119,97,1)"],["0.375","rgba(202,165,59,1)"],["0.416666666666667","rgba(97,162,47,1)"],["0.458333333333333","rgba(68,143,107,1)"],["0.5","rgba(31,120,180,1)"],["0.541666666666667","rgba(23,110,129,1)"],["0.583333333333333","rgba(38,104,108,1)"],["0.625","rgba(96,111,157,1)"],["0.666666666666667","rgba(122,117,162,1)"],["0.708333333333333","rgba(128,127,132,1)"],["0.75","rgba(97,124,154,1)"],["0.791666666666667","rgba(45,123,170,1)"],["0.833333333333333","rgba(67,146,96,1)"],["0.875","rgba(74,132,33,1)"],["0.916666666666667","rgba(86,67,11,1)"],["0.958333333333333","rgba(140,105,53,1)"],["1","rgba(223,194,125,1)"]],"showscale":true,"color":[2,4],"line":{"color":"rgba(255,127,14,1)"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```



Plot PC1 and PC3 of PCAngsd results with estimated individual allele frequencies.

```r
A <- 1
B <- 3

plot(e$vectors[,c(A,B)], col=pop$group_color,  pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Group")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC3-1.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$reef_color,   pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Reef")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC3-2.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$ploidy_color, pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Ploidy")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC3-3.png)<!-- -->



Plot PC1 and PC4 of PCAngsd results with estimated individual allele frequencies.

```r
A <- 1
B <- 4

plot(e$vectors[,c(A,B)], col=pop$group_color,  pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Group")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC4-1.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$reef_color,   pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Reef")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC4-2.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$ploidy_color, pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Ploidy")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC1_PC4-3.png)<!-- -->



Plot PC2 and PC3 of PCAngsd results with estimated individual allele frequencies.

```r
A <- 2
B <- 3

plot(e$vectors[,c(A,B)], col=pop$group_color,  pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Group")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC2_PC3-1.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$reef_color,   pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Reef")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC2_PC3-2.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$ploidy_color, pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Ploidy")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC2_PC3-3.png)<!-- -->



Plot PC2 and PC4 of PCAngsd results with estimated individual allele frequencies.

```r
A <- 2
B <- 4

plot(e$vectors[,c(A,B)], col=pop$group_color,  pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Group")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC2_PC4-1.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$reef_color,   pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Reef")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC2_PC4-2.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$ploidy_color, pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Individual allele frequency - colored by Ploidy")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_IndAllelFreq_PC2_PC4-3.png)<!-- -->





# Without Estimating Individual Allele Frequencies

Load PCAngsd results with WITHOUT estimated individual allele frequencies.

```r
C <- as.matrix(read.table("PCAngsd.angsd.beagle.gz.WithOutIndAlleleFreq.cov"))
e <- eigen(C)
pop <- samples.info
```



Plot PC1 and PC2 of PCAngsd results WITHOUT estimated individual allele frequencies.

```r
A <- 1
B <- 2

plot(e$vectors[,c(A,B)], col=pop$group_color,  pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Without individual allele frequency - colored by Group")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_withoutIndAllelFreq_PC1_PC2-1.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$reef_color,   pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Without individual allele frequency - colored by Reef")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_withoutIndAllelFreq_PC1_PC2-2.png)<!-- -->

```r
plot(e$vectors[,c(A,B)], col=pop$ploidy_color, pch=19, xlab=paste("PC", A, sep=''), ylab=paste("PC", B, sep=''), main="Without individual allele frequency - colored by Ploidy")
```

![](plot_ANGSD_results_files/figure-html/plot_PCAngsd_withoutIndAllelFreq_PC1_PC2-3.png)<!-- -->



Plot PC1 and PC2 of PCAngsd results WITHOUT estimated individual allele frequencies using the `plotly` package. 

```r
df <- as.data.frame(e$vectors)
df <- cbind(df, samples.info)

# Group
p <- plot_ly(df, x=~V1, y=~V2, text=rownames(df), type = "scatter",
             mode="markers", color=~group, colors=label.colors, marker=list(size=11))
p <- layout(p,title="Without individual allele frequency - colored by Group",
            xaxis=list(title="PC1"),
            yaxis=list(title="PC2"))
p
```

```{=html}
<div id="htmlwidget-23f587f6db132e18a0d0" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-23f587f6db132e18a0d0">{"x":{"visdat":{"1488428d5f8a0":["function () ","plotlyVisDat"]},"cur_data":"1488428d5f8a0","attrs":{"1488428d5f8a0":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009","SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747","SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#33a02c","#1f78b4","#01665e","#7570b3","#808080","#1f78b4","#33a02c","#543005","#dfc27d"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Group","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0292450117061691,-0.0312727343293337],"y":[-0.0132485412069704,-0.0139800210928299],"text":["Mcapitata_ATAC_TP11_2302","Mcapitata_HTAC_TP11_2380"],"mode":"markers","marker":{"color":"rgba(31,120,180,1)","size":11,"line":{"color":"rgba(31,120,180,1)"}},"type":"scatter","name":"Group1","textfont":{"color":"rgba(31,120,180,1)"},"error_y":{"color":"rgba(31,120,180,1)"},"error_x":{"color":"rgba(31,120,180,1)"},"line":{"color":"rgba(31,120,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0446426063690368,-0.0438216733551513],"y":[0.0194529170663445,0.0192432959893718],"text":["Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP6_2402"],"mode":"markers","marker":{"color":"rgba(51,160,44,1)","size":11,"line":{"color":"rgba(51,160,44,1)"}},"type":"scatter","name":"Group2","textfont":{"color":"rgba(51,160,44,1)"},"error_y":{"color":"rgba(51,160,44,1)"},"error_x":{"color":"rgba(51,160,44,1)"},"line":{"color":"rgba(51,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.101947155132633,0.169925374298519,0.151445130101844,0.118926735652195,0.140757365595538,0.123406412242653,0.100232185621671,0.151188360915327,0.155688253938549],"y":[-0.130753842904154,-0.211321869131487,-0.196011935514122,-0.12345961880238,-0.138581199128536,-0.128053118706528,-0.128096276620108,-0.198250962083971,-0.197472331674013],"text":["SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747"],"mode":"markers","marker":{"color":"rgba(84,48,5,1)","size":11,"line":{"color":"rgba(84,48,5,1)"}},"type":"scatter","name":"SRA-Kiholo","textfont":{"color":"rgba(84,48,5,1)"},"error_y":{"color":"rgba(84,48,5,1)"},"error_x":{"color":"rgba(84,48,5,1)"},"line":{"color":"rgba(84,48,5,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.168994877569605,0.21509720342966,0.254568411526922,0.138856275324994,0.171702481377677,0.195859199686892,0.180539848936965,0.16424092655158,0.149246060042661,0.121203952823388,0.192417055323507,0.178056638666907,0.23000314893794,0.224335537594743,0.259111940141535,0.158023529046395,0.200859784931242,0.183815309324472],"y":[0.295345772294402,0.152155892960803,0.398763226164787,-0.0992143817114808,-0.167647580548066,0.0450558662621251,-0.0898255568135078,-0.102820139293431,-0.0547096829421381,-0.0800307901453591,-0.133405739687102,-0.128401681801111,0.366329662334079,0.158704545595299,0.404324290130019,-0.108716729952454,-0.184837002934645,0.0471436397166563],"text":["SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"color":"rgba(223,194,125,1)","size":11,"line":{"color":"rgba(223,194,125,1)"}},"type":"scatter","name":"SRA-Waiopae","textfont":{"color":"rgba(223,194,125,1)"},"error_y":{"color":"rgba(223,194,125,1)"},"error_x":{"color":"rgba(223,194,125,1)"},"line":{"color":"rgba(223,194,125,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0259984393575133,-0.0353283383245678,-0.0323376640238664,-0.0410654394546692,-0.024818382982012,-0.0235045784679117,-0.0293737093915074,-0.0334676804562726,-0.0364525034473797,-0.0353034318908815,-0.038112810285706,-0.0389494957653096,-0.0329575914573595,-0.0394449374274202,-0.0436701159998899,-0.036332483373355,-0.0333414861309239,-0.0339195879748338,-0.0247921975595541,-0.0281402178675904,-0.0269030375756818,-0.0306404187791171,-0.0351041069986101,-0.0308549605672989,-0.0329544729411579,-0.0292574405201576,-0.0326241022925872,-0.0447606336039882,-0.0347188028316388,-0.0280177223372142,-0.0387156994742788,-0.034731950916434,-0.0434408995316462,-0.032692920164491,-0.0400282426133999,-0.0417337881359568,-0.0370454974061555,-0.0304936991285582,-0.0345190745745488,-0.0400429811263431,-0.0384514501961568,-0.0457030920108826,-0.0311424153153421,-0.0295116367933613,-0.0377975735649241,-0.0311168811525216,-0.036073349213048,-0.0398474851566358,-0.0259358590739191,-0.0402493856659681,-0.0452789332863361,-0.0431807410559234,-0.036626004446606,-0.0370520983792864,-0.0309172669884343,-0.0353843096320885,-0.0382516942893998,-0.0298324067335494,-0.0319312162365303,-0.0448675828280251,-0.0281309662843891,-0.0376159759943272,-0.0412446884480399,-0.0399493733930138,-0.0323077123296758,-0.0346922932378603,-0.0347757908668086,-0.0323721513550889,-0.0429674931170672,-0.0473494065380055,-0.0104889376097138,-0.041164967807276,-0.0339820363209966,-0.031368639082219,-0.0334591162056912,-0.0262667679351678,-0.0318843306819061,-0.0382554432390951,-0.0350690329703446,-0.0424267530018416,-0.0381188865102576,-0.0304163527300572,-0.0348077090040011,-0.0205173239207203,-0.0376882267183699,-0.0402284760493714,-0.0405122397483136,-0.038279637555649,-0.0338310961146996,-0.0317292179162127,-0.0310620019777845,-0.0453449334759705,-0.0319617373876199,-0.0353802090917302,-0.0352108274111176,-0.0298078151036683,-0.0388331480283056,-0.034304902419877,-0.0334669177230852,-0.0455368930964151,-0.0417314134844289,-0.0377806154522496,-0.0340804430072334,-0.0311981113407001,-0.0417398345418761,-0.0355274227864432,-0.0317196716412144,-0.0327391327721698,-0.0374051471455735,-0.0351359895664081,-0.0416876422450806,-0.040824021862402,-0.0298857743800826,-0.0400879968629934,-0.0308552264537443,-0.0392288773251657,-0.0349431404266913,-0.0368868925529881,-0.0495706249695658,-0.0340432555309663,-0.0316155425076258,-0.0414191317913502,-0.03872879426252,-0.0268278657682002,-0.0453006255320472,-0.0386128238993237,-0.0245392315759141,-0.0390200030270299],"y":[0.00752789801107825,0.0086215729500388,0.0069021802979296,-0.00371152102354799,-0.000697769477866483,0.00402462237005923,0.00740082660305173,0.00882482895190268,0.00371171966129259,-0.00186750481337353,0.0103064901567753,0.00777338915666925,0.0208961644592763,-0.00152714811752618,0.0113391097714042,0.00891407969436935,0.00451228042831412,0.00210165612485316,0.0123735078487645,-0.00680121588758422,0.00269733257968474,0.00945045541406599,0.0121572380563792,-0.00343873854075542,0.00167275413307834,-0.001080879835058,0.0124717396585872,0.0122609393978024,0.000261492778772618,-0.000653872480533765,0.0040734750542804,0.00598768670966988,0.0128910671476036,-0.00166332417244114,0.00252221481742771,0.0109235427064767,0.00982664005711287,0.000799730749840862,-0.00176506228597524,0.00739500248127063,0.0140393765263796,-0.000666676311828346,-0.00740119524011746,0.00473295591764761,0.00185949292475981,0.0091681429556099,0.0019450778832346,0.00947956277800932,-0.00285686968435638,0.0117106593918568,0.00684506168998994,0.012444694807492,0.00999648517240926,0.00992521138326167,0.0102220592693368,0.00460787537043577,0.00816602244194009,0.000673671608306822,-0.0104454792593248,0.00479853733796544,-0.0117125282918918,0.00618920896918417,0.0210155148565302,0.0138902090716542,0.0127221641878477,0.011215012688699,-0.00978118197248325,0.00831680934603369,0.00868792229378588,0.00568786750382434,-0.0110622926807172,0.0060824885581495,0.00199651447173441,0.00658577111237137,0.00635245914950481,-0.00420618113072376,0.0180438290013235,0.0147113976353928,0.00902944719518959,0.00264298440373561,0.00892950422410608,-0.00302924325427964,-0.000496565143664624,0.00673634758945498,0.00991678192004881,0.00626057772248626,0.000152537396615141,0.0130500696695487,-0.0119778351543156,0.0097603154605501,0.00264683161542478,0.00374011002691165,0.00158844631698647,0.0082660743482269,0.00198074219032276,0.0146984362206223,-0.0138922531290385,0.00655896064561861,0.007907944782137,0.0105221320713179,0.00959425609874879,0.014039578065848,0.00510976465980565,-0.000369487216460847,0.0112162997085473,0.00185025766104319,0.0229341965185012,-0.00273864713828892,0.0110126232167368,0.000500250439251407,0.0180986407976146,0.00132292460985308,-0.00920344431416447,-0.00181861048194446,0.0108432929135499,0.0162731219635841,0.0121706893305822,-0.000774896722339206,0.0158125512291628,0.0106575189293047,0.00480025396540249,0.00471743216220084,0.00783451889028733,0.00133554972561204,0.00241266999768364,0.0160110514242019,-0.00155570482177069,0.0120720381022763],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"color":"rgba(128,128,128,1)","size":11,"line":{"color":"rgba(128,128,128,1)"}},"type":"scatter","name":"Ungrouped","textfont":{"color":"rgba(128,128,128,1)"},"error_y":{"color":"rgba(128,128,128,1)"},"error_x":{"color":"rgba(128,128,128,1)"},"line":{"color":"rgba(128,128,128,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
# Reef
p <- plot_ly(df, x=~V1, y=~V2, text=rownames(df), type = "scatter",
             mode="markers", color=~reef, colors=label.colors, marker=list(size=11))
p <- layout(p,title="Without individual allele frequency - colored by Reef",
            xaxis=list(title="PC1"),
            yaxis=list(title="PC2"))
p
```

```{=html}
<div id="htmlwidget-b04e1dbac71850082b5f" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-b04e1dbac71850082b5f">{"x":{"visdat":{"1488477b5e551":["function () ","plotlyVisDat"]},"cur_data":"1488477b5e551","attrs":{"1488477b5e551":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009","SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747","SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#33a02c","#1f78b4","#01665e","#7570b3","#808080","#1f78b4","#33a02c","#543005","#dfc27d"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Reef","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0323376640238664,-0.0235045784679117,-0.0334676804562726,-0.038112810285706,-0.0329575914573595,-0.0436701159998899,-0.036332483373355,-0.0447606336039882,-0.0400429811263431,-0.0377975735649241,-0.0353843096320885,-0.0382516942893998,-0.0376159759943272,-0.0346922932378603,-0.0347757908668086,-0.0104889376097138,-0.0424267530018416,-0.0348077090040011,-0.0455368930964151,-0.0417314134844289,-0.0311981113407001,-0.0417398345418761,-0.0355274227864432,-0.0327391327721698,-0.0374051471455735,-0.0453006255320472],"y":[0.0069021802979296,0.00402462237005923,0.00882482895190268,0.0103064901567753,0.0208961644592763,0.0113391097714042,0.00891407969436935,0.0122609393978024,0.00739500248127063,0.00185949292475981,0.00460787537043577,0.00816602244194009,0.00618920896918417,0.011215012688699,-0.00978118197248325,-0.0110622926807172,0.00264298440373561,-0.000496565143664624,0.0105221320713179,0.00959425609874879,-0.000369487216460847,0.0112162997085473,0.00185025766104319,-0.00273864713828892,0.0110126232167368,0.00241266999768364],"text":["Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP9_1121","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP9_2862","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1997","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP8_2067"],"mode":"markers","marker":{"color":"rgba(179,88,6,1)","size":11,"line":{"color":"rgba(179,88,6,1)"}},"type":"scatter","name":"HIMB","textfont":{"color":"rgba(179,88,6,1)"},"error_y":{"color":"rgba(179,88,6,1)"},"error_x":{"color":"rgba(179,88,6,1)"},"line":{"color":"rgba(179,88,6,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.101947155132633,0.169925374298519,0.151445130101844,0.118926735652195,0.140757365595538,0.123406412242653,0.100232185621671,0.151188360915327,0.155688253938549],"y":[-0.130753842904154,-0.211321869131487,-0.196011935514122,-0.12345961880238,-0.138581199128536,-0.128053118706528,-0.128096276620108,-0.198250962083971,-0.197472331674013],"text":["SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747"],"mode":"markers","marker":{"color":"rgba(51,160,44,1)","size":11,"line":{"color":"rgba(51,160,44,1)"}},"type":"scatter","name":"Kiholo","textfont":{"color":"rgba(51,160,44,1)"},"error_y":{"color":"rgba(51,160,44,1)"},"error_x":{"color":"rgba(51,160,44,1)"},"line":{"color":"rgba(51,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0306404187791171,-0.0399493733930138,-0.0429674931170672,-0.0304163527300572,-0.0319617373876199,-0.0334669177230852,-0.0308552264537443,-0.0390200030270299],"y":[0.00945045541406599,0.0138902090716542,0.00868792229378588,-0.00302924325427964,0.00158844631698647,0.007907944782137,0.0108432929135499,0.0120720381022763],"text":["Mcapitata_ATAC_TP7_1058","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP9_1306","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"color":"rgba(241,163,64,1)","size":11,"line":{"color":"rgba(241,163,64,1)"}},"type":"scatter","name":"Lilipuna.Fringe","textfont":{"color":"rgba(241,163,64,1)"},"error_y":{"color":"rgba(241,163,64,1)"},"error_x":{"color":"rgba(241,163,64,1)"},"line":{"color":"rgba(241,163,64,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0353283383245678,-0.0389494957653096,-0.0387156994742788,-0.0417337881359568,-0.0304936991285582,-0.0402493856659681,-0.0431807410559234,-0.0412446884480399,-0.0323077123296758,-0.0323721513550889,-0.0381188865102576,-0.0405122397483136,-0.0353802090917302,-0.0298078151036683,-0.0388331480283056,-0.034304902419877,-0.040824021862402,-0.0392288773251657,-0.0495706249695658,-0.03872879426252,-0.0386128238993237,-0.0245392315759141],"y":[0.0086215729500388,0.00777338915666925,0.0040734750542804,0.0109235427064767,0.000799730749840862,0.0117106593918568,0.012444694807492,0.0210155148565302,0.0127221641878477,0.00831680934603369,0.00892950422410608,0.000152537396615141,0.0082660743482269,0.0146984362206223,-0.0138922531290385,0.00655896064561861,0.00132292460985308,0.0162731219635841,0.0158125512291628,0.00783451889028733,0.0160110514242019,-0.00155570482177069],"text":["Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP3_1548","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP9_1467","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331"],"mode":"markers","marker":{"color":"rgba(254,224,182,1)","size":11,"line":{"color":"rgba(254,224,182,1)"}},"type":"scatter","name":"Reef.11.13","textfont":{"color":"rgba(254,224,182,1)"},"error_y":{"color":"rgba(254,224,182,1)"},"error_x":{"color":"rgba(254,224,182,1)"},"line":{"color":"rgba(254,224,182,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0292450117061691,-0.0446426063690368,-0.0333414861309239,-0.0438216733551513,-0.034731950916434,-0.0457030920108826,-0.0311424153153421,-0.036073349213048,-0.0370520983792864,-0.0312727343293337,-0.0473494065380055,-0.041164967807276,-0.0376882267183699,-0.0310620019777845,-0.0453449334759705,-0.0352108274111176,-0.0377806154522496,-0.0351359895664081,-0.0400879968629934,-0.0414191317913502],"y":[-0.0132485412069704,0.0194529170663445,0.00451228042831412,0.0192432959893718,0.00598768670966988,-0.000666676311828346,-0.00740119524011746,0.0019450778832346,0.00992521138326167,-0.0139800210928299,0.00568786750382434,0.0060824885581495,0.00991678192004881,0.00264683161542478,0.00374011002691165,0.00198074219032276,0.014039578065848,0.000500250439251407,-0.00181861048194446,0.00471743216220084],"text":["Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP6_2402","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP6_2555","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP7_2419"],"mode":"markers","marker":{"color":"rgba(216,218,235,1)","size":11,"line":{"color":"rgba(216,218,235,1)"}},"type":"scatter","name":"Reef.18","textfont":{"color":"rgba(216,218,235,1)"},"error_y":{"color":"rgba(216,218,235,1)"},"error_x":{"color":"rgba(216,218,235,1)"},"line":{"color":"rgba(216,218,235,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0259984393575133,-0.0410654394546692,-0.024818382982012,-0.0353034318908815,-0.0394449374274202,-0.0351041069986101,-0.0308549605672989,-0.0329544729411579,-0.0292574405201576,-0.0347188028316388,-0.0280177223372142,-0.0434408995316462,-0.032692920164491,-0.0400282426133999,-0.0345190745745488,-0.0311168811525216,-0.0452789332863361,-0.036626004446606,-0.0298324067335494,-0.0334591162056912,-0.0382554432390951,-0.0350690329703446,-0.0205173239207203,-0.0402284760493714,-0.0338310961146996,-0.0340804430072334,-0.0349431404266913,-0.0368868925529881],"y":[0.00752789801107825,-0.00371152102354799,-0.000697769477866483,-0.00186750481337353,-0.00152714811752618,0.0121572380563792,-0.00343873854075542,0.00167275413307834,-0.001080879835058,0.000261492778772618,-0.000653872480533765,0.0128910671476036,-0.00166332417244114,0.00252221481742771,-0.00176506228597524,0.0091681429556099,0.00684506168998994,0.00999648517240926,0.000673671608306822,0.00635245914950481,0.0147113976353928,0.00902944719518959,0.00673634758945498,0.00626057772248626,-0.0119778351543156,0.00510976465980565,0.0121706893305822,-0.000774896722339206],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP8_1260","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1722","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317"],"mode":"markers","marker":{"color":"rgba(153,142,195,1)","size":11,"line":{"color":"rgba(153,142,195,1)"}},"type":"scatter","name":"Reef.35.36","textfont":{"color":"rgba(153,142,195,1)"},"error_y":{"color":"rgba(153,142,195,1)"},"error_x":{"color":"rgba(153,142,195,1)"},"line":{"color":"rgba(153,142,195,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0293737093915074,-0.0364525034473797,-0.0339195879748338,-0.0247921975595541,-0.0281402178675904,-0.0269030375756818,-0.0326241022925872,-0.0370454974061555,-0.0384514501961568,-0.0295116367933613,-0.0398474851566358,-0.0259358590739191,-0.0309172669884343,-0.0319312162365303,-0.0448675828280251,-0.0281309662843891,-0.0339820363209966,-0.031368639082219,-0.0262667679351678,-0.0318843306819061,-0.038279637555649,-0.0317292179162127,-0.0317196716412144,-0.0416876422450806,-0.0298857743800826,-0.0340432555309663,-0.0316155425076258,-0.0268278657682002],"y":[0.00740082660305173,0.00371171966129259,0.00210165612485316,0.0123735078487645,-0.00680121588758422,0.00269733257968474,0.0124717396585872,0.00982664005711287,0.0140393765263796,0.00473295591764761,0.00947956277800932,-0.00285686968435638,0.0102220592693368,-0.0104454792593248,0.00479853733796544,-0.0117125282918918,0.00199651447173441,0.00658577111237137,-0.00420618113072376,0.0180438290013235,0.0130500696695487,0.0097603154605501,0.0229341965185012,0.0180986407976146,-0.00920344431416447,0.0106575189293047,0.00480025396540249,0.00133554972561204],"text":["Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP8_1779","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP8_1235","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP8_1246"],"mode":"markers","marker":{"color":"rgba(84,39,136,1)","size":11,"line":{"color":"rgba(84,39,136,1)"}},"type":"scatter","name":"Reef.42.43","textfont":{"color":"rgba(84,39,136,1)"},"error_y":{"color":"rgba(84,39,136,1)"},"error_x":{"color":"rgba(84,39,136,1)"},"line":{"color":"rgba(84,39,136,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.168994877569605,0.21509720342966,0.254568411526922,0.138856275324994,0.171702481377677,0.195859199686892,0.180539848936965,0.16424092655158,0.149246060042661,0.121203952823388,0.192417055323507,0.178056638666907,0.23000314893794,0.224335537594743,0.259111940141535,0.158023529046395,0.200859784931242,0.183815309324472],"y":[0.295345772294402,0.152155892960803,0.398763226164787,-0.0992143817114808,-0.167647580548066,0.0450558662621251,-0.0898255568135078,-0.102820139293431,-0.0547096829421381,-0.0800307901453591,-0.133405739687102,-0.128401681801111,0.366329662334079,0.158704545595299,0.404324290130019,-0.108716729952454,-0.184837002934645,0.0471436397166563],"text":["SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"color":"rgba(31,120,180,1)","size":11,"line":{"color":"rgba(31,120,180,1)"}},"type":"scatter","name":"Waiopae","textfont":{"color":"rgba(31,120,180,1)"},"error_y":{"color":"rgba(31,120,180,1)"},"error_x":{"color":"rgba(31,120,180,1)"},"line":{"color":"rgba(31,120,180,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
```

```r
# Ploidy
p <- plot_ly(df, x=~V1, y=~V2, text=rownames(df), type = "scatter",
             mode="markers", color=~ploidy, colors=label.colors, marker=list(size=11))
p <- layout(p,title="Without individual allele frequency - colored by Ploidy",
            xaxis=list(title="PC1"),
            yaxis=list(title="PC2"))
p
```

```{=html}
<div id="htmlwidget-dc3fd17b5c43cdc60655" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-dc3fd17b5c43cdc60655">{"x":{"visdat":{"1488446e2663d":["function () ","plotlyVisDat"]},"cur_data":"1488446e2663d","attrs":{"1488446e2663d":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009","SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747","SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#33a02c","#1f78b4","#01665e","#7570b3","#808080","#1f78b4","#33a02c","#543005","#dfc27d"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Ploidy","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0259984393575133,-0.0353283383245678,-0.0323376640238664,-0.0410654394546692,-0.024818382982012,-0.0235045784679117,-0.0292450117061691,-0.0293737093915074,-0.0334676804562726,-0.0364525034473797,-0.0353034318908815,-0.0446426063690368,-0.038112810285706,-0.0389494957653096,-0.0329575914573595,-0.0394449374274202,-0.0436701159998899,-0.036332483373355,-0.0333414861309239,-0.0339195879748338,-0.0247921975595541,-0.0281402178675904,-0.0269030375756818,-0.0438216733551513,-0.0306404187791171,-0.0351041069986101,-0.0308549605672989,-0.0329544729411579,-0.0292574405201576,-0.0326241022925872,-0.0447606336039882,-0.0347188028316388,-0.0280177223372142,-0.0387156994742788,-0.034731950916434,-0.0434408995316462,-0.032692920164491,-0.0400282426133999,-0.0417337881359568,-0.0370454974061555,-0.0304936991285582,-0.0345190745745488,-0.0400429811263431,-0.0384514501961568,-0.0457030920108826,-0.0311424153153421,-0.0295116367933613,-0.0377975735649241,-0.0311168811525216,-0.036073349213048,-0.0398474851566358,-0.0259358590739191,-0.0402493856659681,-0.0452789332863361,-0.0431807410559234,-0.036626004446606,-0.0370520983792864,-0.0309172669884343,-0.0353843096320885,-0.0382516942893998,-0.0298324067335494,-0.0319312162365303,-0.0448675828280251,-0.0281309662843891,-0.0376159759943272,-0.0412446884480399,-0.0399493733930138,-0.0323077123296758,-0.0346922932378603,-0.0347757908668086,-0.0323721513550889,-0.0312727343293337,-0.0429674931170672,-0.0473494065380055,-0.0104889376097138,-0.041164967807276,-0.0339820363209966,-0.031368639082219,-0.0334591162056912,-0.0262667679351678,-0.0318843306819061,-0.0382554432390951,-0.0350690329703446,-0.0424267530018416,-0.0381188865102576,-0.0304163527300572,-0.0348077090040011,-0.0205173239207203,-0.0376882267183699,-0.0402284760493714,-0.0405122397483136,-0.038279637555649,-0.0338310961146996,-0.0317292179162127,-0.0310620019777845,-0.0453449334759705,-0.0319617373876199,-0.0353802090917302,-0.0352108274111176,-0.0298078151036683,-0.0388331480283056,-0.034304902419877,-0.0334669177230852,-0.0455368930964151,-0.0417314134844289,-0.0377806154522496,-0.0340804430072334,-0.0311981113407001,-0.0417398345418761,-0.0355274227864432,-0.0317196716412144,-0.0327391327721698,-0.0374051471455735,-0.0351359895664081,-0.0416876422450806,-0.040824021862402,-0.0298857743800826,-0.0400879968629934,-0.0308552264537443,-0.0392288773251657,-0.0349431404266913,-0.0368868925529881,-0.0495706249695658,-0.0340432555309663,-0.0316155425076258,-0.0414191317913502,-0.03872879426252,-0.0268278657682002,-0.0453006255320472,-0.0386128238993237,-0.0245392315759141,-0.0390200030270299,0.101947155132633,0.169925374298519,0.151445130101844,0.118926735652195,0.140757365595538,0.123406412242653,0.100232185621671,0.151188360915327,0.155688253938549,0.168994877569605,0.21509720342966,0.254568411526922,0.138856275324994,0.171702481377677,0.195859199686892,0.180539848936965,0.16424092655158,0.149246060042661,0.121203952823388,0.192417055323507,0.178056638666907,0.23000314893794,0.224335537594743,0.259111940141535,0.158023529046395,0.200859784931242,0.183815309324472],"y":[0.00752789801107825,0.0086215729500388,0.0069021802979296,-0.00371152102354799,-0.000697769477866483,0.00402462237005923,-0.0132485412069704,0.00740082660305173,0.00882482895190268,0.00371171966129259,-0.00186750481337353,0.0194529170663445,0.0103064901567753,0.00777338915666925,0.0208961644592763,-0.00152714811752618,0.0113391097714042,0.00891407969436935,0.00451228042831412,0.00210165612485316,0.0123735078487645,-0.00680121588758422,0.00269733257968474,0.0192432959893718,0.00945045541406599,0.0121572380563792,-0.00343873854075542,0.00167275413307834,-0.001080879835058,0.0124717396585872,0.0122609393978024,0.000261492778772618,-0.000653872480533765,0.0040734750542804,0.00598768670966988,0.0128910671476036,-0.00166332417244114,0.00252221481742771,0.0109235427064767,0.00982664005711287,0.000799730749840862,-0.00176506228597524,0.00739500248127063,0.0140393765263796,-0.000666676311828346,-0.00740119524011746,0.00473295591764761,0.00185949292475981,0.0091681429556099,0.0019450778832346,0.00947956277800932,-0.00285686968435638,0.0117106593918568,0.00684506168998994,0.012444694807492,0.00999648517240926,0.00992521138326167,0.0102220592693368,0.00460787537043577,0.00816602244194009,0.000673671608306822,-0.0104454792593248,0.00479853733796544,-0.0117125282918918,0.00618920896918417,0.0210155148565302,0.0138902090716542,0.0127221641878477,0.011215012688699,-0.00978118197248325,0.00831680934603369,-0.0139800210928299,0.00868792229378588,0.00568786750382434,-0.0110622926807172,0.0060824885581495,0.00199651447173441,0.00658577111237137,0.00635245914950481,-0.00420618113072376,0.0180438290013235,0.0147113976353928,0.00902944719518959,0.00264298440373561,0.00892950422410608,-0.00302924325427964,-0.000496565143664624,0.00673634758945498,0.00991678192004881,0.00626057772248626,0.000152537396615141,0.0130500696695487,-0.0119778351543156,0.0097603154605501,0.00264683161542478,0.00374011002691165,0.00158844631698647,0.0082660743482269,0.00198074219032276,0.0146984362206223,-0.0138922531290385,0.00655896064561861,0.007907944782137,0.0105221320713179,0.00959425609874879,0.014039578065848,0.00510976465980565,-0.000369487216460847,0.0112162997085473,0.00185025766104319,0.0229341965185012,-0.00273864713828892,0.0110126232167368,0.000500250439251407,0.0180986407976146,0.00132292460985308,-0.00920344431416447,-0.00181861048194446,0.0108432929135499,0.0162731219635841,0.0121706893305822,-0.000774896722339206,0.0158125512291628,0.0106575189293047,0.00480025396540249,0.00471743216220084,0.00783451889028733,0.00133554972561204,0.00241266999768364,0.0160110514242019,-0.00155570482177069,0.0120720381022763,-0.130753842904154,-0.211321869131487,-0.196011935514122,-0.12345961880238,-0.138581199128536,-0.128053118706528,-0.128096276620108,-0.198250962083971,-0.197472331674013,0.295345772294402,0.152155892960803,0.398763226164787,-0.0992143817114808,-0.167647580548066,0.0450558662621251,-0.0898255568135078,-0.102820139293431,-0.0547096829421381,-0.0800307901453591,-0.133405739687102,-0.128401681801111,0.366329662334079,0.158704545595299,0.404324290130019,-0.108716729952454,-0.184837002934645,0.0471436397166563],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009","SRR5453739","SRR5453740","SRR5453741","SRR5453742","SRR5453743","SRR5453744","SRR5453745","SRR5453746","SRR5453747","SRR5453748","SRR5453749","SRR5453750","SRR5453751","SRR5453752","SRR5453753","SRR5453754","SRR5453755","SRR5453756","SRR5453757","SRR5453758","SRR5453759","SRR5453760","SRR5453761","SRR5453762","SRR5453763","SRR5453764","SRR5453765"],"mode":"markers","marker":{"colorbar":{"title":"ploidy","ticklen":2},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(214,189,189,1)"],["0.0833333333333333","rgba(244,201,153,1)"],["0.125","rgba(201,121,56,1)"],["0.166666666666667","rgba(201,130,88,1)"],["0.208333333333333","rgba(217,207,216,1)"],["0.25","rgba(152,126,185,1)"],["0.291666666666667","rgba(102,48,132,1)"],["0.333333333333333","rgba(195,119,97,1)"],["0.375","rgba(202,165,59,1)"],["0.416666666666667","rgba(97,162,47,1)"],["0.458333333333333","rgba(68,143,107,1)"],["0.5","rgba(31,120,180,1)"],["0.541666666666667","rgba(23,110,129,1)"],["0.583333333333333","rgba(38,104,108,1)"],["0.625","rgba(96,111,157,1)"],["0.666666666666667","rgba(122,117,162,1)"],["0.708333333333333","rgba(128,127,132,1)"],["0.75","rgba(97,124,154,1)"],["0.791666666666667","rgba(45,123,170,1)"],["0.833333333333333","rgba(67,146,96,1)"],["0.875","rgba(74,132,33,1)"],["0.916666666666667","rgba(86,67,11,1)"],["0.958333333333333","rgba(140,105,53,1)"],["1","rgba(223,194,125,1)"]],"showscale":false,"color":[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"size":11,"line":{"colorbar":{"title":"","ticklen":2},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(214,189,189,1)"],["0.0833333333333333","rgba(244,201,153,1)"],["0.125","rgba(201,121,56,1)"],["0.166666666666667","rgba(201,130,88,1)"],["0.208333333333333","rgba(217,207,216,1)"],["0.25","rgba(152,126,185,1)"],["0.291666666666667","rgba(102,48,132,1)"],["0.333333333333333","rgba(195,119,97,1)"],["0.375","rgba(202,165,59,1)"],["0.416666666666667","rgba(97,162,47,1)"],["0.458333333333333","rgba(68,143,107,1)"],["0.5","rgba(31,120,180,1)"],["0.541666666666667","rgba(23,110,129,1)"],["0.583333333333333","rgba(38,104,108,1)"],["0.625","rgba(96,111,157,1)"],["0.666666666666667","rgba(122,117,162,1)"],["0.708333333333333","rgba(128,127,132,1)"],["0.75","rgba(97,124,154,1)"],["0.791666666666667","rgba(45,123,170,1)"],["0.833333333333333","rgba(67,146,96,1)"],["0.875","rgba(74,132,33,1)"],["0.916666666666667","rgba(86,67,11,1)"],["0.958333333333333","rgba(140,105,53,1)"],["1","rgba(223,194,125,1)"]],"showscale":false,"color":[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}},"type":"scatter","xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0495706249695658,0.259111940141535],"y":[-0.211321869131487,0.404324290130019],"type":"scatter","mode":"markers","opacity":0,"hoverinfo":"none","showlegend":false,"marker":{"colorbar":{"title":"ploidy","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(214,189,189,1)"],["0.0833333333333333","rgba(244,201,153,1)"],["0.125","rgba(201,121,56,1)"],["0.166666666666667","rgba(201,130,88,1)"],["0.208333333333333","rgba(217,207,216,1)"],["0.25","rgba(152,126,185,1)"],["0.291666666666667","rgba(102,48,132,1)"],["0.333333333333333","rgba(195,119,97,1)"],["0.375","rgba(202,165,59,1)"],["0.416666666666667","rgba(97,162,47,1)"],["0.458333333333333","rgba(68,143,107,1)"],["0.5","rgba(31,120,180,1)"],["0.541666666666667","rgba(23,110,129,1)"],["0.583333333333333","rgba(38,104,108,1)"],["0.625","rgba(96,111,157,1)"],["0.666666666666667","rgba(122,117,162,1)"],["0.708333333333333","rgba(128,127,132,1)"],["0.75","rgba(97,124,154,1)"],["0.791666666666667","rgba(45,123,170,1)"],["0.833333333333333","rgba(67,146,96,1)"],["0.875","rgba(74,132,33,1)"],["0.916666666666667","rgba(86,67,11,1)"],["0.958333333333333","rgba(140,105,53,1)"],["1","rgba(223,194,125,1)"]],"showscale":true,"color":[2,4],"line":{"color":"rgba(255,127,14,1)"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
## [21] rmarkdown_2.16     labeling_0.4.2     stringr_1.4.1      htmlwidgets_1.5.4 
## [25] munsell_0.5.0      compiler_4.1.2     xfun_0.33          pkgconfig_2.0.3   
## [29] htmltools_0.5.3    tidyselect_1.1.2   gridExtra_2.3      fansi_1.0.3       
## [33] withr_2.5.0        MASS_7.3-58.1      bitops_1.0-7       grid_4.1.2        
## [37] nlme_3.1-159       jsonlite_1.8.0     gtable_0.3.1       lifecycle_1.0.2   
## [41] DBI_1.1.3          magrittr_2.0.3     scales_1.2.1       KernSmooth_2.23-20
## [45] cli_3.4.0          stringi_1.7.8      cachem_1.0.6       farver_2.1.1      
## [49] bslib_0.4.0        ellipsis_0.3.2     generics_0.1.3     vctrs_0.4.1       
## [53] tools_4.1.2        glue_1.6.2         purrr_0.3.4        crosstalk_1.2.0   
## [57] parallel_4.1.2     fastmap_1.1.0      yaml_2.3.5         colorspace_2.0-3  
## [61] caTools_1.18.2     sass_0.4.2
```
