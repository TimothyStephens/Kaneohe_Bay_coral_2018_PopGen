---
title: "Plot `vcf_clone_detect.py` results for *M. capitata* RNA-seq samples"
author: "Timothy Stephens"
date: "28/09/2022"
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
samples.info <- read.table("../../samples_Mcapitata.annotations.txt", header=T, comment.char='')
rownames(samples.info) <- samples.info$sample
samples.info
```

```
##                                            sample   species treatment timepoint
## Mcapitata_ATAC_TP11_2302 Mcapitata_ATAC_TP11_2302 Mcapitata      ATAC      TP11
## Mcapitata_HTAC_TP11_2380 Mcapitata_HTAC_TP11_2380 Mcapitata      HTAC      TP11
## Mcapitata_ATAC_TP6_2402   Mcapitata_ATAC_TP6_2402 Mcapitata      ATAC       TP6
## Mcapitata_ATAC_TP12_2403 Mcapitata_ATAC_TP12_2403 Mcapitata      ATAC      TP12
## Mcapitata_ATAC_TP1_1037   Mcapitata_ATAC_TP1_1037 Mcapitata      ATAC       TP1
## Mcapitata_ATAC_TP1_1600   Mcapitata_ATAC_TP1_1600 Mcapitata      ATAC       TP1
## Mcapitata_ATAC_TP1_1652   Mcapitata_ATAC_TP1_1652 Mcapitata      ATAC       TP1
## Mcapitata_ATAC_TP3_1101   Mcapitata_ATAC_TP3_1101 Mcapitata      ATAC       TP3
## Mcapitata_ATAC_TP3_1548   Mcapitata_ATAC_TP3_1548 Mcapitata      ATAC       TP3
## Mcapitata_ATAC_TP3_1628   Mcapitata_ATAC_TP3_1628 Mcapitata      ATAC       TP3
## Mcapitata_ATAC_TP4_1108   Mcapitata_ATAC_TP4_1108 Mcapitata      ATAC       TP4
## Mcapitata_ATAC_TP4_1609   Mcapitata_ATAC_TP4_1609 Mcapitata      ATAC       TP4
## Mcapitata_ATAC_TP4_1651   Mcapitata_ATAC_TP4_1651 Mcapitata      ATAC       TP4
## Mcapitata_ATAC_TP5_1196   Mcapitata_ATAC_TP5_1196 Mcapitata      ATAC       TP5
## Mcapitata_ATAC_TP5_1610   Mcapitata_ATAC_TP5_1610 Mcapitata      ATAC       TP5
## Mcapitata_ATAC_TP5_1776   Mcapitata_ATAC_TP5_1776 Mcapitata      ATAC       TP5
## Mcapitata_ATAC_TP6_1114   Mcapitata_ATAC_TP6_1114 Mcapitata      ATAC       TP6
## Mcapitata_ATAC_TP6_1611   Mcapitata_ATAC_TP6_1611 Mcapitata      ATAC       TP6
## Mcapitata_ATAC_TP7_1058   Mcapitata_ATAC_TP7_1058 Mcapitata      ATAC       TP7
## Mcapitata_ATAC_TP7_1455   Mcapitata_ATAC_TP7_1455 Mcapitata      ATAC       TP7
## Mcapitata_ATAC_TP7_1499   Mcapitata_ATAC_TP7_1499 Mcapitata      ATAC       TP7
## Mcapitata_ATAC_TP8_1083   Mcapitata_ATAC_TP8_1083 Mcapitata      ATAC       TP8
## Mcapitata_ATAC_TP8_1436   Mcapitata_ATAC_TP8_1436 Mcapitata      ATAC       TP8
## Mcapitata_ATAC_TP8_1779   Mcapitata_ATAC_TP8_1779 Mcapitata      ATAC       TP8
## Mcapitata_ATAC_TP9_1121   Mcapitata_ATAC_TP9_1121 Mcapitata      ATAC       TP9
## Mcapitata_ATAC_TP9_1420   Mcapitata_ATAC_TP9_1420 Mcapitata      ATAC       TP9
## Mcapitata_ATAC_TP9_1580   Mcapitata_ATAC_TP9_1580 Mcapitata      ATAC       TP9
## Mcapitata_ATAC_TP10_1095 Mcapitata_ATAC_TP10_1095 Mcapitata      ATAC      TP10
## Mcapitata_ATAC_TP10_1561 Mcapitata_ATAC_TP10_1561 Mcapitata      ATAC      TP10
## Mcapitata_ATAC_TP10_1631 Mcapitata_ATAC_TP10_1631 Mcapitata      ATAC      TP10
## Mcapitata_ATAC_TP11_1076 Mcapitata_ATAC_TP11_1076 Mcapitata      ATAC      TP11
## Mcapitata_ATAC_TP11_1644 Mcapitata_ATAC_TP11_1644 Mcapitata      ATAC      TP11
## Mcapitata_ATAC_TP12_1120 Mcapitata_ATAC_TP12_1120 Mcapitata      ATAC      TP12
## Mcapitata_ATAC_TP12_1452 Mcapitata_ATAC_TP12_1452 Mcapitata      ATAC      TP12
## Mcapitata_ATHC_TP1_1218   Mcapitata_ATHC_TP1_1218 Mcapitata      ATHC       TP1
## Mcapitata_ATHC_TP1_1826   Mcapitata_ATHC_TP1_1826 Mcapitata      ATHC       TP1
## Mcapitata_ATHC_TP1_2068   Mcapitata_ATHC_TP1_2068 Mcapitata      ATHC       TP1
## Mcapitata_ATHC_TP3_1544   Mcapitata_ATHC_TP3_1544 Mcapitata      ATHC       TP3
## Mcapitata_ATHC_TP3_2731   Mcapitata_ATHC_TP3_2731 Mcapitata      ATHC       TP3
## Mcapitata_ATHC_TP3_2866   Mcapitata_ATHC_TP3_2866 Mcapitata      ATHC       TP3
## Mcapitata_ATHC_TP4_1221   Mcapitata_ATHC_TP4_1221 Mcapitata      ATHC       TP4
## Mcapitata_ATHC_TP4_2561   Mcapitata_ATHC_TP4_2561 Mcapitata      ATHC       TP4
## Mcapitata_ATHC_TP4_2734   Mcapitata_ATHC_TP4_2734 Mcapitata      ATHC       TP4
## Mcapitata_ATHC_TP5_1229   Mcapitata_ATHC_TP5_1229 Mcapitata      ATHC       TP5
## Mcapitata_ATHC_TP5_1706   Mcapitata_ATHC_TP5_1706 Mcapitata      ATHC       TP5
## Mcapitata_ATHC_TP5_2986   Mcapitata_ATHC_TP5_2986 Mcapitata      ATHC       TP5
## Mcapitata_ATHC_TP6_1212   Mcapitata_ATHC_TP6_1212 Mcapitata      ATHC       TP6
## Mcapitata_ATHC_TP6_2016   Mcapitata_ATHC_TP6_2016 Mcapitata      ATHC       TP6
## Mcapitata_ATHC_TP6_2555   Mcapitata_ATHC_TP6_2555 Mcapitata      ATHC       TP6
## Mcapitata_ATHC_TP7_1223   Mcapitata_ATHC_TP7_1223 Mcapitata      ATHC       TP7
## Mcapitata_ATHC_TP7_2860   Mcapitata_ATHC_TP7_2860 Mcapitata      ATHC       TP7
## Mcapitata_ATHC_TP7_2875   Mcapitata_ATHC_TP7_2875 Mcapitata      ATHC       TP7
## Mcapitata_ATHC_TP8_1260   Mcapitata_ATHC_TP8_1260 Mcapitata      ATHC       TP8
## Mcapitata_ATHC_TP8_2735   Mcapitata_ATHC_TP8_2735 Mcapitata      ATHC       TP8
## Mcapitata_ATHC_TP8_2753   Mcapitata_ATHC_TP8_2753 Mcapitata      ATHC       TP8
## Mcapitata_ATHC_TP9_1148   Mcapitata_ATHC_TP9_1148 Mcapitata      ATHC       TP9
## Mcapitata_ATHC_TP9_2862   Mcapitata_ATHC_TP9_2862 Mcapitata      ATHC       TP9
## Mcapitata_ATHC_TP9_2995   Mcapitata_ATHC_TP9_2995 Mcapitata      ATHC       TP9
## Mcapitata_ATHC_TP10_1204 Mcapitata_ATHC_TP10_1204 Mcapitata      ATHC      TP10
## Mcapitata_ATHC_TP10_2554 Mcapitata_ATHC_TP10_2554 Mcapitata      ATHC      TP10
## Mcapitata_ATHC_TP10_2737 Mcapitata_ATHC_TP10_2737 Mcapitata      ATHC      TP10
## Mcapitata_ATHC_TP11_1237 Mcapitata_ATHC_TP11_1237 Mcapitata      ATHC      TP11
## Mcapitata_ATHC_TP11_2188 Mcapitata_ATHC_TP11_2188 Mcapitata      ATHC      TP11
## Mcapitata_ATHC_TP11_2756 Mcapitata_ATHC_TP11_2756 Mcapitata      ATHC      TP11
## Mcapitata_ATHC_TP12_1154 Mcapitata_ATHC_TP12_1154 Mcapitata      ATHC      TP12
## Mcapitata_ATHC_TP12_2736 Mcapitata_ATHC_TP12_2736 Mcapitata      ATHC      TP12
## Mcapitata_ATHC_TP12_2990 Mcapitata_ATHC_TP12_2990 Mcapitata      ATHC      TP12
## Mcapitata_HTAC_TP1_1579   Mcapitata_HTAC_TP1_1579 Mcapitata      HTAC       TP1
## Mcapitata_HTAC_TP1_2153   Mcapitata_HTAC_TP1_2153 Mcapitata      HTAC       TP1
## Mcapitata_HTAC_TP1_2183   Mcapitata_HTAC_TP1_2183 Mcapitata      HTAC       TP1
## Mcapitata_HTAC_TP3_1289   Mcapitata_HTAC_TP3_1289 Mcapitata      HTAC       TP3
## Mcapitata_HTAC_TP3_1751   Mcapitata_HTAC_TP3_1751 Mcapitata      HTAC       TP3
## Mcapitata_HTAC_TP3_2021   Mcapitata_HTAC_TP3_2021 Mcapitata      HTAC       TP3
## Mcapitata_HTAC_TP4_1269   Mcapitata_HTAC_TP4_1269 Mcapitata      HTAC       TP4
## Mcapitata_HTAC_TP4_1481   Mcapitata_HTAC_TP4_1481 Mcapitata      HTAC       TP4
## Mcapitata_HTAC_TP4_2000   Mcapitata_HTAC_TP4_2000 Mcapitata      HTAC       TP4
## Mcapitata_HTAC_TP5_1321   Mcapitata_HTAC_TP5_1321 Mcapitata      HTAC       TP5
## Mcapitata_HTAC_TP5_1583   Mcapitata_HTAC_TP5_1583 Mcapitata      HTAC       TP5
## Mcapitata_HTAC_TP5_1997   Mcapitata_HTAC_TP5_1997 Mcapitata      HTAC       TP5
## Mcapitata_HTAC_TP6_1496   Mcapitata_HTAC_TP6_1496 Mcapitata      HTAC       TP6
## Mcapitata_HTAC_TP6_1588   Mcapitata_HTAC_TP6_1588 Mcapitata      HTAC       TP6
## Mcapitata_HTAC_TP6_1705   Mcapitata_HTAC_TP6_1705 Mcapitata      HTAC       TP6
## Mcapitata_HTAC_TP7_1278   Mcapitata_HTAC_TP7_1278 Mcapitata      HTAC       TP7
## Mcapitata_HTAC_TP7_1645   Mcapitata_HTAC_TP7_1645 Mcapitata      HTAC       TP7
## Mcapitata_HTAC_TP7_1722   Mcapitata_HTAC_TP7_1722 Mcapitata      HTAC       TP7
## Mcapitata_HTAC_TP8_1235   Mcapitata_HTAC_TP8_1235 Mcapitata      HTAC       TP8
## Mcapitata_HTAC_TP8_2386   Mcapitata_HTAC_TP8_2386 Mcapitata      HTAC       TP8
## Mcapitata_HTAC_TP8_2410   Mcapitata_HTAC_TP8_2410 Mcapitata      HTAC       TP8
## Mcapitata_HTAC_TP9_1306   Mcapitata_HTAC_TP9_1306 Mcapitata      HTAC       TP9
## Mcapitata_HTAC_TP9_1467   Mcapitata_HTAC_TP9_1467 Mcapitata      HTAC       TP9
## Mcapitata_HTAC_TP9_2412   Mcapitata_HTAC_TP9_2412 Mcapitata      HTAC       TP9
## Mcapitata_HTAC_TP10_1315 Mcapitata_HTAC_TP10_1315 Mcapitata      HTAC      TP10
## Mcapitata_HTAC_TP10_1478 Mcapitata_HTAC_TP10_1478 Mcapitata      HTAC      TP10
## Mcapitata_HTAC_TP10_1754 Mcapitata_HTAC_TP10_1754 Mcapitata      HTAC      TP10
## Mcapitata_HTAC_TP11_1248 Mcapitata_HTAC_TP11_1248 Mcapitata      HTAC      TP11
## Mcapitata_HTAC_TP11_1562 Mcapitata_HTAC_TP11_1562 Mcapitata      HTAC      TP11
## Mcapitata_HTAC_TP12_1632 Mcapitata_HTAC_TP12_1632 Mcapitata      HTAC      TP12
## Mcapitata_HTAC_TP12_1729 Mcapitata_HTAC_TP12_1729 Mcapitata      HTAC      TP12
## Mcapitata_HTAC_TP12_2007 Mcapitata_HTAC_TP12_2007 Mcapitata      HTAC      TP12
## Mcapitata_HTHC_TP1_1145   Mcapitata_HTHC_TP1_1145 Mcapitata      HTHC       TP1
## Mcapitata_HTHC_TP1_1323   Mcapitata_HTHC_TP1_1323 Mcapitata      HTHC       TP1
## Mcapitata_HTHC_TP1_2081   Mcapitata_HTHC_TP1_2081 Mcapitata      HTHC       TP1
## Mcapitata_HTHC_TP3_1128   Mcapitata_HTHC_TP3_1128 Mcapitata      HTHC       TP3
## Mcapitata_HTHC_TP3_1277   Mcapitata_HTHC_TP3_1277 Mcapitata      HTHC       TP3
## Mcapitata_HTHC_TP3_2518   Mcapitata_HTHC_TP3_2518 Mcapitata      HTHC       TP3
## Mcapitata_HTHC_TP4_1124   Mcapitata_HTHC_TP4_1124 Mcapitata      HTHC       TP4
## Mcapitata_HTHC_TP4_1328   Mcapitata_HTHC_TP4_1328 Mcapitata      HTHC       TP4
## Mcapitata_HTHC_TP4_2204   Mcapitata_HTHC_TP4_2204 Mcapitata      HTHC       TP4
## Mcapitata_HTHC_TP5_1345   Mcapitata_HTHC_TP5_1345 Mcapitata      HTHC       TP5
## Mcapitata_HTHC_TP5_1449   Mcapitata_HTHC_TP5_1449 Mcapitata      HTHC       TP5
## Mcapitata_HTHC_TP5_1694   Mcapitata_HTHC_TP5_1694 Mcapitata      HTHC       TP5
## Mcapitata_HTHC_TP6_1164   Mcapitata_HTHC_TP6_1164 Mcapitata      HTHC       TP6
## Mcapitata_HTHC_TP6_1317   Mcapitata_HTHC_TP6_1317 Mcapitata      HTHC       TP6
## Mcapitata_HTHC_TP6_1604   Mcapitata_HTHC_TP6_1604 Mcapitata      HTHC       TP6
## Mcapitata_HTHC_TP7_1126   Mcapitata_HTHC_TP7_1126 Mcapitata      HTHC       TP7
## Mcapitata_HTHC_TP7_1250   Mcapitata_HTHC_TP7_1250 Mcapitata      HTHC       TP7
## Mcapitata_HTHC_TP7_2419   Mcapitata_HTHC_TP7_2419 Mcapitata      HTHC       TP7
## Mcapitata_HTHC_TP8_1082   Mcapitata_HTHC_TP8_1082 Mcapitata      HTHC       TP8
## Mcapitata_HTHC_TP8_1246   Mcapitata_HTHC_TP8_1246 Mcapitata      HTHC       TP8
## Mcapitata_HTHC_TP8_2067   Mcapitata_HTHC_TP8_2067 Mcapitata      HTHC       TP8
## Mcapitata_HTHC_TP9_1078   Mcapitata_HTHC_TP9_1078 Mcapitata      HTHC       TP9
## Mcapitata_HTHC_TP9_1331   Mcapitata_HTHC_TP9_1331 Mcapitata      HTHC       TP9
## Mcapitata_HTHC_TP9_2009   Mcapitata_HTHC_TP9_2009 Mcapitata      HTHC       TP9
## Mcapitata_HTHC_TP10_1074 Mcapitata_HTHC_TP10_1074 Mcapitata      HTHC      TP10
## Mcapitata_HTHC_TP10_1332 Mcapitata_HTHC_TP10_1332 Mcapitata      HTHC      TP10
## Mcapitata_HTHC_TP10_1689 Mcapitata_HTHC_TP10_1689 Mcapitata      HTHC      TP10
## Mcapitata_HTHC_TP11_1178 Mcapitata_HTHC_TP11_1178 Mcapitata      HTHC      TP11
## Mcapitata_HTHC_TP11_1270 Mcapitata_HTHC_TP11_1270 Mcapitata      HTHC      TP11
## Mcapitata_HTHC_TP11_2511 Mcapitata_HTHC_TP11_2511 Mcapitata      HTHC      TP11
## Mcapitata_HTHC_TP12_1140 Mcapitata_HTHC_TP12_1140 Mcapitata      HTHC      TP12
## Mcapitata_HTHC_TP12_1274 Mcapitata_HTHC_TP12_1274 Mcapitata      HTHC      TP12
## Mcapitata_HTHC_TP12_2190 Mcapitata_HTHC_TP12_2190 Mcapitata      HTHC      TP12
##                          plugid            reef reef_color ploidy ploidy_color
## Mcapitata_ATAC_TP11_2302   2302         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP11_2380   2380         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATAC_TP6_2402    2402         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATAC_TP12_2403   2403         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATAC_TP1_1037    1037      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP1_1600    1600      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP1_1652    1652            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP3_1101    1101            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP3_1548    1548      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATAC_TP3_1628    1628            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP4_1108    1108      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP4_1609    1609            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP4_1651    1651            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP5_1196    1196         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATAC_TP5_1610    1610      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP5_1776    1776      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP6_1114    1114      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP6_1611    1611      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP7_1058    1058 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_ATAC_TP7_1455    1455      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP7_1499    1499      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP8_1083    1083      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP8_1436    1436      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP8_1779    1779      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP9_1121    1121            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP9_1420    1420      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP9_1580    1580      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP10_1095   1095      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP10_1561   1561      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATAC_TP10_1631   1631            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP11_1076   1076      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP11_1644   1644            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP12_1120   1120      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP12_1452   1452      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP1_1218    1218      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP1_1826    1826      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP1_2068    2068      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP3_1544    1544         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATHC_TP3_2731    2731      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP3_2866    2866            HIMB    #b35806      2      #01665e
## Mcapitata_ATHC_TP4_1221    1221      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP4_2561    2561         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATHC_TP4_2734    2734      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP5_1229    1229      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP5_1706    1706      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP5_2986    2986      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP6_1212    1212      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP6_2016    2016      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP6_2555    2555         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATHC_TP7_1223    1223      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP7_2860    2860            HIMB    #b35806      2      #01665e
## Mcapitata_ATHC_TP7_2875    2875            HIMB    #b35806      2      #01665e
## Mcapitata_ATHC_TP8_1260    1260      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP8_2735    2735      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP8_2753    2753      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP9_1148    1148      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP9_2862    2862            HIMB    #b35806      2      #01665e
## Mcapitata_ATHC_TP9_2995    2995      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP10_1204   1204      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP10_2554   2554         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATHC_TP10_2737   2737      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP11_1237   1237      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP11_2188   2188      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP11_2756   2756      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP12_1154   1154            HIMB    #b35806      2      #01665e
## Mcapitata_ATHC_TP12_2736   2736      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP12_2990   2990         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP1_1579    1579 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTAC_TP1_2153    2153         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP1_2183    2183      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTAC_TP3_1289    1289      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTAC_TP3_1751    1751      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTAC_TP3_2021    2021      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTAC_TP4_1269    1269      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTAC_TP4_1481    1481      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTAC_TP4_2000    2000            HIMB    #b35806      2      #01665e
## Mcapitata_HTAC_TP5_1321    1321      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTAC_TP5_1583    1583 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTAC_TP5_1997    1997            HIMB    #b35806      2      #01665e
## Mcapitata_HTAC_TP6_1496    1496      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTAC_TP6_1588    1588         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP6_1705    1705      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTAC_TP7_1278    1278      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTAC_TP7_1645    1645      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTAC_TP7_1722    1722      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTAC_TP8_1235    1235      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTAC_TP8_2386    2386         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP8_2410    2410         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP9_1306    1306 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTAC_TP9_1467    1467      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTAC_TP9_2412    2412         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP10_1315   1315 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTAC_TP10_1478   1478      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTAC_TP10_1754   1754            HIMB    #b35806      2      #01665e
## Mcapitata_HTAC_TP11_1248   1248            HIMB    #b35806      2      #01665e
## Mcapitata_HTAC_TP11_1562   1562      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTAC_TP12_1632   1632            HIMB    #b35806      4      #7570b3
## Mcapitata_HTAC_TP12_1729   1729         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP12_2007   2007      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTHC_TP1_1145    1145            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP1_1323    1323      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTHC_TP1_2081    2081            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP3_1128    1128            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP3_1277    1277            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP3_2518    2518         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTHC_TP4_1124    1124      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTHC_TP4_1328    1328      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP4_2204    2204      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTHC_TP5_1345    1345         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTHC_TP5_1449    1449 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTHC_TP5_1694    1694      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP6_1164    1164      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTHC_TP6_1317    1317      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTHC_TP6_1604    1604      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP7_1126    1126      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTHC_TP7_1250    1250      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTHC_TP7_2419    2419         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTHC_TP8_1082    1082      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP8_1246    1246      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTHC_TP8_2067    2067            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP9_1078    1078      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP9_1331    1331      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP9_2009    2009 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTHC_TP10_1074   1074      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP10_1332   1332      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP10_1689   1689      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP11_1178   1178 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTHC_TP11_1270   1270            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP11_2511   2511         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTHC_TP12_1140   1140            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP12_1274   1274            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP12_2190   2190      Reef.42.43    #542788      2      #01665e
##                              group group_color
## Mcapitata_ATAC_TP11_2302    Group1     #1f78b4
## Mcapitata_HTAC_TP11_2380    Group1     #1f78b4
## Mcapitata_ATAC_TP6_2402     Group2     #33a02c
## Mcapitata_ATAC_TP12_2403    Group2     #33a02c
## Mcapitata_ATAC_TP1_1037  Ungrouped     #808080
## Mcapitata_ATAC_TP1_1600  Ungrouped     #808080
## Mcapitata_ATAC_TP1_1652  Ungrouped     #808080
## Mcapitata_ATAC_TP3_1101  Ungrouped     #808080
## Mcapitata_ATAC_TP3_1548  Ungrouped     #808080
## Mcapitata_ATAC_TP3_1628  Ungrouped     #808080
## Mcapitata_ATAC_TP4_1108  Ungrouped     #808080
## Mcapitata_ATAC_TP4_1609  Ungrouped     #808080
## Mcapitata_ATAC_TP4_1651  Ungrouped     #808080
## Mcapitata_ATAC_TP5_1196  Ungrouped     #808080
## Mcapitata_ATAC_TP5_1610  Ungrouped     #808080
## Mcapitata_ATAC_TP5_1776  Ungrouped     #808080
## Mcapitata_ATAC_TP6_1114  Ungrouped     #808080
## Mcapitata_ATAC_TP6_1611  Ungrouped     #808080
## Mcapitata_ATAC_TP7_1058  Ungrouped     #808080
## Mcapitata_ATAC_TP7_1455  Ungrouped     #808080
## Mcapitata_ATAC_TP7_1499  Ungrouped     #808080
## Mcapitata_ATAC_TP8_1083  Ungrouped     #808080
## Mcapitata_ATAC_TP8_1436  Ungrouped     #808080
## Mcapitata_ATAC_TP8_1779  Ungrouped     #808080
## Mcapitata_ATAC_TP9_1121  Ungrouped     #808080
## Mcapitata_ATAC_TP9_1420  Ungrouped     #808080
## Mcapitata_ATAC_TP9_1580  Ungrouped     #808080
## Mcapitata_ATAC_TP10_1095 Ungrouped     #808080
## Mcapitata_ATAC_TP10_1561 Ungrouped     #808080
## Mcapitata_ATAC_TP10_1631 Ungrouped     #808080
## Mcapitata_ATAC_TP11_1076 Ungrouped     #808080
## Mcapitata_ATAC_TP11_1644 Ungrouped     #808080
## Mcapitata_ATAC_TP12_1120 Ungrouped     #808080
## Mcapitata_ATAC_TP12_1452 Ungrouped     #808080
## Mcapitata_ATHC_TP1_1218  Ungrouped     #808080
## Mcapitata_ATHC_TP1_1826  Ungrouped     #808080
## Mcapitata_ATHC_TP1_2068  Ungrouped     #808080
## Mcapitata_ATHC_TP3_1544  Ungrouped     #808080
## Mcapitata_ATHC_TP3_2731  Ungrouped     #808080
## Mcapitata_ATHC_TP3_2866  Ungrouped     #808080
## Mcapitata_ATHC_TP4_1221  Ungrouped     #808080
## Mcapitata_ATHC_TP4_2561  Ungrouped     #808080
## Mcapitata_ATHC_TP4_2734  Ungrouped     #808080
## Mcapitata_ATHC_TP5_1229  Ungrouped     #808080
## Mcapitata_ATHC_TP5_1706  Ungrouped     #808080
## Mcapitata_ATHC_TP5_2986  Ungrouped     #808080
## Mcapitata_ATHC_TP6_1212  Ungrouped     #808080
## Mcapitata_ATHC_TP6_2016  Ungrouped     #808080
## Mcapitata_ATHC_TP6_2555  Ungrouped     #808080
## Mcapitata_ATHC_TP7_1223  Ungrouped     #808080
## Mcapitata_ATHC_TP7_2860  Ungrouped     #808080
## Mcapitata_ATHC_TP7_2875  Ungrouped     #808080
## Mcapitata_ATHC_TP8_1260  Ungrouped     #808080
## Mcapitata_ATHC_TP8_2735  Ungrouped     #808080
## Mcapitata_ATHC_TP8_2753  Ungrouped     #808080
## Mcapitata_ATHC_TP9_1148  Ungrouped     #808080
## Mcapitata_ATHC_TP9_2862  Ungrouped     #808080
## Mcapitata_ATHC_TP9_2995  Ungrouped     #808080
## Mcapitata_ATHC_TP10_1204 Ungrouped     #808080
## Mcapitata_ATHC_TP10_2554 Ungrouped     #808080
## Mcapitata_ATHC_TP10_2737 Ungrouped     #808080
## Mcapitata_ATHC_TP11_1237 Ungrouped     #808080
## Mcapitata_ATHC_TP11_2188 Ungrouped     #808080
## Mcapitata_ATHC_TP11_2756 Ungrouped     #808080
## Mcapitata_ATHC_TP12_1154 Ungrouped     #808080
## Mcapitata_ATHC_TP12_2736 Ungrouped     #808080
## Mcapitata_ATHC_TP12_2990 Ungrouped     #808080
## Mcapitata_HTAC_TP1_1579  Ungrouped     #808080
## Mcapitata_HTAC_TP1_2153  Ungrouped     #808080
## Mcapitata_HTAC_TP1_2183  Ungrouped     #808080
## Mcapitata_HTAC_TP3_1289  Ungrouped     #808080
## Mcapitata_HTAC_TP3_1751  Ungrouped     #808080
## Mcapitata_HTAC_TP3_2021  Ungrouped     #808080
## Mcapitata_HTAC_TP4_1269  Ungrouped     #808080
## Mcapitata_HTAC_TP4_1481  Ungrouped     #808080
## Mcapitata_HTAC_TP4_2000  Ungrouped     #808080
## Mcapitata_HTAC_TP5_1321  Ungrouped     #808080
## Mcapitata_HTAC_TP5_1583  Ungrouped     #808080
## Mcapitata_HTAC_TP5_1997  Ungrouped     #808080
## Mcapitata_HTAC_TP6_1496  Ungrouped     #808080
## Mcapitata_HTAC_TP6_1588  Ungrouped     #808080
## Mcapitata_HTAC_TP6_1705  Ungrouped     #808080
## Mcapitata_HTAC_TP7_1278  Ungrouped     #808080
## Mcapitata_HTAC_TP7_1645  Ungrouped     #808080
## Mcapitata_HTAC_TP7_1722  Ungrouped     #808080
## Mcapitata_HTAC_TP8_1235  Ungrouped     #808080
## Mcapitata_HTAC_TP8_2386  Ungrouped     #808080
## Mcapitata_HTAC_TP8_2410  Ungrouped     #808080
## Mcapitata_HTAC_TP9_1306  Ungrouped     #808080
## Mcapitata_HTAC_TP9_1467  Ungrouped     #808080
## Mcapitata_HTAC_TP9_2412  Ungrouped     #808080
## Mcapitata_HTAC_TP10_1315 Ungrouped     #808080
## Mcapitata_HTAC_TP10_1478 Ungrouped     #808080
## Mcapitata_HTAC_TP10_1754 Ungrouped     #808080
## Mcapitata_HTAC_TP11_1248 Ungrouped     #808080
## Mcapitata_HTAC_TP11_1562 Ungrouped     #808080
## Mcapitata_HTAC_TP12_1632 Ungrouped     #808080
## Mcapitata_HTAC_TP12_1729 Ungrouped     #808080
## Mcapitata_HTAC_TP12_2007 Ungrouped     #808080
## Mcapitata_HTHC_TP1_1145  Ungrouped     #808080
## Mcapitata_HTHC_TP1_1323  Ungrouped     #808080
## Mcapitata_HTHC_TP1_2081  Ungrouped     #808080
## Mcapitata_HTHC_TP3_1128  Ungrouped     #808080
## Mcapitata_HTHC_TP3_1277  Ungrouped     #808080
## Mcapitata_HTHC_TP3_2518  Ungrouped     #808080
## Mcapitata_HTHC_TP4_1124  Ungrouped     #808080
## Mcapitata_HTHC_TP4_1328  Ungrouped     #808080
## Mcapitata_HTHC_TP4_2204  Ungrouped     #808080
## Mcapitata_HTHC_TP5_1345  Ungrouped     #808080
## Mcapitata_HTHC_TP5_1449  Ungrouped     #808080
## Mcapitata_HTHC_TP5_1694  Ungrouped     #808080
## Mcapitata_HTHC_TP6_1164  Ungrouped     #808080
## Mcapitata_HTHC_TP6_1317  Ungrouped     #808080
## Mcapitata_HTHC_TP6_1604  Ungrouped     #808080
## Mcapitata_HTHC_TP7_1126  Ungrouped     #808080
## Mcapitata_HTHC_TP7_1250  Ungrouped     #808080
## Mcapitata_HTHC_TP7_2419  Ungrouped     #808080
## Mcapitata_HTHC_TP8_1082  Ungrouped     #808080
## Mcapitata_HTHC_TP8_1246  Ungrouped     #808080
## Mcapitata_HTHC_TP8_2067  Ungrouped     #808080
## Mcapitata_HTHC_TP9_1078  Ungrouped     #808080
## Mcapitata_HTHC_TP9_1331  Ungrouped     #808080
## Mcapitata_HTHC_TP9_2009  Ungrouped     #808080
## Mcapitata_HTHC_TP10_1074 Ungrouped     #808080
## Mcapitata_HTHC_TP10_1332 Ungrouped     #808080
## Mcapitata_HTHC_TP10_1689 Ungrouped     #808080
## Mcapitata_HTHC_TP11_1178 Ungrouped     #808080
## Mcapitata_HTHC_TP11_1270 Ungrouped     #808080
## Mcapitata_HTHC_TP11_2511 Ungrouped     #808080
## Mcapitata_HTHC_TP12_1140 Ungrouped     #808080
## Mcapitata_HTHC_TP12_1274 Ungrouped     #808080
## Mcapitata_HTHC_TP12_2190 Ungrouped     #808080
```





# Load pairwise similarity scores (simScore) and cluster samples

Load pairwise percent shared SNP results and convert if from one line per pair format into a symetric matrix uisng the `xtabs` function.

```r
pairwise_percent_shared <- read.table("GVCFall.filtered.recode.vcf.gz.allelic_similarity.full.tsv", sep='\t', header=T)
pairwise_percent_shared.matrix <- xtabs(match_perc ~ ind1 + ind2, data=pairwise_percent_shared)
```



Process the simScore data ready for plotting (i.e., cluster samples based on computed distance and construct a dendrogram of the relationship).

```r
# Make a copy of input matrix
tmp.mat <- pairwise_percent_shared.matrix

# Generate clusters from the allelic similarity scores matrix
tmp.mat.dist   <- dist(tmp.mat, method="euclidean")
tmp.mat.hclust <- hclust(tmp.mat.dist, method="complete")

# Get sample colors for sidebar in the same order as the matrix
Rcols <- samples.info[rownames(tmp.mat),]$group_color
Ccols <- samples.info[colnames(tmp.mat),]$group_color

# Get min and max values of matrix
tmp.mat.min <- min(tmp.mat); tmp.mat.min
```

```
## [1] 79.16
```

```r
tmp.mat.max <- max(tmp.mat); tmp.mat.max
```

```
## [1] 100
```



Plot clustering results and write dendrogram to file.

```r
simScore.dendro <- as.dendrogram(tmp.mat.hclust)
write.dendrogram(simScore.dendro, "cluster_dendrogram.tre", edges=TRUE)
plot(tmp.mat.hclust)
```

![](plot_vcf_clone_detect_results_files/figure-html/plot_clusters-1.png)<!-- -->



Write ordered matrix to file incase we want it later.

```r
tmp.mat.order <- dendro_data(simScore.dendro, type="rectangle")[["labels"]][["label"]]
tmp.mat.ordered <- tmp.mat[
  order(match(rownames(tmp.mat), rev(tmp.mat.order))),
  order(match(colnames(tmp.mat), tmp.mat.order))
  ]
write.table(as.matrix(tmp.mat.ordered), 
            "GVCFall.filtered.recode.vcf.gz.allelic_similarity.full.matrix.tsv", 
            sep='\t', quote=FALSE, 
            row.names=TRUE, col.names=TRUE)
```





# Plot pairwise similarity scores

Plot simScores heatmap.

```r
par(cex.main=0.7, cex.lab=0.9, cex.axis=0.9) # Change the font size of the legend title, labels and axis
heatmap.2(x=tmp.mat,
          Rowv = simScore.dendro,
          Colv = simScore.dendro,
          col=cividis(length(seq(75, 100, 0.1))-1),
          breaks=seq(75, 100, 0.1),
          trace="none", # Dont draw "trace" line
          key=TRUE, # Show color-key
          margins=c(8,8), # Numeric vector of length 2 containing the margins for column and row names, respectively.
          ColSideColors=Ccols,
          RowSideColors=Rcols,
          offsetRow=0,
          offsetCol=0,
          cexRow=cexSize, # Positive numbers, used as cex.axis in for the row or column axis labeling.
          cexCol=cexSize,
          cellnote=formatC(tmp.mat, big.mark=","),
          notecex=0.1,
          )
```

![](plot_vcf_clone_detect_results_files/figure-html/plot_simScore_heatmap-1.png)<!-- -->



Plot simScores heatmap without cell notes/values.

```r
par(cex.main=0.7, cex.lab=0.9, cex.axis=0.9) # Change the font size of the legend title, labels and axis
heatmap.2(x=tmp.mat,
          Rowv = simScore.dendro,
          Colv = simScore.dendro,
          col=cividis(length(seq(75, 100, 0.1))-1),
          breaks=seq(75, 100, 0.1),
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
          )
```

![](plot_vcf_clone_detect_results_files/figure-html/plot_simScore_heatmap_noCellNotes-1.png)<!-- -->



Plot simScores heatmap with specific colors applied to selected ranges of values (i.e., value ranges picked from the histogram on the top left).

```r
par(cex.main=0.7, cex.lab=0.9, cex.axis=0.9) # Change the font size of the legend title, labels and axis
heatmap.2(x=tmp.mat,
          Rowv = simScore.dendro,
          Colv = simScore.dendro,
          col=c(rep(cividis(5)[1], length(seq(75.1, 80.6,  0.1))),
                rep(cividis(5)[2], length(seq(80.7, 84.3,  0.1))),
                rep(cividis(5)[3], length(seq(84.4, 91.0,  0.1))),
                rep(cividis(5)[4], length(seq(91.1, 94.0,  0.1))),
                rep(cividis(5)[5], length(seq(94.1, 100.0, 0.1)))
                ),
          breaks=seq(75, 100, 0.1),
          trace="none", # Don't draw "trace" line
          key=TRUE, # Show color-key
          margins=c(8,8), # Numeric vector of length 2 containing the margins for column and row names, respectively.
          ColSideColors=Ccols,
          RowSideColors=Rcols,
          offsetRow=0,
          offsetCol=0,
          cexRow=cexSize, # Positive numbers, used as cex.axis in for the row or column axis labeling.
          cexCol=cexSize,
          cellnote=formatC(tmp.mat, big.mark=","),
          notecex=0.1,
          )
```

![](plot_vcf_clone_detect_results_files/figure-html/plot_simScore_heatmap_withCutoff-1.png)<!-- -->



Plot simScores heatmap with specific colors applied to selected ranges of values (i.e., value ranges picked from the histogram on the top left) without cell notes/values.

```r
par(cex.main=0.7, cex.lab=0.9, cex.axis=0.9) # Change the font size of the legend title, labels and axis
heatmap.2(x=tmp.mat,
          Rowv = simScore.dendro,
          Colv = simScore.dendro,
          col=c(rep(cividis(5)[1], length(seq(75.1, 80.6,  0.1))),
                rep(cividis(5)[2], length(seq(80.7, 84.3,  0.1))),
                rep(cividis(5)[3], length(seq(84.4, 91.0,  0.1))),
                rep(cividis(5)[4], length(seq(91.1, 94.0,  0.1))),
                rep(cividis(5)[5], length(seq(94.1, 100.0, 0.1)))
                ),
          breaks=seq(75, 100, 0.1),
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
          )
```

![](plot_vcf_clone_detect_results_files/figure-html/plot_simScore_heatmap_withCutoff_noCellNotes-1.png)<!-- -->





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
