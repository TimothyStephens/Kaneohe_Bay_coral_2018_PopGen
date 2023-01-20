---
title: "Plot `ANGSD` results for *M. capitata* RNA-seq samples"
author: "Timothy Stephens"
date: "20/09/2022"
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
samples.info <- read.table("../../samples_Mcapitata.annotations.txt", header=T, comment.char='')
rownames(samples.info) <- samples.info$sample

samples.order <- read.table("bam.filelist.labels", header=F)
samples.info <- samples.info[samples.order$V1,]
samples.info
```

```
##                                            sample   species treatment timepoint
## Mcapitata_ATAC_TP10_1095 Mcapitata_ATAC_TP10_1095 Mcapitata      ATAC      TP10
## Mcapitata_ATAC_TP10_1561 Mcapitata_ATAC_TP10_1561 Mcapitata      ATAC      TP10
## Mcapitata_ATAC_TP10_1631 Mcapitata_ATAC_TP10_1631 Mcapitata      ATAC      TP10
## Mcapitata_ATAC_TP1_1037   Mcapitata_ATAC_TP1_1037 Mcapitata      ATAC       TP1
## Mcapitata_ATAC_TP11_1076 Mcapitata_ATAC_TP11_1076 Mcapitata      ATAC      TP11
## Mcapitata_ATAC_TP11_1644 Mcapitata_ATAC_TP11_1644 Mcapitata      ATAC      TP11
## Mcapitata_ATAC_TP11_2302 Mcapitata_ATAC_TP11_2302 Mcapitata      ATAC      TP11
## Mcapitata_ATAC_TP1_1600   Mcapitata_ATAC_TP1_1600 Mcapitata      ATAC       TP1
## Mcapitata_ATAC_TP1_1652   Mcapitata_ATAC_TP1_1652 Mcapitata      ATAC       TP1
## Mcapitata_ATAC_TP12_1120 Mcapitata_ATAC_TP12_1120 Mcapitata      ATAC      TP12
## Mcapitata_ATAC_TP12_1452 Mcapitata_ATAC_TP12_1452 Mcapitata      ATAC      TP12
## Mcapitata_ATAC_TP12_2403 Mcapitata_ATAC_TP12_2403 Mcapitata      ATAC      TP12
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
## Mcapitata_ATAC_TP6_2402   Mcapitata_ATAC_TP6_2402 Mcapitata      ATAC       TP6
## Mcapitata_ATAC_TP7_1058   Mcapitata_ATAC_TP7_1058 Mcapitata      ATAC       TP7
## Mcapitata_ATAC_TP7_1455   Mcapitata_ATAC_TP7_1455 Mcapitata      ATAC       TP7
## Mcapitata_ATAC_TP7_1499   Mcapitata_ATAC_TP7_1499 Mcapitata      ATAC       TP7
## Mcapitata_ATAC_TP8_1083   Mcapitata_ATAC_TP8_1083 Mcapitata      ATAC       TP8
## Mcapitata_ATAC_TP8_1436   Mcapitata_ATAC_TP8_1436 Mcapitata      ATAC       TP8
## Mcapitata_ATAC_TP8_1779   Mcapitata_ATAC_TP8_1779 Mcapitata      ATAC       TP8
## Mcapitata_ATAC_TP9_1121   Mcapitata_ATAC_TP9_1121 Mcapitata      ATAC       TP9
## Mcapitata_ATAC_TP9_1420   Mcapitata_ATAC_TP9_1420 Mcapitata      ATAC       TP9
## Mcapitata_ATAC_TP9_1580   Mcapitata_ATAC_TP9_1580 Mcapitata      ATAC       TP9
## Mcapitata_ATHC_TP10_1204 Mcapitata_ATHC_TP10_1204 Mcapitata      ATHC      TP10
## Mcapitata_ATHC_TP10_2554 Mcapitata_ATHC_TP10_2554 Mcapitata      ATHC      TP10
## Mcapitata_ATHC_TP10_2737 Mcapitata_ATHC_TP10_2737 Mcapitata      ATHC      TP10
## Mcapitata_ATHC_TP11_1237 Mcapitata_ATHC_TP11_1237 Mcapitata      ATHC      TP11
## Mcapitata_ATHC_TP11_2188 Mcapitata_ATHC_TP11_2188 Mcapitata      ATHC      TP11
## Mcapitata_ATHC_TP1_1218   Mcapitata_ATHC_TP1_1218 Mcapitata      ATHC       TP1
## Mcapitata_ATHC_TP11_2756 Mcapitata_ATHC_TP11_2756 Mcapitata      ATHC      TP11
## Mcapitata_ATHC_TP1_1826   Mcapitata_ATHC_TP1_1826 Mcapitata      ATHC       TP1
## Mcapitata_ATHC_TP1_2068   Mcapitata_ATHC_TP1_2068 Mcapitata      ATHC       TP1
## Mcapitata_ATHC_TP12_1154 Mcapitata_ATHC_TP12_1154 Mcapitata      ATHC      TP12
## Mcapitata_ATHC_TP12_2736 Mcapitata_ATHC_TP12_2736 Mcapitata      ATHC      TP12
## Mcapitata_ATHC_TP12_2990 Mcapitata_ATHC_TP12_2990 Mcapitata      ATHC      TP12
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
## Mcapitata_HTAC_TP10_1315 Mcapitata_HTAC_TP10_1315 Mcapitata      HTAC      TP10
## Mcapitata_HTAC_TP10_1478 Mcapitata_HTAC_TP10_1478 Mcapitata      HTAC      TP10
## Mcapitata_HTAC_TP10_1754 Mcapitata_HTAC_TP10_1754 Mcapitata      HTAC      TP10
## Mcapitata_HTAC_TP11_1248 Mcapitata_HTAC_TP11_1248 Mcapitata      HTAC      TP11
## Mcapitata_HTAC_TP11_1562 Mcapitata_HTAC_TP11_1562 Mcapitata      HTAC      TP11
## Mcapitata_HTAC_TP11_2380 Mcapitata_HTAC_TP11_2380 Mcapitata      HTAC      TP11
## Mcapitata_HTAC_TP1_1579   Mcapitata_HTAC_TP1_1579 Mcapitata      HTAC       TP1
## Mcapitata_HTAC_TP1_2153   Mcapitata_HTAC_TP1_2153 Mcapitata      HTAC       TP1
## Mcapitata_HTAC_TP12_1632 Mcapitata_HTAC_TP12_1632 Mcapitata      HTAC      TP12
## Mcapitata_HTAC_TP12_1729 Mcapitata_HTAC_TP12_1729 Mcapitata      HTAC      TP12
## Mcapitata_HTAC_TP1_2183   Mcapitata_HTAC_TP1_2183 Mcapitata      HTAC       TP1
## Mcapitata_HTAC_TP12_2007 Mcapitata_HTAC_TP12_2007 Mcapitata      HTAC      TP12
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
## Mcapitata_HTHC_TP10_1074 Mcapitata_HTHC_TP10_1074 Mcapitata      HTHC      TP10
## Mcapitata_HTHC_TP10_1332 Mcapitata_HTHC_TP10_1332 Mcapitata      HTHC      TP10
## Mcapitata_HTHC_TP10_1689 Mcapitata_HTHC_TP10_1689 Mcapitata      HTHC      TP10
## Mcapitata_HTHC_TP11_1178 Mcapitata_HTHC_TP11_1178 Mcapitata      HTHC      TP11
## Mcapitata_HTHC_TP11_1270 Mcapitata_HTHC_TP11_1270 Mcapitata      HTHC      TP11
## Mcapitata_HTHC_TP1_1145   Mcapitata_HTHC_TP1_1145 Mcapitata      HTHC       TP1
## Mcapitata_HTHC_TP11_2511 Mcapitata_HTHC_TP11_2511 Mcapitata      HTHC      TP11
## Mcapitata_HTHC_TP1_1323   Mcapitata_HTHC_TP1_1323 Mcapitata      HTHC       TP1
## Mcapitata_HTHC_TP1_2081   Mcapitata_HTHC_TP1_2081 Mcapitata      HTHC       TP1
## Mcapitata_HTHC_TP12_1140 Mcapitata_HTHC_TP12_1140 Mcapitata      HTHC      TP12
## Mcapitata_HTHC_TP12_1274 Mcapitata_HTHC_TP12_1274 Mcapitata      HTHC      TP12
## Mcapitata_HTHC_TP12_2190 Mcapitata_HTHC_TP12_2190 Mcapitata      HTHC      TP12
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
##                          plugid            reef reef_color ploidy ploidy_color
## Mcapitata_ATAC_TP10_1095   1095      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP10_1561   1561      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATAC_TP10_1631   1631            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP1_1037    1037      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP11_1076   1076      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP11_1644   1644            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP11_2302   2302         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATAC_TP1_1600    1600      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP1_1652    1652            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP12_1120   1120      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP12_1452   1452      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP12_2403   2403         Reef.18    #d8daeb      2      #01665e
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
## Mcapitata_ATAC_TP6_2402    2402         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATAC_TP7_1058    1058 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_ATAC_TP7_1455    1455      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP7_1499    1499      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP8_1083    1083      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP8_1436    1436      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP8_1779    1779      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATAC_TP9_1121    1121            HIMB    #b35806      2      #01665e
## Mcapitata_ATAC_TP9_1420    1420      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATAC_TP9_1580    1580      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP10_1204   1204      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP10_2554   2554         Reef.18    #d8daeb      2      #01665e
## Mcapitata_ATHC_TP10_2737   2737      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP11_1237   1237      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP11_2188   2188      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP1_1218    1218      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP11_2756   2756      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP1_1826    1826      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_ATHC_TP1_2068    2068      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_ATHC_TP12_1154   1154            HIMB    #b35806      2      #01665e
## Mcapitata_ATHC_TP12_2736   2736      Reef.42.43    #542788      2      #01665e
## Mcapitata_ATHC_TP12_2990   2990         Reef.18    #d8daeb      2      #01665e
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
## Mcapitata_HTAC_TP10_1315   1315 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTAC_TP10_1478   1478      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTAC_TP10_1754   1754            HIMB    #b35806      2      #01665e
## Mcapitata_HTAC_TP11_1248   1248            HIMB    #b35806      2      #01665e
## Mcapitata_HTAC_TP11_1562   1562      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTAC_TP11_2380   2380         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP1_1579    1579 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTAC_TP1_2153    2153         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP12_1632   1632            HIMB    #b35806      4      #7570b3
## Mcapitata_HTAC_TP12_1729   1729         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTAC_TP1_2183    2183      Reef.42.43    #542788      2      #01665e
## Mcapitata_HTAC_TP12_2007   2007      Reef.42.43    #542788      2      #01665e
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
## Mcapitata_HTHC_TP10_1074   1074      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP10_1332   1332      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP10_1689   1689      Reef.11.13    #fee0b6      2      #01665e
## Mcapitata_HTHC_TP11_1178   1178 Lilipuna.Fringe    #f1a340      2      #01665e
## Mcapitata_HTHC_TP11_1270   1270            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP1_1145    1145            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP11_2511   2511         Reef.18    #d8daeb      2      #01665e
## Mcapitata_HTHC_TP1_1323    1323      Reef.35.36    #998ec3      2      #01665e
## Mcapitata_HTHC_TP1_2081    2081            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP12_1140   1140            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP12_1274   1274            HIMB    #b35806      2      #01665e
## Mcapitata_HTHC_TP12_2190   2190      Reef.42.43    #542788      2      #01665e
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
##                              group group_color
## Mcapitata_ATAC_TP10_1095 Ungrouped     #808080
## Mcapitata_ATAC_TP10_1561 Ungrouped     #808080
## Mcapitata_ATAC_TP10_1631 Ungrouped     #808080
## Mcapitata_ATAC_TP1_1037  Ungrouped     #808080
## Mcapitata_ATAC_TP11_1076 Ungrouped     #808080
## Mcapitata_ATAC_TP11_1644 Ungrouped     #808080
## Mcapitata_ATAC_TP11_2302    Group1     #1f78b4
## Mcapitata_ATAC_TP1_1600  Ungrouped     #808080
## Mcapitata_ATAC_TP1_1652  Ungrouped     #808080
## Mcapitata_ATAC_TP12_1120 Ungrouped     #808080
## Mcapitata_ATAC_TP12_1452 Ungrouped     #808080
## Mcapitata_ATAC_TP12_2403    Group2     #33a02c
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
## Mcapitata_ATAC_TP6_2402     Group2     #33a02c
## Mcapitata_ATAC_TP7_1058  Ungrouped     #808080
## Mcapitata_ATAC_TP7_1455  Ungrouped     #808080
## Mcapitata_ATAC_TP7_1499  Ungrouped     #808080
## Mcapitata_ATAC_TP8_1083  Ungrouped     #808080
## Mcapitata_ATAC_TP8_1436  Ungrouped     #808080
## Mcapitata_ATAC_TP8_1779  Ungrouped     #808080
## Mcapitata_ATAC_TP9_1121  Ungrouped     #808080
## Mcapitata_ATAC_TP9_1420  Ungrouped     #808080
## Mcapitata_ATAC_TP9_1580  Ungrouped     #808080
## Mcapitata_ATHC_TP10_1204 Ungrouped     #808080
## Mcapitata_ATHC_TP10_2554 Ungrouped     #808080
## Mcapitata_ATHC_TP10_2737 Ungrouped     #808080
## Mcapitata_ATHC_TP11_1237 Ungrouped     #808080
## Mcapitata_ATHC_TP11_2188 Ungrouped     #808080
## Mcapitata_ATHC_TP1_1218  Ungrouped     #808080
## Mcapitata_ATHC_TP11_2756 Ungrouped     #808080
## Mcapitata_ATHC_TP1_1826  Ungrouped     #808080
## Mcapitata_ATHC_TP1_2068  Ungrouped     #808080
## Mcapitata_ATHC_TP12_1154 Ungrouped     #808080
## Mcapitata_ATHC_TP12_2736 Ungrouped     #808080
## Mcapitata_ATHC_TP12_2990 Ungrouped     #808080
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
## Mcapitata_HTAC_TP10_1315 Ungrouped     #808080
## Mcapitata_HTAC_TP10_1478 Ungrouped     #808080
## Mcapitata_HTAC_TP10_1754 Ungrouped     #808080
## Mcapitata_HTAC_TP11_1248 Ungrouped     #808080
## Mcapitata_HTAC_TP11_1562 Ungrouped     #808080
## Mcapitata_HTAC_TP11_2380    Group1     #1f78b4
## Mcapitata_HTAC_TP1_1579  Ungrouped     #808080
## Mcapitata_HTAC_TP1_2153  Ungrouped     #808080
## Mcapitata_HTAC_TP12_1632 Ungrouped     #808080
## Mcapitata_HTAC_TP12_1729 Ungrouped     #808080
## Mcapitata_HTAC_TP1_2183  Ungrouped     #808080
## Mcapitata_HTAC_TP12_2007 Ungrouped     #808080
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
## Mcapitata_HTHC_TP10_1074 Ungrouped     #808080
## Mcapitata_HTHC_TP10_1332 Ungrouped     #808080
## Mcapitata_HTHC_TP10_1689 Ungrouped     #808080
## Mcapitata_HTHC_TP11_1178 Ungrouped     #808080
## Mcapitata_HTHC_TP11_1270 Ungrouped     #808080
## Mcapitata_HTHC_TP1_1145  Ungrouped     #808080
## Mcapitata_HTHC_TP11_2511 Ungrouped     #808080
## Mcapitata_HTHC_TP1_1323  Ungrouped     #808080
## Mcapitata_HTHC_TP1_2081  Ungrouped     #808080
## Mcapitata_HTHC_TP12_1140 Ungrouped     #808080
## Mcapitata_HTHC_TP12_1274 Ungrouped     #808080
## Mcapitata_HTHC_TP12_2190 Ungrouped     #808080
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
## Lilipuna.Fringe               2               4       Ungrouped          Group1 
##       "#f1a340"       "#01665e"       "#7570b3"       "#808080"       "#1f78b4" 
##          Group2 
##       "#33a02c"
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
<div id="htmlwidget-6cc3a7674f607a1eca7c" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-6cc3a7674f607a1eca7c">{"x":{"visdat":{"6aac70a0c7d8":["function () ","plotlyVisDat"]},"cur_data":"6aac70a0c7d8","attrs":{"6aac70a0c7d8":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#01665e","#7570b3","#808080","#1f78b4","#33a02c"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Group","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[0.0889479242601857,0.09092310626426],"y":[0.624161687513457,0.6400298931495],"text":["Mcapitata_ATAC_TP11_2302","Mcapitata_HTAC_TP11_2380"],"mode":"markers","marker":{"color":"rgba(31,120,180,1)","size":11,"line":{"color":"rgba(31,120,180,1)"}},"type":"scatter","name":"Group1","textfont":{"color":"rgba(31,120,180,1)"},"error_y":{"color":"rgba(31,120,180,1)"},"error_x":{"color":"rgba(31,120,180,1)"},"line":{"color":"rgba(31,120,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0861859616461228,0.0868350795172354],"y":[0.0711207334406863,0.0681824583681851],"text":["Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP6_2402"],"mode":"markers","marker":{"color":"rgba(51,160,44,1)","size":11,"line":{"color":"rgba(51,160,44,1)"}},"type":"scatter","name":"Group2","textfont":{"color":"rgba(51,160,44,1)"},"error_y":{"color":"rgba(51,160,44,1)"},"error_x":{"color":"rgba(51,160,44,1)"},"line":{"color":"rgba(51,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0527365612881002,-0.0543805285080944,-0.00907460413367875,-0.108355785998507,-0.184213341100805,-0.09247053150177,-0.0247710199682509,0.0737660204352675,0.0305108926840005,-0.181398012019601,0.0640025596993812,0.0813334102681023,0.0771224347278394,0.0167366928372258,0.0776292147189089,0.0840127169455604,-0.0471570905146239,0.0478644795430772,-0.00660194602055129,-0.18566221088566,-0.00463215467427667,-0.11475992676886,0.0451460193046511,-0.0377292805428004,0.0370209170823294,-0.076054291666041,-0.0840093440919089,0.0819253182981672,0.00841569917074762,-0.121275529558162,0.0757892873201796,-0.039364275195964,0.0443645816586402,0.0521742878051731,0.0335350141464862,0.0772693395311863,0.00806901299075618,-0.0544886318526593,-0.120349004908409,0.0238526192209174,0.0499514472470783,0.0627267185076457,-0.0934958642123163,-0.109174213067995,-0.081741172459687,-0.172330720860554,0.00577993899172363,0.0727944933326374,-0.178785641167139,0.0473742984601693,0.0474169582811036,0.0826103699803756,0.0742501702741995,-0.0193070210245759,-0.0281893042508799,-0.0171424521578102,-0.0219232146403328,-0.182385135307786,-0.190141964042346,0.0228396535887499,-0.194554416711963,-0.0553505359001984,0.0271483413810599,0.0377669896009464,-0.182101165617165,-0.0514686046024263,0.0262533863504204,0.0593661769285165,0.048051770850318,-0.0858196694078485,-0.0666308127115562,0.0380295080167764,-0.0257800588075927,0.0593212572156575,0.056580795038582,-0.111317891839658,0.0597130986768337,0.0810841136170262,0.078315179811192,0.0843554989402621,0.0767524161030123,0.0629908583367122,0.0304870341283838,-0.184642276080596,0.0348462076557378,-0.123475310281597,-0.0144399711474971,0.0584263533905763,-0.177874095946366,-0.0335608317599685,-0.0977293163059628,0.0750646398407285,-0.0428332888893134,0.0180466549711767,-0.0522181139283315,-0.00188240556089456,0.0843715829375297,-0.175812475881633,0.0582941363354656,0.0610190772396278,0.0810302876141942,0.0678058360054623,-0.0464311699208948,-0.0154821162724179,0.0844028176591667,0.076217743987794,0.0746329019393565,0.0192466635492942,-0.148185437720757,-0.00133257992576647,0.0835529279053348,0.0797693555356296,-0.181691570831866,-0.0339855820773435,0.0717528035142179,0.0624017028157107,-0.0253404367117615,-0.0435859602429304,0.057950054067745,0.00273161791659655,-0.0280008526333427,0.0703487517815327,0.0642400891175913,-0.190128968724382,-0.121514855459929,0.0326902830355003,-0.18238390661973,0.064160595280505],"y":[-0.00129394712661751,-0.099253864441658,-0.0677862338325736,0.0106092659710857,0.0350613702890352,0.0108406176006241,-0.0393489960045606,-0.0434730155011276,-0.0486052863217551,0.0423161333600887,-0.0246606157733388,-0.0283568223973203,-0.0170237192652383,0.0379070457562692,-0.0271528131228215,-0.00776735123801566,-0.0653785120100204,-0.0127646151727145,0.000522193719633659,0.0308035819997685,-0.0932694070682311,0.00342298387817902,-0.0269286376173691,-0.0198539620722995,0.0462846279831741,-0.0232353516740996,-0.00468233377167214,-0.0425430729640562,-0.0267189617934759,-0.0121931580060471,0.0161388449708949,-0.0834109736615799,-0.0716461497676098,0.0112869165502641,0.0310764060440811,0.00381159017030018,-0.0048065342778704,0.0293346129855834,0.0317137435534603,-0.0133279032335982,0.0106814793667712,-0.00649484068778893,0.017646801784634,0.0112129969305108,-0.00463148786766793,0.0548375916802918,-0.0121702474032001,-0.00828199100678087,0.0320035211768657,-0.0693068896888515,-0.028607606689892,-0.0285463027251396,-0.0560344535419861,0.0354989561225742,0.0107934925980308,-0.059188050509131,-0.0426899738250583,0.00929038404341888,0.0335116622411095,-0.0370281245465273,0.0166554266637875,-0.0423966191744478,-0.0390744296967309,-0.0153207147773141,-0.00136168173519049,-0.080703770118225,0.0317710040187157,0.0344666796364336,-0.00924352547484828,0.0249471994260217,0.131824708899737,0.0333855833733984,-0.0715328433810549,-0.0486518407837484,0.026244979109336,-0.00116436412085055,-0.0305087614353352,-0.00342660522872357,-0.0174759885209238,-0.0316248404458193,-0.0432771189126409,-0.00189813523674944,-0.0172863569305958,0.0276114594074663,0.00638963490906601,0.00240450879198285,-0.0849103620090995,-0.0318945348790096,0.0317663436279504,-0.051082378072563,-0.026642684505122,-0.0239439291119711,-0.0341900627039689,0.000662442161083721,-0.00364263003490803,-0.0770997057952562,-0.0125907705377221,0.0205787854420928,-0.0253866874132107,-0.0226120185914045,-0.0430391357122137,0.00180964460188482,-0.00460762866530396,-0.0386094593532937,-0.00243298272101401,0.00335668264190462,0.00191510677534522,0.012692619507298,0.0390714705666226,-0.0607652756439441,-0.0493224878961267,-0.00430516400169643,0.0178646934184598,-0.0492293284686281,-0.000224638128483849,-0.0547390734333192,-0.0333923686761083,-0.0335425975062624,-0.0654890327418421,-0.0580619305499708,-0.00417108001498296,-0.0367082735260535,-0.0389804605821278,0.020228903478627,-0.0192470929749976,-0.0161671327980797,0.0125969209106861,-0.0755536703887436],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"color":"rgba(128,128,128,1)","size":11,"line":{"color":"rgba(128,128,128,1)"}},"type":"scatter","name":"Ungrouped","textfont":{"color":"rgba(128,128,128,1)"},"error_y":{"color":"rgba(128,128,128,1)"},"error_x":{"color":"rgba(128,128,128,1)"},"line":{"color":"rgba(128,128,128,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-981efe28cc89d461ca21" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-981efe28cc89d461ca21">{"x":{"visdat":{"6aac4ba06fac":["function () ","plotlyVisDat"]},"cur_data":"6aac4ba06fac","attrs":{"6aac4ba06fac":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#01665e","#7570b3","#808080","#1f78b4","#33a02c"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Reef","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.00907460413367875,-0.09247053150177,0.0737660204352675,0.0640025596993812,0.0771224347278394,0.0776292147189089,0.0840127169455604,0.0819253182981672,0.0238526192209174,-0.081741172459687,-0.0171424521578102,-0.0219232146403328,-0.0553505359001984,-0.0514686046024263,0.0262533863504204,-0.0666308127115562,0.0843554989402621,0.0304870341283838,0.0610190772396278,0.0810302876141942,-0.0154821162724179,0.0844028176591667,0.076217743987794,0.0192466635492942,-0.148185437720757,-0.121514855459929],"y":[-0.0677862338325736,0.0108406176006241,-0.0434730155011276,-0.0246606157733388,-0.0170237192652383,-0.0271528131228215,-0.00776735123801566,-0.0425430729640562,-0.0133279032335982,-0.00463148786766793,-0.059188050509131,-0.0426899738250583,-0.0423966191744478,-0.080703770118225,0.0317710040187157,0.131824708899737,-0.0316248404458193,-0.0172863569305958,-0.0226120185914045,-0.0430391357122137,-0.0386094593532937,-0.00243298272101401,0.00335668264190462,0.012692619507298,0.0390714705666226,-0.0192470929749976],"text":["Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP9_1121","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP9_2862","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1997","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP8_2067"],"mode":"markers","marker":{"color":"rgba(179,88,6,1)","size":11,"line":{"color":"rgba(179,88,6,1)"}},"type":"scatter","name":"HIMB","textfont":{"color":"rgba(179,88,6,1)"},"error_y":{"color":"rgba(179,88,6,1)"},"error_x":{"color":"rgba(179,88,6,1)"},"line":{"color":"rgba(179,88,6,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.11475992676886,0.0377669896009464,0.048051770850318,0.0629908583367122,-0.0428332888893134,0.0582941363354656,0.0717528035142179,0.064160595280505],"y":[0.00342298387817902,-0.0153207147773141,-0.00924352547484828,-0.00189813523674944,-0.0341900627039689,-0.0253866874132107,-0.000224638128483849,-0.0755536703887436],"text":["Mcapitata_ATAC_TP7_1058","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP9_1306","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"color":"rgba(241,163,64,1)","size":11,"line":{"color":"rgba(241,163,64,1)"}},"type":"scatter","name":"Lilipuna.Fringe","textfont":{"color":"rgba(241,163,64,1)"},"error_y":{"color":"rgba(241,163,64,1)"},"error_x":{"color":"rgba(241,163,64,1)"},"line":{"color":"rgba(241,163,64,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0543805285080944,0.0813334102681023,0.0757892873201796,0.0772693395311863,-0.0544886318526593,0.0473742984601693,0.0826103699803756,0.0271483413810599,-0.182101165617165,0.0593661769285165,0.0767524161030123,-0.0144399711474971,0.0180466549711767,-0.00188240556089456,0.0843715829375297,-0.175812475881633,0.0797693555356296,0.0624017028157107,0.057950054067745,0.0642400891175913,0.0326902830355003,-0.18238390661973],"y":[-0.099253864441658,-0.0283568223973203,0.0161388449708949,0.00381159017030018,0.0293346129855834,-0.0693068896888515,-0.0285463027251396,-0.0390744296967309,-0.00136168173519049,0.0344666796364336,-0.0432771189126409,-0.0849103620090995,0.000662442161083721,-0.0770997057952562,-0.0125907705377221,0.0205787854420928,-0.00430516400169643,-0.0547390734333192,-0.0654890327418421,-0.0389804605821278,-0.0161671327980797,0.0125969209106861],"text":["Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP3_1548","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP9_1467","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331"],"mode":"markers","marker":{"color":"rgba(254,224,182,1)","size":11,"line":{"color":"rgba(254,224,182,1)"}},"type":"scatter","name":"Reef.11.13","textfont":{"color":"rgba(254,224,182,1)"},"error_y":{"color":"rgba(254,224,182,1)"},"error_x":{"color":"rgba(254,224,182,1)"},"line":{"color":"rgba(254,224,182,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0889479242601857,0.0861859616461228,-0.0471570905146239,0.0868350795172354,-0.039364275195964,0.0627267185076457,-0.0934958642123163,0.00577993899172363,-0.0193070210245759,0.09092310626426,-0.0858196694078485,0.0380295080167764,0.0348462076557378,-0.0977293163059628,0.0750646398407285,-0.0522181139283315,0.0678058360054623,-0.00133257992576647,-0.0339855820773435,0.0703487517815327],"y":[0.624161687513457,0.0711207334406863,-0.0653785120100204,0.0681824583681851,-0.0834109736615799,-0.00649484068778893,0.017646801784634,-0.0121702474032001,0.0354989561225742,0.6400298931495,0.0249471994260217,0.0333855833733984,0.00638963490906601,-0.026642684505122,-0.0239439291119711,-0.00364263003490803,0.00180964460188482,-0.0607652756439441,-0.0492293284686281,-0.0367082735260535],"text":["Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP6_2402","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP6_2555","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP7_2419"],"mode":"markers","marker":{"color":"rgba(216,218,235,1)","size":11,"line":{"color":"rgba(216,218,235,1)"}},"type":"scatter","name":"Reef.18","textfont":{"color":"rgba(216,218,235,1)"},"error_y":{"color":"rgba(216,218,235,1)"},"error_x":{"color":"rgba(216,218,235,1)"},"line":{"color":"rgba(216,218,235,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0527365612881002,-0.108355785998507,-0.184213341100805,-0.181398012019601,0.0167366928372258,0.0451460193046511,-0.0377292805428004,0.0370209170823294,-0.076054291666041,0.00841569917074762,-0.121275529558162,0.0443645816586402,0.0521742878051731,0.0335350141464862,-0.120349004908409,-0.172330720860554,0.0474169582811036,0.0742501702741995,-0.182385135307786,0.056580795038582,0.0810841136170262,0.078315179811192,-0.184642276080596,-0.123475310281597,-0.177874095946366,-0.0464311699208948,-0.0253404367117615,-0.0435859602429304],"y":[-0.00129394712661751,0.0106092659710857,0.0350613702890352,0.0423161333600887,0.0379070457562692,-0.0269286376173691,-0.0198539620722995,0.0462846279831741,-0.0232353516740996,-0.0267189617934759,-0.0121931580060471,-0.0716461497676098,0.0112869165502641,0.0310764060440811,0.0317137435534603,0.0548375916802918,-0.028607606689892,-0.0560344535419861,0.00929038404341888,0.026244979109336,-0.00342660522872357,-0.0174759885209238,0.0276114594074663,0.00240450879198285,0.0317663436279504,-0.00460762866530396,-0.0333923686761083,-0.0335425975062624],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP8_1260","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1722","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317"],"mode":"markers","marker":{"color":"rgba(153,142,195,1)","size":11,"line":{"color":"rgba(153,142,195,1)"}},"type":"scatter","name":"Reef.35.36","textfont":{"color":"rgba(153,142,195,1)"},"error_y":{"color":"rgba(153,142,195,1)"},"error_x":{"color":"rgba(153,142,195,1)"},"line":{"color":"rgba(153,142,195,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0247710199682509,0.0305108926840005,0.0478644795430772,-0.00660194602055129,-0.18566221088566,-0.00463215467427667,-0.0840093440919089,0.00806901299075618,0.0499514472470783,-0.109174213067995,0.0727944933326374,-0.178785641167139,-0.0281893042508799,-0.190141964042346,0.0228396535887499,-0.194554416711963,-0.0257800588075927,0.0593212572156575,-0.111317891839658,0.0597130986768337,0.0584263533905763,-0.0335608317599685,0.0746329019393565,0.0835529279053348,-0.181691570831866,0.00273161791659655,-0.0280008526333427,-0.190128968724382],"y":[-0.0393489960045606,-0.0486052863217551,-0.0127646151727145,0.000522193719633659,0.0308035819997685,-0.0932694070682311,-0.00468233377167214,-0.0048065342778704,0.0106814793667712,0.0112129969305108,-0.00828199100678087,0.0320035211768657,0.0107934925980308,0.0335116622411095,-0.0370281245465273,0.0166554266637875,-0.0715328433810549,-0.0486518407837484,-0.00116436412085055,-0.0305087614353352,-0.0318945348790096,-0.051082378072563,0.00191510677534522,-0.0493224878961267,0.0178646934184598,-0.0580619305499708,-0.00417108001498296,0.020228903478627],"text":["Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP8_1779","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP8_1235","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP8_1246"],"mode":"markers","marker":{"color":"rgba(84,39,136,1)","size":11,"line":{"color":"rgba(84,39,136,1)"}},"type":"scatter","name":"Reef.42.43","textfont":{"color":"rgba(84,39,136,1)"},"error_y":{"color":"rgba(84,39,136,1)"},"error_x":{"color":"rgba(84,39,136,1)"},"line":{"color":"rgba(84,39,136,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-54dab9e22688bc5d929e" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-54dab9e22688bc5d929e">{"x":{"visdat":{"6aac1203abfe":["function () ","plotlyVisDat"]},"cur_data":"6aac1203abfe","attrs":{"6aac1203abfe":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#01665e","#7570b3","#808080","#1f78b4","#33a02c"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Ploidy","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0527365612881002,-0.0543805285080944,-0.00907460413367875,-0.108355785998507,-0.184213341100805,-0.09247053150177,0.0889479242601857,-0.0247710199682509,0.0737660204352675,0.0305108926840005,-0.181398012019601,0.0861859616461228,0.0640025596993812,0.0813334102681023,0.0771224347278394,0.0167366928372258,0.0776292147189089,0.0840127169455604,-0.0471570905146239,0.0478644795430772,-0.00660194602055129,-0.18566221088566,-0.00463215467427667,0.0868350795172354,-0.11475992676886,0.0451460193046511,-0.0377292805428004,0.0370209170823294,-0.076054291666041,-0.0840093440919089,0.0819253182981672,0.00841569917074762,-0.121275529558162,0.0757892873201796,-0.039364275195964,0.0443645816586402,0.0521742878051731,0.0335350141464862,0.0772693395311863,0.00806901299075618,-0.0544886318526593,-0.120349004908409,0.0238526192209174,0.0499514472470783,0.0627267185076457,-0.0934958642123163,-0.109174213067995,-0.081741172459687,-0.172330720860554,0.00577993899172363,0.0727944933326374,-0.178785641167139,0.0473742984601693,0.0474169582811036,0.0826103699803756,0.0742501702741995,-0.0193070210245759,-0.0281893042508799,-0.0171424521578102,-0.0219232146403328,-0.182385135307786,-0.190141964042346,0.0228396535887499,-0.194554416711963,-0.0553505359001984,0.0271483413810599,0.0377669896009464,-0.182101165617165,-0.0514686046024263,0.0262533863504204,0.0593661769285165,0.09092310626426,0.048051770850318,-0.0858196694078485,-0.0666308127115562,0.0380295080167764,-0.0257800588075927,0.0593212572156575,0.056580795038582,-0.111317891839658,0.0597130986768337,0.0810841136170262,0.078315179811192,0.0843554989402621,0.0767524161030123,0.0629908583367122,0.0304870341283838,-0.184642276080596,0.0348462076557378,-0.123475310281597,-0.0144399711474971,0.0584263533905763,-0.177874095946366,-0.0335608317599685,-0.0977293163059628,0.0750646398407285,-0.0428332888893134,0.0180466549711767,-0.0522181139283315,-0.00188240556089456,0.0843715829375297,-0.175812475881633,0.0582941363354656,0.0610190772396278,0.0810302876141942,0.0678058360054623,-0.0464311699208948,-0.0154821162724179,0.0844028176591667,0.076217743987794,0.0746329019393565,0.0192466635492942,-0.148185437720757,-0.00133257992576647,0.0835529279053348,0.0797693555356296,-0.181691570831866,-0.0339855820773435,0.0717528035142179,0.0624017028157107,-0.0253404367117615,-0.0435859602429304,0.057950054067745,0.00273161791659655,-0.0280008526333427,0.0703487517815327,0.0642400891175913,-0.190128968724382,-0.121514855459929,0.0326902830355003,-0.18238390661973,0.064160595280505],"y":[-0.00129394712661751,-0.099253864441658,-0.0677862338325736,0.0106092659710857,0.0350613702890352,0.0108406176006241,0.624161687513457,-0.0393489960045606,-0.0434730155011276,-0.0486052863217551,0.0423161333600887,0.0711207334406863,-0.0246606157733388,-0.0283568223973203,-0.0170237192652383,0.0379070457562692,-0.0271528131228215,-0.00776735123801566,-0.0653785120100204,-0.0127646151727145,0.000522193719633659,0.0308035819997685,-0.0932694070682311,0.0681824583681851,0.00342298387817902,-0.0269286376173691,-0.0198539620722995,0.0462846279831741,-0.0232353516740996,-0.00468233377167214,-0.0425430729640562,-0.0267189617934759,-0.0121931580060471,0.0161388449708949,-0.0834109736615799,-0.0716461497676098,0.0112869165502641,0.0310764060440811,0.00381159017030018,-0.0048065342778704,0.0293346129855834,0.0317137435534603,-0.0133279032335982,0.0106814793667712,-0.00649484068778893,0.017646801784634,0.0112129969305108,-0.00463148786766793,0.0548375916802918,-0.0121702474032001,-0.00828199100678087,0.0320035211768657,-0.0693068896888515,-0.028607606689892,-0.0285463027251396,-0.0560344535419861,0.0354989561225742,0.0107934925980308,-0.059188050509131,-0.0426899738250583,0.00929038404341888,0.0335116622411095,-0.0370281245465273,0.0166554266637875,-0.0423966191744478,-0.0390744296967309,-0.0153207147773141,-0.00136168173519049,-0.080703770118225,0.0317710040187157,0.0344666796364336,0.6400298931495,-0.00924352547484828,0.0249471994260217,0.131824708899737,0.0333855833733984,-0.0715328433810549,-0.0486518407837484,0.026244979109336,-0.00116436412085055,-0.0305087614353352,-0.00342660522872357,-0.0174759885209238,-0.0316248404458193,-0.0432771189126409,-0.00189813523674944,-0.0172863569305958,0.0276114594074663,0.00638963490906601,0.00240450879198285,-0.0849103620090995,-0.0318945348790096,0.0317663436279504,-0.051082378072563,-0.026642684505122,-0.0239439291119711,-0.0341900627039689,0.000662442161083721,-0.00364263003490803,-0.0770997057952562,-0.0125907705377221,0.0205787854420928,-0.0253866874132107,-0.0226120185914045,-0.0430391357122137,0.00180964460188482,-0.00460762866530396,-0.0386094593532937,-0.00243298272101401,0.00335668264190462,0.00191510677534522,0.012692619507298,0.0390714705666226,-0.0607652756439441,-0.0493224878961267,-0.00430516400169643,0.0178646934184598,-0.0492293284686281,-0.000224638128483849,-0.0547390734333192,-0.0333923686761083,-0.0335425975062624,-0.0654890327418421,-0.0580619305499708,-0.00417108001498296,-0.0367082735260535,-0.0389804605821278,0.020228903478627,-0.0192470929749976,-0.0161671327980797,0.0125969209106861,-0.0755536703887436],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"colorbar":{"title":"ploidy","ticklen":2},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(197,175,191,1)"],["0.0833333333333333","rgba(238,210,185,1)"],["0.125","rgba(239,189,139,1)"],["0.166666666666667","rgba(208,132,70,1)"],["0.208333333333333","rgba(185,99,33,1)"],["0.25","rgba(209,152,123,1)"],["0.291666666666667","rgba(217,207,216,1)"],["0.333333333333333","rgba(173,156,202,1)"],["0.375","rgba(119,83,161,1)"],["0.416666666666667","rgba(118,58,128,1)"],["0.458333333333333","rgba(183,108,103,1)"],["0.5","rgba(241,163,64,1)"],["0.541666666666667","rgba(161,139,80,1)"],["0.583333333333333","rgba(73,113,91,1)"],["0.625","rgba(48,106,115,1)"],["0.666666666666667","rgba(89,110,150,1)"],["0.708333333333333","rgba(119,113,175,1)"],["0.75","rgba(124,120,154,1)"],["0.791666666666667","rgba(128,127,132,1)"],["0.833333333333333","rgba(109,125,145,1)"],["0.875","rgba(73,122,167,1)"],["0.916666666666667","rgba(54,127,159,1)"],["0.958333333333333","rgba(68,143,107,1)"],["1","rgba(51,160,44,1)"]],"showscale":false,"color":[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"size":11,"line":{"colorbar":{"title":"","ticklen":2},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(197,175,191,1)"],["0.0833333333333333","rgba(238,210,185,1)"],["0.125","rgba(239,189,139,1)"],["0.166666666666667","rgba(208,132,70,1)"],["0.208333333333333","rgba(185,99,33,1)"],["0.25","rgba(209,152,123,1)"],["0.291666666666667","rgba(217,207,216,1)"],["0.333333333333333","rgba(173,156,202,1)"],["0.375","rgba(119,83,161,1)"],["0.416666666666667","rgba(118,58,128,1)"],["0.458333333333333","rgba(183,108,103,1)"],["0.5","rgba(241,163,64,1)"],["0.541666666666667","rgba(161,139,80,1)"],["0.583333333333333","rgba(73,113,91,1)"],["0.625","rgba(48,106,115,1)"],["0.666666666666667","rgba(89,110,150,1)"],["0.708333333333333","rgba(119,113,175,1)"],["0.75","rgba(124,120,154,1)"],["0.791666666666667","rgba(128,127,132,1)"],["0.833333333333333","rgba(109,125,145,1)"],["0.875","rgba(73,122,167,1)"],["0.916666666666667","rgba(54,127,159,1)"],["0.958333333333333","rgba(68,143,107,1)"],["1","rgba(51,160,44,1)"]],"showscale":false,"color":[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}},"type":"scatter","xaxis":"x","yaxis":"y","frame":null},{"x":[-0.194554416711963,0.09092310626426],"y":[-0.099253864441658,0.6400298931495],"type":"scatter","mode":"markers","opacity":0,"hoverinfo":"none","showlegend":false,"marker":{"colorbar":{"title":"ploidy","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(197,175,191,1)"],["0.0833333333333333","rgba(238,210,185,1)"],["0.125","rgba(239,189,139,1)"],["0.166666666666667","rgba(208,132,70,1)"],["0.208333333333333","rgba(185,99,33,1)"],["0.25","rgba(209,152,123,1)"],["0.291666666666667","rgba(217,207,216,1)"],["0.333333333333333","rgba(173,156,202,1)"],["0.375","rgba(119,83,161,1)"],["0.416666666666667","rgba(118,58,128,1)"],["0.458333333333333","rgba(183,108,103,1)"],["0.5","rgba(241,163,64,1)"],["0.541666666666667","rgba(161,139,80,1)"],["0.583333333333333","rgba(73,113,91,1)"],["0.625","rgba(48,106,115,1)"],["0.666666666666667","rgba(89,110,150,1)"],["0.708333333333333","rgba(119,113,175,1)"],["0.75","rgba(124,120,154,1)"],["0.791666666666667","rgba(128,127,132,1)"],["0.833333333333333","rgba(109,125,145,1)"],["0.875","rgba(73,122,167,1)"],["0.916666666666667","rgba(54,127,159,1)"],["0.958333333333333","rgba(68,143,107,1)"],["1","rgba(51,160,44,1)"]],"showscale":true,"color":[2,4],"line":{"color":"rgba(255,127,14,1)"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-8e6071ba625830f821c9" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-8e6071ba625830f821c9">{"x":{"visdat":{"6aac49423a66":["function () ","plotlyVisDat"]},"cur_data":"6aac49423a66","attrs":{"6aac49423a66":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#01665e","#7570b3","#808080","#1f78b4","#33a02c"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Group","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.29916825882356,-0.311412109391442],"y":[-0.603312761641301,-0.618994536873599],"text":["Mcapitata_ATAC_TP11_2302","Mcapitata_HTAC_TP11_2380"],"mode":"markers","marker":{"color":"rgba(31,120,180,1)","size":11,"line":{"color":"rgba(31,120,180,1)"}},"type":"scatter","name":"Group1","textfont":{"color":"rgba(31,120,180,1)"},"error_y":{"color":"rgba(31,120,180,1)"},"error_x":{"color":"rgba(31,120,180,1)"},"line":{"color":"rgba(31,120,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0877707720041646,-0.0886764913873593],"y":[0.0786018638669425,0.0838365249211458],"text":["Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP6_2402"],"mode":"markers","marker":{"color":"rgba(51,160,44,1)","size":11,"line":{"color":"rgba(51,160,44,1)"}},"type":"scatter","name":"Group2","textfont":{"color":"rgba(51,160,44,1)"},"error_y":{"color":"rgba(51,160,44,1)"},"error_x":{"color":"rgba(51,160,44,1)"},"line":{"color":"rgba(51,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.00091609669894894,0.0057316484867194,-0.0105700763651247,0.047172829393583,0.233010886890424,0.00958335504394506,-0.00561145103017648,-0.041918969221494,-0.0194337796562408,0.194280928185142,-0.0303155168719377,-0.0546676963802766,-0.044746807733397,-0.0195570165045021,-0.0586547242329451,-0.0732566449833663,0.00036046664647183,-0.0288777729511739,-0.0142695891540695,0.186131135831846,0.000876834989257695,0.0406236917291812,-0.0335337281519618,0.000227459747879847,-0.0208875075605548,0.0310344703048944,0.0214318787600169,-0.0427067969509806,-0.0249180342121015,0.048618366078666,-0.057719681146009,0.00187366690365705,-0.01853787853239,-0.0371700515924536,-0.0331827606899027,-0.048029876281193,-0.0134260614890732,0.0100830936652151,0.0547843724192226,-0.017501684704163,-0.0307544384939403,-0.031732537662001,0.0158763249521742,0.0655583256823185,0.0147918420317416,0.116639439533855,-0.0203796247081006,-0.0381665694771089,0.153762728277716,-0.021749262006007,-0.0301064916737387,-0.0566063504433607,-0.0376841310124589,-0.00471694459757565,-0.0107051829644248,-0.00496610718356012,-0.00866138913224288,0.142782673006459,0.334377076176709,-0.0187520387181668,0.379576096673653,0.00488825905761097,-0.027011620949135,-0.0211398793427561,0.127852916658528,0.0134630806597076,-0.0202483481513357,-0.0349570714741235,-0.031699657931949,0.0448089124300323,-0.00964106400981089,-0.0273036038056401,-0.00535631199475873,-0.0374871976528879,-0.0329908221650068,0.0437539855029125,-0.0279209803485453,-0.0473582403205568,-0.055164981219857,-0.0600125424249579,-0.0505175487013237,-0.0424092491421702,-0.0238010943829872,0.155503992686468,-0.0241554683326098,0.0519061939600896,-0.00281074078588498,-0.035857990332428,0.142081103926687,0.00799813584663422,0.0202708093189836,-0.0415371147020249,-0.00360048641721173,-0.0168732968158609,-7.69682520422544e-06,-0.0148161068344636,-0.0247655694038798,0.143300275500161,-0.0264015214760327,-0.037098482286972,-0.0657298552521841,-0.041016146287813,0.00331356398351555,-0.00951519283132867,-0.0743136508329666,-0.055966080634118,-0.0497354378087559,-0.0250640118128057,0.0786741422259724,-0.0159658446013581,-0.0496383086357143,-0.0679627996973122,0.149234928006918,-0.00198864868764515,-0.0381531123703957,-0.0334478612393855,-0.00112439920417344,0.0115875991108607,-0.0344434627069261,-0.0133225708178531,-0.00218700526442425,-0.0387214209138188,-0.02964996570632,0.303291381964901,0.0652596872069457,-0.0114591684313822,0.159799013574546,-0.0339701564660326],"y":[-0.00832993190021194,0.0461352226456102,0.0311843684064767,-0.00711323617161305,-0.103196754508555,0.00321042279272735,0.0183396247096725,0.0600278694611149,0.0420759673417728,-0.10994987877627,0.0389441706127317,0.0491267840840507,0.0360451274237058,-0.0126723884186522,0.0432047507121745,0.0384809042771993,0.0272609052859637,0.0220516375645498,0.00161896729366038,-0.0953409758467297,0.0562677535888514,-0.00289907191150793,0.0336117918233875,0.00474818913842732,-0.0189412330465925,0.00777783935323498,-0.00560455446671439,0.060975523384182,0.00741404306408246,1.7415366610263e-05,-0.00115199832843433,0.0339865206138905,0.0554617018734445,-0.000873131476494058,-0.0175686036673802,0.0196149013441324,0.00302330727480495,-0.0259791690417513,-0.0411832188380884,0.0194331248163251,0.00439842344278094,0.0240947374483589,-0.0209714438007427,-0.0125690753015572,-0.00471702949172137,-0.0981872107158194,-0.00139669569176224,0.0320630821813491,-0.078225652249938,0.0596288788278913,0.0361169651139819,0.0449005976500665,0.0590855425335113,-0.00446182397694678,-0.005628370067735,0.0285645537675899,0.0176792232445731,-0.0488904358009087,-0.137771219487387,0.030517675468135,-0.131252495650675,0.0157842757729852,0.0279598636964525,0.0174809838138155,-0.0249220441003576,0.0381657192380996,-0.0194841687781322,-0.0115901617600932,0.0211981728961165,-0.0175015584639848,-0.031659056735442,-0.0123248547993839,0.024079714295579,0.0546135675390365,-0.0160091541359515,-0.00855886294764185,0.0369558406883803,0.0234958500802668,0.0379747916419227,0.0612314872870572,0.0550836623574365,0.0159277691593418,0.0151551300849081,-0.0786190167790306,-0.000560474919667497,-0.00941371188152452,0.040699116155725,0.0460348608516633,-0.0699986505516042,0.0200394164302417,0.0122954277726383,0.0254015745332463,0.0202491632468052,0.00985032036569114,-0.00338622620122989,0.0403673490600918,0.0342678631101315,-0.0547520082363424,0.0290573178707946,0.0370831251755702,0.0708594396607094,0.0162570761396918,0.00216791482621174,0.021063112672686,0.0323793313493692,0.0160070057720461,0.0193447234661261,-0.00618680322604863,-0.0574707505371626,0.0282851565723481,0.0760948703269063,0.0269621035218558,-0.0637665768189749,0.0165924921811498,0.0115518723685748,0.0574764465955154,0.0126790463220343,0.0168750821660283,0.0647548658065089,0.0249520630878216,-0.000734122349481701,0.0411354599057743,0.0357461541407483,-0.113233836869369,0.0326659165408808,0.0277455693337356,-0.055199630036227,0.0822490257367456],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"color":"rgba(128,128,128,1)","size":11,"line":{"color":"rgba(128,128,128,1)"}},"type":"scatter","name":"Ungrouped","textfont":{"color":"rgba(128,128,128,1)"},"error_y":{"color":"rgba(128,128,128,1)"},"error_x":{"color":"rgba(128,128,128,1)"},"line":{"color":"rgba(128,128,128,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-5009c207c348a575d918" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-5009c207c348a575d918">{"x":{"visdat":{"6aac1b08241d":["function () ","plotlyVisDat"]},"cur_data":"6aac1b08241d","attrs":{"6aac1b08241d":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#01665e","#7570b3","#808080","#1f78b4","#33a02c"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Reef","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0105700763651247,0.00958335504394506,-0.041918969221494,-0.0303155168719377,-0.044746807733397,-0.0586547242329451,-0.0732566449833663,-0.0427067969509806,-0.017501684704163,0.0147918420317416,-0.00496610718356012,-0.00866138913224288,0.00488825905761097,0.0134630806597076,-0.0202483481513357,-0.00964106400981089,-0.0600125424249579,-0.0238010943829872,-0.037098482286972,-0.0657298552521841,-0.00951519283132867,-0.0743136508329666,-0.055966080634118,-0.0250640118128057,0.0786741422259724,0.0652596872069457],"y":[0.0311843684064767,0.00321042279272735,0.0600278694611149,0.0389441706127317,0.0360451274237058,0.0432047507121745,0.0384809042771993,0.060975523384182,0.0194331248163251,-0.00471702949172137,0.0285645537675899,0.0176792232445731,0.0157842757729852,0.0381657192380996,-0.0194841687781322,-0.031659056735442,0.0612314872870572,0.0151551300849081,0.0370831251755702,0.0708594396607094,0.021063112672686,0.0323793313493692,0.0160070057720461,-0.00618680322604863,-0.0574707505371626,0.0326659165408808],"text":["Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP9_1121","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP9_2862","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1997","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP8_2067"],"mode":"markers","marker":{"color":"rgba(179,88,6,1)","size":11,"line":{"color":"rgba(179,88,6,1)"}},"type":"scatter","name":"HIMB","textfont":{"color":"rgba(179,88,6,1)"},"error_y":{"color":"rgba(179,88,6,1)"},"error_x":{"color":"rgba(179,88,6,1)"},"line":{"color":"rgba(179,88,6,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0406236917291812,-0.0211398793427561,-0.031699657931949,-0.0424092491421702,-0.00360048641721173,-0.0264015214760327,-0.0381531123703957,-0.0339701564660326],"y":[-0.00289907191150793,0.0174809838138155,0.0211981728961165,0.0159277691593418,0.0202491632468052,0.0290573178707946,0.0115518723685748,0.0822490257367456],"text":["Mcapitata_ATAC_TP7_1058","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP9_1306","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"color":"rgba(241,163,64,1)","size":11,"line":{"color":"rgba(241,163,64,1)"}},"type":"scatter","name":"Lilipuna.Fringe","textfont":{"color":"rgba(241,163,64,1)"},"error_y":{"color":"rgba(241,163,64,1)"},"error_x":{"color":"rgba(241,163,64,1)"},"line":{"color":"rgba(241,163,64,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0057316484867194,-0.0546676963802766,-0.057719681146009,-0.048029876281193,0.0100830936652151,-0.021749262006007,-0.0566063504433607,-0.027011620949135,0.127852916658528,-0.0349570714741235,-0.0505175487013237,-0.00281074078588498,-0.0168732968158609,-0.0148161068344636,-0.0247655694038798,0.143300275500161,-0.0679627996973122,-0.0334478612393855,-0.0344434627069261,-0.02964996570632,-0.0114591684313822,0.159799013574546],"y":[0.0461352226456102,0.0491267840840507,-0.00115199832843433,0.0196149013441324,-0.0259791690417513,0.0596288788278913,0.0449005976500665,0.0279598636964525,-0.0249220441003576,-0.0115901617600932,0.0550836623574365,0.040699116155725,0.00985032036569114,0.0403673490600918,0.0342678631101315,-0.0547520082363424,0.0269621035218558,0.0574764465955154,0.0647548658065089,0.0357461541407483,0.0277455693337356,-0.055199630036227],"text":["Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP3_1548","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP9_1467","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331"],"mode":"markers","marker":{"color":"rgba(254,224,182,1)","size":11,"line":{"color":"rgba(254,224,182,1)"}},"type":"scatter","name":"Reef.11.13","textfont":{"color":"rgba(254,224,182,1)"},"error_y":{"color":"rgba(254,224,182,1)"},"error_x":{"color":"rgba(254,224,182,1)"},"line":{"color":"rgba(254,224,182,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.29916825882356,-0.0877707720041646,0.00036046664647183,-0.0886764913873593,0.00187366690365705,-0.031732537662001,0.0158763249521742,-0.0203796247081006,-0.00471694459757565,-0.311412109391442,0.0448089124300323,-0.0273036038056401,-0.0241554683326098,0.0202708093189836,-0.0415371147020249,-7.69682520422544e-06,-0.041016146287813,-0.0159658446013581,-0.00198864868764515,-0.0387214209138188],"y":[-0.603312761641301,0.0786018638669425,0.0272609052859637,0.0838365249211458,0.0339865206138905,0.0240947374483589,-0.0209714438007427,-0.00139669569176224,-0.00446182397694678,-0.618994536873599,-0.0175015584639848,-0.0123248547993839,-0.000560474919667497,0.0122954277726383,0.0254015745332463,-0.00338622620122989,0.0162570761396918,0.0282851565723481,0.0165924921811498,0.0411354599057743],"text":["Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP6_2402","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP6_2555","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP7_2419"],"mode":"markers","marker":{"color":"rgba(216,218,235,1)","size":11,"line":{"color":"rgba(216,218,235,1)"}},"type":"scatter","name":"Reef.18","textfont":{"color":"rgba(216,218,235,1)"},"error_y":{"color":"rgba(216,218,235,1)"},"error_x":{"color":"rgba(216,218,235,1)"},"line":{"color":"rgba(216,218,235,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.00091609669894894,0.047172829393583,0.233010886890424,0.194280928185142,-0.0195570165045021,-0.0335337281519618,0.000227459747879847,-0.0208875075605548,0.0310344703048944,-0.0249180342121015,0.048618366078666,-0.01853787853239,-0.0371700515924536,-0.0331827606899027,0.0547843724192226,0.116639439533855,-0.0301064916737387,-0.0376841310124589,0.142782673006459,-0.0329908221650068,-0.0473582403205568,-0.055164981219857,0.155503992686468,0.0519061939600896,0.142081103926687,0.00331356398351555,-0.00112439920417344,0.0115875991108607],"y":[-0.00832993190021194,-0.00711323617161305,-0.103196754508555,-0.10994987877627,-0.0126723884186522,0.0336117918233875,0.00474818913842732,-0.0189412330465925,0.00777783935323498,0.00741404306408246,1.7415366610263e-05,0.0554617018734445,-0.000873131476494058,-0.0175686036673802,-0.0411832188380884,-0.0981872107158194,0.0361169651139819,0.0590855425335113,-0.0488904358009087,-0.0160091541359515,0.0234958500802668,0.0379747916419227,-0.0786190167790306,-0.00941371188152452,-0.0699986505516042,0.00216791482621174,0.0126790463220343,0.0168750821660283],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP8_1260","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1722","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317"],"mode":"markers","marker":{"color":"rgba(153,142,195,1)","size":11,"line":{"color":"rgba(153,142,195,1)"}},"type":"scatter","name":"Reef.35.36","textfont":{"color":"rgba(153,142,195,1)"},"error_y":{"color":"rgba(153,142,195,1)"},"error_x":{"color":"rgba(153,142,195,1)"},"line":{"color":"rgba(153,142,195,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.00561145103017648,-0.0194337796562408,-0.0288777729511739,-0.0142695891540695,0.186131135831846,0.000876834989257695,0.0214318787600169,-0.0134260614890732,-0.0307544384939403,0.0655583256823185,-0.0381665694771089,0.153762728277716,-0.0107051829644248,0.334377076176709,-0.0187520387181668,0.379576096673653,-0.00535631199475873,-0.0374871976528879,0.0437539855029125,-0.0279209803485453,-0.035857990332428,0.00799813584663422,-0.0497354378087559,-0.0496383086357143,0.149234928006918,-0.0133225708178531,-0.00218700526442425,0.303291381964901],"y":[0.0183396247096725,0.0420759673417728,0.0220516375645498,0.00161896729366038,-0.0953409758467297,0.0562677535888514,-0.00560455446671439,0.00302330727480495,0.00439842344278094,-0.0125690753015572,0.0320630821813491,-0.078225652249938,-0.005628370067735,-0.137771219487387,0.030517675468135,-0.131252495650675,0.024079714295579,0.0546135675390365,-0.00855886294764185,0.0369558406883803,0.0460348608516633,0.0200394164302417,0.0193447234661261,0.0760948703269063,-0.0637665768189749,0.0249520630878216,-0.000734122349481701,-0.113233836869369],"text":["Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP8_1779","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP8_1235","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP8_1246"],"mode":"markers","marker":{"color":"rgba(84,39,136,1)","size":11,"line":{"color":"rgba(84,39,136,1)"}},"type":"scatter","name":"Reef.42.43","textfont":{"color":"rgba(84,39,136,1)"},"error_y":{"color":"rgba(84,39,136,1)"},"error_x":{"color":"rgba(84,39,136,1)"},"line":{"color":"rgba(84,39,136,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-fe76de657b235835a81a" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-fe76de657b235835a81a">{"x":{"visdat":{"6aac337af9c8":["function () ","plotlyVisDat"]},"cur_data":"6aac337af9c8","attrs":{"6aac337af9c8":{"x":{},"y":{},"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#998ec3","#fee0b6","#b35806","#d8daeb","#542788","#f1a340","#01665e","#7570b3","#808080","#1f78b4","#33a02c"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Ploidy","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[0.00091609669894894,0.0057316484867194,-0.0105700763651247,0.047172829393583,0.233010886890424,0.00958335504394506,-0.29916825882356,-0.00561145103017648,-0.041918969221494,-0.0194337796562408,0.194280928185142,-0.0877707720041646,-0.0303155168719377,-0.0546676963802766,-0.044746807733397,-0.0195570165045021,-0.0586547242329451,-0.0732566449833663,0.00036046664647183,-0.0288777729511739,-0.0142695891540695,0.186131135831846,0.000876834989257695,-0.0886764913873593,0.0406236917291812,-0.0335337281519618,0.000227459747879847,-0.0208875075605548,0.0310344703048944,0.0214318787600169,-0.0427067969509806,-0.0249180342121015,0.048618366078666,-0.057719681146009,0.00187366690365705,-0.01853787853239,-0.0371700515924536,-0.0331827606899027,-0.048029876281193,-0.0134260614890732,0.0100830936652151,0.0547843724192226,-0.017501684704163,-0.0307544384939403,-0.031732537662001,0.0158763249521742,0.0655583256823185,0.0147918420317416,0.116639439533855,-0.0203796247081006,-0.0381665694771089,0.153762728277716,-0.021749262006007,-0.0301064916737387,-0.0566063504433607,-0.0376841310124589,-0.00471694459757565,-0.0107051829644248,-0.00496610718356012,-0.00866138913224288,0.142782673006459,0.334377076176709,-0.0187520387181668,0.379576096673653,0.00488825905761097,-0.027011620949135,-0.0211398793427561,0.127852916658528,0.0134630806597076,-0.0202483481513357,-0.0349570714741235,-0.311412109391442,-0.031699657931949,0.0448089124300323,-0.00964106400981089,-0.0273036038056401,-0.00535631199475873,-0.0374871976528879,-0.0329908221650068,0.0437539855029125,-0.0279209803485453,-0.0473582403205568,-0.055164981219857,-0.0600125424249579,-0.0505175487013237,-0.0424092491421702,-0.0238010943829872,0.155503992686468,-0.0241554683326098,0.0519061939600896,-0.00281074078588498,-0.035857990332428,0.142081103926687,0.00799813584663422,0.0202708093189836,-0.0415371147020249,-0.00360048641721173,-0.0168732968158609,-7.69682520422544e-06,-0.0148161068344636,-0.0247655694038798,0.143300275500161,-0.0264015214760327,-0.037098482286972,-0.0657298552521841,-0.041016146287813,0.00331356398351555,-0.00951519283132867,-0.0743136508329666,-0.055966080634118,-0.0497354378087559,-0.0250640118128057,0.0786741422259724,-0.0159658446013581,-0.0496383086357143,-0.0679627996973122,0.149234928006918,-0.00198864868764515,-0.0381531123703957,-0.0334478612393855,-0.00112439920417344,0.0115875991108607,-0.0344434627069261,-0.0133225708178531,-0.00218700526442425,-0.0387214209138188,-0.02964996570632,0.303291381964901,0.0652596872069457,-0.0114591684313822,0.159799013574546,-0.0339701564660326],"y":[-0.00832993190021194,0.0461352226456102,0.0311843684064767,-0.00711323617161305,-0.103196754508555,0.00321042279272735,-0.603312761641301,0.0183396247096725,0.0600278694611149,0.0420759673417728,-0.10994987877627,0.0786018638669425,0.0389441706127317,0.0491267840840507,0.0360451274237058,-0.0126723884186522,0.0432047507121745,0.0384809042771993,0.0272609052859637,0.0220516375645498,0.00161896729366038,-0.0953409758467297,0.0562677535888514,0.0838365249211458,-0.00289907191150793,0.0336117918233875,0.00474818913842732,-0.0189412330465925,0.00777783935323498,-0.00560455446671439,0.060975523384182,0.00741404306408246,1.7415366610263e-05,-0.00115199832843433,0.0339865206138905,0.0554617018734445,-0.000873131476494058,-0.0175686036673802,0.0196149013441324,0.00302330727480495,-0.0259791690417513,-0.0411832188380884,0.0194331248163251,0.00439842344278094,0.0240947374483589,-0.0209714438007427,-0.0125690753015572,-0.00471702949172137,-0.0981872107158194,-0.00139669569176224,0.0320630821813491,-0.078225652249938,0.0596288788278913,0.0361169651139819,0.0449005976500665,0.0590855425335113,-0.00446182397694678,-0.005628370067735,0.0285645537675899,0.0176792232445731,-0.0488904358009087,-0.137771219487387,0.030517675468135,-0.131252495650675,0.0157842757729852,0.0279598636964525,0.0174809838138155,-0.0249220441003576,0.0381657192380996,-0.0194841687781322,-0.0115901617600932,-0.618994536873599,0.0211981728961165,-0.0175015584639848,-0.031659056735442,-0.0123248547993839,0.024079714295579,0.0546135675390365,-0.0160091541359515,-0.00855886294764185,0.0369558406883803,0.0234958500802668,0.0379747916419227,0.0612314872870572,0.0550836623574365,0.0159277691593418,0.0151551300849081,-0.0786190167790306,-0.000560474919667497,-0.00941371188152452,0.040699116155725,0.0460348608516633,-0.0699986505516042,0.0200394164302417,0.0122954277726383,0.0254015745332463,0.0202491632468052,0.00985032036569114,-0.00338622620122989,0.0403673490600918,0.0342678631101315,-0.0547520082363424,0.0290573178707946,0.0370831251755702,0.0708594396607094,0.0162570761396918,0.00216791482621174,0.021063112672686,0.0323793313493692,0.0160070057720461,0.0193447234661261,-0.00618680322604863,-0.0574707505371626,0.0282851565723481,0.0760948703269063,0.0269621035218558,-0.0637665768189749,0.0165924921811498,0.0115518723685748,0.0574764465955154,0.0126790463220343,0.0168750821660283,0.0647548658065089,0.0249520630878216,-0.000734122349481701,0.0411354599057743,0.0357461541407483,-0.113233836869369,0.0326659165408808,0.0277455693337356,-0.055199630036227,0.0822490257367456],"text":["Mcapitata_ATAC_TP10_1095","Mcapitata_ATAC_TP10_1561","Mcapitata_ATAC_TP10_1631","Mcapitata_ATAC_TP1_1037","Mcapitata_ATAC_TP11_1076","Mcapitata_ATAC_TP11_1644","Mcapitata_ATAC_TP11_2302","Mcapitata_ATAC_TP1_1600","Mcapitata_ATAC_TP1_1652","Mcapitata_ATAC_TP12_1120","Mcapitata_ATAC_TP12_1452","Mcapitata_ATAC_TP12_2403","Mcapitata_ATAC_TP3_1101","Mcapitata_ATAC_TP3_1548","Mcapitata_ATAC_TP3_1628","Mcapitata_ATAC_TP4_1108","Mcapitata_ATAC_TP4_1609","Mcapitata_ATAC_TP4_1651","Mcapitata_ATAC_TP5_1196","Mcapitata_ATAC_TP5_1610","Mcapitata_ATAC_TP5_1776","Mcapitata_ATAC_TP6_1114","Mcapitata_ATAC_TP6_1611","Mcapitata_ATAC_TP6_2402","Mcapitata_ATAC_TP7_1058","Mcapitata_ATAC_TP7_1455","Mcapitata_ATAC_TP7_1499","Mcapitata_ATAC_TP8_1083","Mcapitata_ATAC_TP8_1436","Mcapitata_ATAC_TP8_1779","Mcapitata_ATAC_TP9_1121","Mcapitata_ATAC_TP9_1420","Mcapitata_ATAC_TP9_1580","Mcapitata_ATHC_TP10_1204","Mcapitata_ATHC_TP10_2554","Mcapitata_ATHC_TP10_2737","Mcapitata_ATHC_TP11_1237","Mcapitata_ATHC_TP11_2188","Mcapitata_ATHC_TP1_1218","Mcapitata_ATHC_TP11_2756","Mcapitata_ATHC_TP1_1826","Mcapitata_ATHC_TP1_2068","Mcapitata_ATHC_TP12_1154","Mcapitata_ATHC_TP12_2736","Mcapitata_ATHC_TP12_2990","Mcapitata_ATHC_TP3_1544","Mcapitata_ATHC_TP3_2731","Mcapitata_ATHC_TP3_2866","Mcapitata_ATHC_TP4_1221","Mcapitata_ATHC_TP4_2561","Mcapitata_ATHC_TP4_2734","Mcapitata_ATHC_TP5_1229","Mcapitata_ATHC_TP5_1706","Mcapitata_ATHC_TP5_2986","Mcapitata_ATHC_TP6_1212","Mcapitata_ATHC_TP6_2016","Mcapitata_ATHC_TP6_2555","Mcapitata_ATHC_TP7_1223","Mcapitata_ATHC_TP7_2860","Mcapitata_ATHC_TP7_2875","Mcapitata_ATHC_TP8_1260","Mcapitata_ATHC_TP8_2735","Mcapitata_ATHC_TP8_2753","Mcapitata_ATHC_TP9_1148","Mcapitata_ATHC_TP9_2862","Mcapitata_ATHC_TP9_2995","Mcapitata_HTAC_TP10_1315","Mcapitata_HTAC_TP10_1478","Mcapitata_HTAC_TP10_1754","Mcapitata_HTAC_TP11_1248","Mcapitata_HTAC_TP11_1562","Mcapitata_HTAC_TP11_2380","Mcapitata_HTAC_TP1_1579","Mcapitata_HTAC_TP1_2153","Mcapitata_HTAC_TP12_1632","Mcapitata_HTAC_TP12_1729","Mcapitata_HTAC_TP1_2183","Mcapitata_HTAC_TP12_2007","Mcapitata_HTAC_TP3_1289","Mcapitata_HTAC_TP3_1751","Mcapitata_HTAC_TP3_2021","Mcapitata_HTAC_TP4_1269","Mcapitata_HTAC_TP4_1481","Mcapitata_HTAC_TP4_2000","Mcapitata_HTAC_TP5_1321","Mcapitata_HTAC_TP5_1583","Mcapitata_HTAC_TP5_1997","Mcapitata_HTAC_TP6_1496","Mcapitata_HTAC_TP6_1588","Mcapitata_HTAC_TP6_1705","Mcapitata_HTAC_TP7_1278","Mcapitata_HTAC_TP7_1645","Mcapitata_HTAC_TP7_1722","Mcapitata_HTAC_TP8_1235","Mcapitata_HTAC_TP8_2386","Mcapitata_HTAC_TP8_2410","Mcapitata_HTAC_TP9_1306","Mcapitata_HTAC_TP9_1467","Mcapitata_HTAC_TP9_2412","Mcapitata_HTHC_TP10_1074","Mcapitata_HTHC_TP10_1332","Mcapitata_HTHC_TP10_1689","Mcapitata_HTHC_TP11_1178","Mcapitata_HTHC_TP11_1270","Mcapitata_HTHC_TP1_1145","Mcapitata_HTHC_TP11_2511","Mcapitata_HTHC_TP1_1323","Mcapitata_HTHC_TP1_2081","Mcapitata_HTHC_TP12_1140","Mcapitata_HTHC_TP12_1274","Mcapitata_HTHC_TP12_2190","Mcapitata_HTHC_TP3_1128","Mcapitata_HTHC_TP3_1277","Mcapitata_HTHC_TP3_2518","Mcapitata_HTHC_TP4_1124","Mcapitata_HTHC_TP4_1328","Mcapitata_HTHC_TP4_2204","Mcapitata_HTHC_TP5_1345","Mcapitata_HTHC_TP5_1449","Mcapitata_HTHC_TP5_1694","Mcapitata_HTHC_TP6_1164","Mcapitata_HTHC_TP6_1317","Mcapitata_HTHC_TP6_1604","Mcapitata_HTHC_TP7_1126","Mcapitata_HTHC_TP7_1250","Mcapitata_HTHC_TP7_2419","Mcapitata_HTHC_TP8_1082","Mcapitata_HTHC_TP8_1246","Mcapitata_HTHC_TP8_2067","Mcapitata_HTHC_TP9_1078","Mcapitata_HTHC_TP9_1331","Mcapitata_HTHC_TP9_2009"],"mode":"markers","marker":{"colorbar":{"title":"ploidy","ticklen":2},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(197,175,191,1)"],["0.0833333333333333","rgba(238,210,185,1)"],["0.125","rgba(239,189,139,1)"],["0.166666666666667","rgba(208,132,70,1)"],["0.208333333333333","rgba(185,99,33,1)"],["0.25","rgba(209,152,123,1)"],["0.291666666666667","rgba(217,207,216,1)"],["0.333333333333333","rgba(173,156,202,1)"],["0.375","rgba(119,83,161,1)"],["0.416666666666667","rgba(118,58,128,1)"],["0.458333333333333","rgba(183,108,103,1)"],["0.5","rgba(241,163,64,1)"],["0.541666666666667","rgba(161,139,80,1)"],["0.583333333333333","rgba(73,113,91,1)"],["0.625","rgba(48,106,115,1)"],["0.666666666666667","rgba(89,110,150,1)"],["0.708333333333333","rgba(119,113,175,1)"],["0.75","rgba(124,120,154,1)"],["0.791666666666667","rgba(128,127,132,1)"],["0.833333333333333","rgba(109,125,145,1)"],["0.875","rgba(73,122,167,1)"],["0.916666666666667","rgba(54,127,159,1)"],["0.958333333333333","rgba(68,143,107,1)"],["1","rgba(51,160,44,1)"]],"showscale":false,"color":[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"size":11,"line":{"colorbar":{"title":"","ticklen":2},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(197,175,191,1)"],["0.0833333333333333","rgba(238,210,185,1)"],["0.125","rgba(239,189,139,1)"],["0.166666666666667","rgba(208,132,70,1)"],["0.208333333333333","rgba(185,99,33,1)"],["0.25","rgba(209,152,123,1)"],["0.291666666666667","rgba(217,207,216,1)"],["0.333333333333333","rgba(173,156,202,1)"],["0.375","rgba(119,83,161,1)"],["0.416666666666667","rgba(118,58,128,1)"],["0.458333333333333","rgba(183,108,103,1)"],["0.5","rgba(241,163,64,1)"],["0.541666666666667","rgba(161,139,80,1)"],["0.583333333333333","rgba(73,113,91,1)"],["0.625","rgba(48,106,115,1)"],["0.666666666666667","rgba(89,110,150,1)"],["0.708333333333333","rgba(119,113,175,1)"],["0.75","rgba(124,120,154,1)"],["0.791666666666667","rgba(128,127,132,1)"],["0.833333333333333","rgba(109,125,145,1)"],["0.875","rgba(73,122,167,1)"],["0.916666666666667","rgba(54,127,159,1)"],["0.958333333333333","rgba(68,143,107,1)"],["1","rgba(51,160,44,1)"]],"showscale":false,"color":[2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,4,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}},"type":"scatter","xaxis":"x","yaxis":"y","frame":null},{"x":[-0.311412109391442,0.379576096673653],"y":[-0.618994536873599,0.0838365249211458],"type":"scatter","mode":"markers","opacity":0,"hoverinfo":"none","showlegend":false,"marker":{"colorbar":{"title":"ploidy","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"cmin":2,"cmax":4,"colorscale":[["0","rgba(153,142,195,1)"],["0.0416666666666667","rgba(197,175,191,1)"],["0.0833333333333333","rgba(238,210,185,1)"],["0.125","rgba(239,189,139,1)"],["0.166666666666667","rgba(208,132,70,1)"],["0.208333333333333","rgba(185,99,33,1)"],["0.25","rgba(209,152,123,1)"],["0.291666666666667","rgba(217,207,216,1)"],["0.333333333333333","rgba(173,156,202,1)"],["0.375","rgba(119,83,161,1)"],["0.416666666666667","rgba(118,58,128,1)"],["0.458333333333333","rgba(183,108,103,1)"],["0.5","rgba(241,163,64,1)"],["0.541666666666667","rgba(161,139,80,1)"],["0.583333333333333","rgba(73,113,91,1)"],["0.625","rgba(48,106,115,1)"],["0.666666666666667","rgba(89,110,150,1)"],["0.708333333333333","rgba(119,113,175,1)"],["0.75","rgba(124,120,154,1)"],["0.791666666666667","rgba(128,127,132,1)"],["0.833333333333333","rgba(109,125,145,1)"],["0.875","rgba(73,122,167,1)"],["0.916666666666667","rgba(54,127,159,1)"],["0.958333333333333","rgba(68,143,107,1)"],["1","rgba(51,160,44,1)"]],"showscale":true,"color":[2,4],"line":{"color":"rgba(255,127,14,1)"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
