---
title: "Plot `ANGSD` results for *P. acuta* RNA-seq samples from this study and SRA"
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
samples.info <- read.table("../../../samples_Pacuta_ALL.annotations.txt", header=T, comment.char='')
rownames(samples.info) <- samples.info$sample

samples.order <- read.table("bam.filelist.labels", header=F)
samples.info <- samples.info[samples.order$V1,]
samples.info
```

```
##                                      sample species          treatment
## Pacuta_ATAC_TP10_1159 Pacuta_ATAC_TP10_1159  Pacuta               ATAC
## Pacuta_ATAC_TP10_1559 Pacuta_ATAC_TP10_1559  Pacuta               ATAC
## Pacuta_ATAC_TP10_1641 Pacuta_ATAC_TP10_1641  Pacuta               ATAC
## Pacuta_ATAC_TP1_1043   Pacuta_ATAC_TP1_1043  Pacuta               ATAC
## Pacuta_ATAC_TP11_1103 Pacuta_ATAC_TP11_1103  Pacuta               ATAC
## Pacuta_ATAC_TP11_1777 Pacuta_ATAC_TP11_1777  Pacuta               ATAC
## Pacuta_ATAC_TP11_2306 Pacuta_ATAC_TP11_2306  Pacuta               ATAC
## Pacuta_ATAC_TP1_1775   Pacuta_ATAC_TP1_1775  Pacuta               ATAC
## Pacuta_ATAC_TP1_2363   Pacuta_ATAC_TP1_2363  Pacuta               ATAC
## Pacuta_ATAC_TP3_1041   Pacuta_ATAC_TP3_1041  Pacuta               ATAC
## Pacuta_ATAC_TP3_1471   Pacuta_ATAC_TP3_1471  Pacuta               ATAC
## Pacuta_ATAC_TP3_1637   Pacuta_ATAC_TP3_1637  Pacuta               ATAC
## Pacuta_ATAC_TP4_1060   Pacuta_ATAC_TP4_1060  Pacuta               ATAC
## Pacuta_ATAC_TP4_1762   Pacuta_ATAC_TP4_1762  Pacuta               ATAC
## Pacuta_ATAC_TP4_2002   Pacuta_ATAC_TP4_2002  Pacuta               ATAC
## Pacuta_ATAC_TP5_1059   Pacuta_ATAC_TP5_1059  Pacuta               ATAC
## Pacuta_ATAC_TP5_1563   Pacuta_ATAC_TP5_1563  Pacuta               ATAC
## Pacuta_ATAC_TP5_1757   Pacuta_ATAC_TP5_1757  Pacuta               ATAC
## Pacuta_ATAC_TP6_1050   Pacuta_ATAC_TP6_1050  Pacuta               ATAC
## Pacuta_ATAC_TP6_1468   Pacuta_ATAC_TP6_1468  Pacuta               ATAC
## Pacuta_ATAC_TP6_1542   Pacuta_ATAC_TP6_1542  Pacuta               ATAC
## Pacuta_ATAC_TP7_1047   Pacuta_ATAC_TP7_1047  Pacuta               ATAC
## Pacuta_ATAC_TP7_1445   Pacuta_ATAC_TP7_1445  Pacuta               ATAC
## Pacuta_ATAC_TP7_2413   Pacuta_ATAC_TP7_2413  Pacuta               ATAC
## Pacuta_ATAC_TP8_1051   Pacuta_ATAC_TP8_1051  Pacuta               ATAC
## Pacuta_ATAC_TP8_1755   Pacuta_ATAC_TP8_1755  Pacuta               ATAC
## Pacuta_ATAC_TP8_2012   Pacuta_ATAC_TP8_2012  Pacuta               ATAC
## Pacuta_ATAC_TP9_1141   Pacuta_ATAC_TP9_1141  Pacuta               ATAC
## Pacuta_ATAC_TP9_1594   Pacuta_ATAC_TP9_1594  Pacuta               ATAC
## Pacuta_ATAC_TP9_2357   Pacuta_ATAC_TP9_2357  Pacuta               ATAC
## Pacuta_ATHC_TP10_1205 Pacuta_ATHC_TP10_1205  Pacuta               ATHC
## Pacuta_ATHC_TP10_2197 Pacuta_ATHC_TP10_2197  Pacuta               ATHC
## Pacuta_ATHC_TP10_2550 Pacuta_ATHC_TP10_2550  Pacuta               ATHC
## Pacuta_ATHC_TP11_1147 Pacuta_ATHC_TP11_1147  Pacuta               ATHC
## Pacuta_ATHC_TP1_1207   Pacuta_ATHC_TP1_1207  Pacuta               ATHC
## Pacuta_ATHC_TP11_2668 Pacuta_ATHC_TP11_2668  Pacuta               ATHC
## Pacuta_ATHC_TP11_2879 Pacuta_ATHC_TP11_2879  Pacuta               ATHC
## Pacuta_ATHC_TP1_2743   Pacuta_ATHC_TP1_2743  Pacuta               ATHC
## Pacuta_ATHC_TP1_2977   Pacuta_ATHC_TP1_2977  Pacuta               ATHC
## Pacuta_ATHC_TP3_1219   Pacuta_ATHC_TP3_1219  Pacuta               ATHC
## Pacuta_ATHC_TP3_2534   Pacuta_ATHC_TP3_2534  Pacuta               ATHC
## Pacuta_ATHC_TP3_2750   Pacuta_ATHC_TP3_2750  Pacuta               ATHC
## Pacuta_ATHC_TP4_1220   Pacuta_ATHC_TP4_1220  Pacuta               ATHC
## Pacuta_ATHC_TP4_2733   Pacuta_ATHC_TP4_2733  Pacuta               ATHC
## Pacuta_ATHC_TP4_2993   Pacuta_ATHC_TP4_2993  Pacuta               ATHC
## Pacuta_ATHC_TP5_1296   Pacuta_ATHC_TP5_1296  Pacuta               ATHC
## Pacuta_ATHC_TP5_2212   Pacuta_ATHC_TP5_2212  Pacuta               ATHC
## Pacuta_ATHC_TP5_2877   Pacuta_ATHC_TP5_2877  Pacuta               ATHC
## Pacuta_ATHC_TP6_1254   Pacuta_ATHC_TP6_1254  Pacuta               ATHC
## Pacuta_ATHC_TP6_2870   Pacuta_ATHC_TP6_2870  Pacuta               ATHC
## Pacuta_ATHC_TP6_2999   Pacuta_ATHC_TP6_2999  Pacuta               ATHC
## Pacuta_ATHC_TP7_1281   Pacuta_ATHC_TP7_1281  Pacuta               ATHC
## Pacuta_ATHC_TP7_2409   Pacuta_ATHC_TP7_2409  Pacuta               ATHC
## Pacuta_ATHC_TP7_2878   Pacuta_ATHC_TP7_2878  Pacuta               ATHC
## Pacuta_ATHC_TP8_1459   Pacuta_ATHC_TP8_1459  Pacuta               ATHC
## Pacuta_ATHC_TP8_2564   Pacuta_ATHC_TP8_2564  Pacuta               ATHC
## Pacuta_ATHC_TP8_2861   Pacuta_ATHC_TP8_2861  Pacuta               ATHC
## Pacuta_ATHC_TP9_1451   Pacuta_ATHC_TP9_1451  Pacuta               ATHC
## Pacuta_ATHC_TP9_2873   Pacuta_ATHC_TP9_2873  Pacuta               ATHC
## Pacuta_ATHC_TP9_2979   Pacuta_ATHC_TP9_2979  Pacuta               ATHC
## Pacuta_HTAC_TP10_1225 Pacuta_HTAC_TP10_1225  Pacuta               HTAC
## Pacuta_HTAC_TP10_1536 Pacuta_HTAC_TP10_1536  Pacuta               HTAC
## Pacuta_HTAC_TP10_2064 Pacuta_HTAC_TP10_2064  Pacuta               HTAC
## Pacuta_HTAC_TP11_1582 Pacuta_HTAC_TP11_1582  Pacuta               HTAC
## Pacuta_HTAC_TP11_1596 Pacuta_HTAC_TP11_1596  Pacuta               HTAC
## Pacuta_HTAC_TP11_1647 Pacuta_HTAC_TP11_1647  Pacuta               HTAC
## Pacuta_HTAC_TP1_1653   Pacuta_HTAC_TP1_1653  Pacuta               HTAC
## Pacuta_HTAC_TP1_2005   Pacuta_HTAC_TP1_2005  Pacuta               HTAC
## Pacuta_HTAC_TP1_2414   Pacuta_HTAC_TP1_2414  Pacuta               HTAC
## Pacuta_HTAC_TP3_1617   Pacuta_HTAC_TP3_1617  Pacuta               HTAC
## Pacuta_HTAC_TP3_1642   Pacuta_HTAC_TP3_1642  Pacuta               HTAC
## Pacuta_HTAC_TP3_2026   Pacuta_HTAC_TP3_2026  Pacuta               HTAC
## Pacuta_HTAC_TP4_1581   Pacuta_HTAC_TP4_1581  Pacuta               HTAC
## Pacuta_HTAC_TP4_1701   Pacuta_HTAC_TP4_1701  Pacuta               HTAC
## Pacuta_HTAC_TP4_1767   Pacuta_HTAC_TP4_1767  Pacuta               HTAC
## Pacuta_HTAC_TP5_1303   Pacuta_HTAC_TP5_1303  Pacuta               HTAC
## Pacuta_HTAC_TP5_1571   Pacuta_HTAC_TP5_1571  Pacuta               HTAC
## Pacuta_HTAC_TP5_1707   Pacuta_HTAC_TP5_1707  Pacuta               HTAC
## Pacuta_HTAC_TP6_1330   Pacuta_HTAC_TP6_1330  Pacuta               HTAC
## Pacuta_HTAC_TP6_1466   Pacuta_HTAC_TP6_1466  Pacuta               HTAC
## Pacuta_HTAC_TP6_1744   Pacuta_HTAC_TP6_1744  Pacuta               HTAC
## Pacuta_HTAC_TP7_1487   Pacuta_HTAC_TP7_1487  Pacuta               HTAC
## Pacuta_HTAC_TP7_1728   Pacuta_HTAC_TP7_1728  Pacuta               HTAC
## Pacuta_HTAC_TP7_2072   Pacuta_HTAC_TP7_2072  Pacuta               HTAC
## Pacuta_HTAC_TP8_1329   Pacuta_HTAC_TP8_1329  Pacuta               HTAC
## Pacuta_HTAC_TP8_1765   Pacuta_HTAC_TP8_1765  Pacuta               HTAC
## Pacuta_HTAC_TP8_2513   Pacuta_HTAC_TP8_2513  Pacuta               HTAC
## Pacuta_HTAC_TP9_1302   Pacuta_HTAC_TP9_1302  Pacuta               HTAC
## Pacuta_HTAC_TP9_1486   Pacuta_HTAC_TP9_1486  Pacuta               HTAC
## Pacuta_HTAC_TP9_1696   Pacuta_HTAC_TP9_1696  Pacuta               HTAC
## Pacuta_HTHC_TP10_1238 Pacuta_HTHC_TP10_1238  Pacuta               HTHC
## Pacuta_HTHC_TP10_1732 Pacuta_HTHC_TP10_1732  Pacuta               HTHC
## Pacuta_HTHC_TP10_2300 Pacuta_HTHC_TP10_2300  Pacuta               HTHC
## Pacuta_HTHC_TP11_1416 Pacuta_HTHC_TP11_1416  Pacuta               HTHC
## Pacuta_HTHC_TP11_2185 Pacuta_HTHC_TP11_2185  Pacuta               HTHC
## Pacuta_HTHC_TP1_1239   Pacuta_HTHC_TP1_1239  Pacuta               HTHC
## Pacuta_HTHC_TP1_1676   Pacuta_HTHC_TP1_1676  Pacuta               HTHC
## Pacuta_HTHC_TP1_2210   Pacuta_HTHC_TP1_2210  Pacuta               HTHC
## Pacuta_HTHC_TP3_1227   Pacuta_HTHC_TP3_1227  Pacuta               HTHC
## Pacuta_HTHC_TP3_1418   Pacuta_HTHC_TP3_1418  Pacuta               HTHC
## Pacuta_HTHC_TP3_2527   Pacuta_HTHC_TP3_2527  Pacuta               HTHC
## Pacuta_HTHC_TP4_1169   Pacuta_HTHC_TP4_1169  Pacuta               HTHC
## Pacuta_HTHC_TP4_1343   Pacuta_HTHC_TP4_1343  Pacuta               HTHC
## Pacuta_HTHC_TP4_2195   Pacuta_HTHC_TP4_2195  Pacuta               HTHC
## Pacuta_HTHC_TP5_1168   Pacuta_HTHC_TP5_1168  Pacuta               HTHC
## Pacuta_HTHC_TP5_1415   Pacuta_HTHC_TP5_1415  Pacuta               HTHC
## Pacuta_HTHC_TP5_2087   Pacuta_HTHC_TP5_2087  Pacuta               HTHC
## Pacuta_HTHC_TP6_1138   Pacuta_HTHC_TP6_1138  Pacuta               HTHC
## Pacuta_HTHC_TP6_1595   Pacuta_HTHC_TP6_1595  Pacuta               HTHC
## Pacuta_HTHC_TP6_1721   Pacuta_HTHC_TP6_1721  Pacuta               HTHC
## Pacuta_HTHC_TP7_1090   Pacuta_HTHC_TP7_1090  Pacuta               HTHC
## Pacuta_HTHC_TP7_1427   Pacuta_HTHC_TP7_1427  Pacuta               HTHC
## Pacuta_HTHC_TP7_1820   Pacuta_HTHC_TP7_1820  Pacuta               HTHC
## Pacuta_HTHC_TP8_1184   Pacuta_HTHC_TP8_1184  Pacuta               HTHC
## Pacuta_HTHC_TP8_1709   Pacuta_HTHC_TP8_1709  Pacuta               HTHC
## Pacuta_HTHC_TP8_2304   Pacuta_HTHC_TP8_2304  Pacuta               HTHC
## Pacuta_HTHC_TP9_1131   Pacuta_HTHC_TP9_1131  Pacuta               HTHC
## Pacuta_HTHC_TP9_2202   Pacuta_HTHC_TP9_2202  Pacuta               HTHC
## Pacuta_HTHC_TP9_2305   Pacuta_HTHC_TP9_2305  Pacuta               HTHC
## SRR6914151                       SRR6914151  Pacuta      sediment-only
## SRR6914609                       SRR6914609  Pacuta      sediment-only
## SRR6914908                       SRR6914908  Pacuta          heat-only
## SRR6934388                       SRR6934388  Pacuta          heat-only
## SRR6934542                       SRR6934542  Pacuta          heat-only
## SRR6935629                       SRR6935629  Pacuta heat_with_sediment
## SRR6942678                       SRR6942678  Pacuta heat_with_sediment
## SRR6942729                       SRR6942729  Pacuta heat_with_sediment
## SRR6951423                       SRR6951423  Pacuta            control
## SRR6951744                       SRR6951744  Pacuta            control
## SRR6952431                       SRR6952431  Pacuta            control
## SRR6963586                       SRR6963586  Pacuta      sediment-only
## SRR6963878                       SRR6963878  Pacuta      not_specified
## SRR6963891                       SRR6963891  Pacuta      sediment-only
## SRR6964364                       SRR6964364  Pacuta          heat-only
## SRR6986864                       SRR6986864  Pacuta          heat-only
## SRR6987146                       SRR6987146  Pacuta          heat-only
## SRR7039808                       SRR7039808  Pacuta heat_with_sediment
## SRR7040514                       SRR7040514  Pacuta heat_with_sediment
## SRR7041301                       SRR7041301  Pacuta heat_with_sediment
## SRR7042978                       SRR7042978  Pacuta            control
## SRR7043013                       SRR7043013  Pacuta            control
## SRR7043704                       SRR7043704  Pacuta            control
## SRR7046161                       SRR7046161  Pacuta      sediment-only
## SRR7055829                       SRR7055829  Pacuta      sediment-only
## SRR7058378                       SRR7058378  Pacuta      sediment-only
## SRR7058566                       SRR7058566  Pacuta          heat-only
## SRR7058606                       SRR7058606  Pacuta          heat-only
## SRR7058616                       SRR7058616  Pacuta          heat-only
## SRR7059699                       SRR7059699  Pacuta heat_with_sediment
## SRR7062295                       SRR7062295  Pacuta heat_with_sediment
## SRR7062350                       SRR7062350  Pacuta heat_with_sediment
##                       timepoint    plugid                             reef
## Pacuta_ATAC_TP10_1159      TP10      1159                       Reef.42.43
## Pacuta_ATAC_TP10_1559      TP10      1559                  Lilipuna.Fringe
## Pacuta_ATAC_TP10_1641      TP10      1641                             HIMB
## Pacuta_ATAC_TP1_1043        TP1      1043                  Lilipuna.Fringe
## Pacuta_ATAC_TP11_1103      TP11      1103                       Reef.42.43
## Pacuta_ATAC_TP11_1777      TP11      1777                             HIMB
## Pacuta_ATAC_TP11_2306      TP11      2306                          Reef.18
## Pacuta_ATAC_TP1_1775        TP1      1775                       Reef.42.43
## Pacuta_ATAC_TP1_2363        TP1      2363                          Reef.18
## Pacuta_ATAC_TP3_1041        TP3      1041                       Reef.11.13
## Pacuta_ATAC_TP3_1471        TP3      1471                       Reef.35.36
## Pacuta_ATAC_TP3_1637        TP3      1637                       Reef.11.13
## Pacuta_ATAC_TP4_1060        TP4      1060                          Reef.18
## Pacuta_ATAC_TP4_1762        TP4      1762                       Reef.42.43
## Pacuta_ATAC_TP4_2002        TP4      2002                       Reef.42.43
## Pacuta_ATAC_TP5_1059        TP5      1059                       Reef.11.13
## Pacuta_ATAC_TP5_1563        TP5      1563                  Lilipuna.Fringe
## Pacuta_ATAC_TP5_1757        TP5      1757                       Reef.42.43
## Pacuta_ATAC_TP6_1050        TP6      1050                          Reef.18
## Pacuta_ATAC_TP6_1468        TP6      1468                       Reef.35.36
## Pacuta_ATAC_TP6_1542        TP6      1542                  Lilipuna.Fringe
## Pacuta_ATAC_TP7_1047        TP7      1047                  Lilipuna.Fringe
## Pacuta_ATAC_TP7_1445        TP7      1445                       Reef.11.13
## Pacuta_ATAC_TP7_2413        TP7      2413                          Reef.18
## Pacuta_ATAC_TP8_1051        TP8      1051                       Reef.11.13
## Pacuta_ATAC_TP8_1755        TP8      1755                       Reef.42.43
## Pacuta_ATAC_TP8_2012        TP8      2012                       Reef.42.43
## Pacuta_ATAC_TP9_1141        TP9      1141                       Reef.42.43
## Pacuta_ATAC_TP9_1594        TP9      1594                       Reef.35.36
## Pacuta_ATAC_TP9_2357        TP9      2357                          Reef.18
## Pacuta_ATHC_TP10_1205      TP10      1205                       Reef.35.36
## Pacuta_ATHC_TP10_2197      TP10      2197                       Reef.35.36
## Pacuta_ATHC_TP10_2550      TP10      2550                          Reef.18
## Pacuta_ATHC_TP11_1147      TP11      1147                             HIMB
## Pacuta_ATHC_TP1_1207        TP1      1207                       Reef.11.13
## Pacuta_ATHC_TP11_2668      TP11      2668                  Lilipuna.Fringe
## Pacuta_ATHC_TP11_2879      TP11      2879                             HIMB
## Pacuta_ATHC_TP1_2743        TP1      2743                  Lilipuna.Fringe
## Pacuta_ATHC_TP1_2977        TP1      2977                       Reef.11.13
## Pacuta_ATHC_TP3_1219        TP3      1219                       Reef.35.36
## Pacuta_ATHC_TP3_2534        TP3      2534                          Reef.18
## Pacuta_ATHC_TP3_2750        TP3      2750                       Reef.42.43
## Pacuta_ATHC_TP4_1220        TP4      1220                  Lilipuna.Fringe
## Pacuta_ATHC_TP4_2733        TP4      2733                       Reef.42.43
## Pacuta_ATHC_TP4_2993        TP4      2993                       Reef.11.13
## Pacuta_ATHC_TP5_1296        TP5      1296                       Reef.11.13
## Pacuta_ATHC_TP5_2212        TP5      2212                       Reef.35.36
## Pacuta_ATHC_TP5_2877        TP5      2877                             HIMB
## Pacuta_ATHC_TP6_1254        TP6      1254                       Reef.42.43
## Pacuta_ATHC_TP6_2870        TP6      2870                             HIMB
## Pacuta_ATHC_TP6_2999        TP6      2999                  Lilipuna.Fringe
## Pacuta_ATHC_TP7_1281        TP7      1281                       Reef.35.36
## Pacuta_ATHC_TP7_2409        TP7      2409                  Lilipuna.Fringe
## Pacuta_ATHC_TP7_2878        TP7      2878                             HIMB
## Pacuta_ATHC_TP8_1459        TP8      1459                          Reef.18
## Pacuta_ATHC_TP8_2564        TP8      2564                          Reef.18
## Pacuta_ATHC_TP8_2861        TP8      2861                             HIMB
## Pacuta_ATHC_TP9_1451        TP9      1451                          Reef.18
## Pacuta_ATHC_TP9_2873        TP9      2873                             HIMB
## Pacuta_ATHC_TP9_2979        TP9      2979                       Reef.35.36
## Pacuta_HTAC_TP10_1225      TP10      1225                             HIMB
## Pacuta_HTAC_TP10_1536      TP10      1536                       Reef.35.36
## Pacuta_HTAC_TP10_2064      TP10      2064                             HIMB
## Pacuta_HTAC_TP11_1582      TP11      1582                  Lilipuna.Fringe
## Pacuta_HTAC_TP11_1596      TP11      1596                       Reef.42.43
## Pacuta_HTAC_TP11_1647      TP11      1647                             HIMB
## Pacuta_HTAC_TP1_1653        TP1      1653                             HIMB
## Pacuta_HTAC_TP1_2005        TP1      2005                          Reef.18
## Pacuta_HTAC_TP1_2414        TP1      2414                          Reef.18
## Pacuta_HTAC_TP3_1617        TP3      1617                       Reef.42.43
## Pacuta_HTAC_TP3_1642        TP3      1642                             HIMB
## Pacuta_HTAC_TP3_2026        TP3      2026                       Reef.42.43
## Pacuta_HTAC_TP4_1581        TP4      1581                  Lilipuna.Fringe
## Pacuta_HTAC_TP4_1701        TP4      1701                          Reef.18
## Pacuta_HTAC_TP4_1767        TP4      1767                       Reef.42.43
## Pacuta_HTAC_TP5_1303        TP5      1303                       Reef.35.36
## Pacuta_HTAC_TP5_1571        TP5      1571                       Reef.35.36
## Pacuta_HTAC_TP5_1707        TP5      1707                       Reef.11.13
## Pacuta_HTAC_TP6_1330        TP6      1330                       Reef.35.36
## Pacuta_HTAC_TP6_1466        TP6      1466                  Lilipuna.Fringe
## Pacuta_HTAC_TP6_1744        TP6      1744                       Reef.35.36
## Pacuta_HTAC_TP7_1487        TP7      1487                  Lilipuna.Fringe
## Pacuta_HTAC_TP7_1728        TP7      1728                          Reef.18
## Pacuta_HTAC_TP7_2072        TP7      2072                             HIMB
## Pacuta_HTAC_TP8_1329        TP8      1329                  Lilipuna.Fringe
## Pacuta_HTAC_TP8_1765        TP8      1765                       Reef.42.43
## Pacuta_HTAC_TP8_2513        TP8      2513                          Reef.18
## Pacuta_HTAC_TP9_1302        TP9      1302                       Reef.35.36
## Pacuta_HTAC_TP9_1486        TP9      1486                       Reef.11.13
## Pacuta_HTAC_TP9_1696        TP9      1696                  Lilipuna.Fringe
## Pacuta_HTHC_TP10_1238      TP10      1238                       Reef.42.43
## Pacuta_HTHC_TP10_1732      TP10      1732                  Lilipuna.Fringe
## Pacuta_HTHC_TP10_2300      TP10      2300                          Reef.18
## Pacuta_HTHC_TP11_1416      TP11      1416                       Reef.11.13
## Pacuta_HTHC_TP11_2185      TP11      2185                             HIMB
## Pacuta_HTHC_TP1_1239        TP1      1239                       Reef.42.43
## Pacuta_HTHC_TP1_1676        TP1      1676                       Reef.42.43
## Pacuta_HTHC_TP1_2210        TP1      2210                       Reef.42.43
## Pacuta_HTHC_TP3_1227        TP3      1227                       Reef.42.43
## Pacuta_HTHC_TP3_1418        TP3      1418                  Lilipuna.Fringe
## Pacuta_HTHC_TP3_2527        TP3      2527                          Reef.18
## Pacuta_HTHC_TP4_1169        TP4      1169                       Reef.35.36
## Pacuta_HTHC_TP4_1343        TP4      1343                  Lilipuna.Fringe
## Pacuta_HTHC_TP4_2195        TP4      2195                       Reef.42.43
## Pacuta_HTHC_TP5_1168        TP5      1168                       Reef.35.36
## Pacuta_HTHC_TP5_1415        TP5      1415                       Reef.11.13
## Pacuta_HTHC_TP5_2087        TP5      2087                             HIMB
## Pacuta_HTHC_TP6_1138        TP6      1138                             HIMB
## Pacuta_HTHC_TP6_1595        TP6      1595                             HIMB
## Pacuta_HTHC_TP6_1721        TP6      1721                  Lilipuna.Fringe
## Pacuta_HTHC_TP7_1090        TP7      1090                  Lilipuna.Fringe
## Pacuta_HTHC_TP7_1427        TP7      1427                  Lilipuna.Fringe
## Pacuta_HTHC_TP7_1820        TP7      1820                       Reef.11.13
## Pacuta_HTHC_TP8_1184        TP8      1184                  Lilipuna.Fringe
## Pacuta_HTHC_TP8_1709        TP8      1709                  Lilipuna.Fringe
## Pacuta_HTHC_TP8_2304        TP8      2304                          Reef.18
## Pacuta_HTHC_TP9_1131        TP9      1131                             HIMB
## Pacuta_HTHC_TP9_2202        TP9      2202                             HIMB
## Pacuta_HTHC_TP9_2305        TP9      2305                          Reef.18
## SRR6914151                 <NA>    RAS1_3 SRA-Singapore-Raffles_Lighthouse
## SRR6914609                 <NA>    RAS1_2 SRA-Singapore-Raffles_Lighthouse
## SRR6914908                 <NA>     RH1_1 SRA-Singapore-Raffles_Lighthouse
## SRR6934388                 <NA>     RH1_2 SRA-Singapore-Raffles_Lighthouse
## SRR6934542                 <NA>     RH1_3 SRA-Singapore-Raffles_Lighthouse
## SRR6935629                 <NA>    RHS1_1 SRA-Singapore-Raffles_Lighthouse
## SRR6942678                 <NA>    RHS1_2 SRA-Singapore-Raffles_Lighthouse
## SRR6942729                 <NA>    RHS1_3 SRA-Singapore-Raffles_Lighthouse
## SRR6951423                 <NA>     KA1_1        SRA-Singapore-Kusu_Island
## SRR6951744                 <NA>     KA1_2        SRA-Singapore-Kusu_Island
## SRR6952431                 <NA>     KA1_3        SRA-Singapore-Kusu_Island
## SRR6963586                 <NA>    KAS1_1        SRA-Singapore-Kusu_Island
## SRR6963878                 <NA>    KAS1_2        SRA-Singapore-Kusu_Island
## SRR6963891                 <NA>    KAS1_3        SRA-Singapore-Kusu_Island
## SRR6964364                 <NA>     KH1_1        SRA-Singapore-Kusu_Island
## SRR6986864                 <NA>     KH1_2        SRA-Singapore-Kusu_Island
## SRR6987146                 <NA>     KH1_3        SRA-Singapore-Kusu_Island
## SRR7039808                 <NA>    KHS1_1        SRA-Singapore-Kusu_Island
## SRR7040514                 <NA>    KHS1_2        SRA-Singapore-Kusu_Island
## SRR7041301                 <NA>    KHS1_3        SRA-Singapore-Kusu_Island
## SRR7042978                 <NA>     SA1_1    SRA-Singapore-St_Johns_Island
## SRR7043013                 <NA>     SA1_2    SRA-Singapore-St_Johns_Island
## SRR7043704                 <NA>     SA1_3    SRA-Singapore-St_Johns_Island
## SRR7046161                 <NA>    SAS1_1    SRA-Singapore-St_Johns_Island
## SRR7055829                 <NA>    SAS1_2    SRA-Singapore-St_Johns_Island
## SRR7058378                 <NA>    SAS1_3    SRA-Singapore-St_Johns_Island
## SRR7058566                 <NA>     SH1_1    SRA-Singapore-St_Johns_Island
## SRR7058606                 <NA>     SH1_2    SRA-Singapore-St_Johns_Island
## SRR7058616                 <NA>     SH1_3    SRA-Singapore-St_Johns_Island
## SRR7059699                 <NA>    SHS1_1    SRA-Singapore-St_Johns_Island
## SRR7062295                 <NA>    SHS1_2    SRA-Singapore-St_Johns_Island
## SRR7062350                 <NA> SHS1_3_B2    SRA-Singapore-St_Johns_Island
##                       reef_color ploidy ploidy_color
## Pacuta_ATAC_TP10_1159    #542788      2      #1b9e77
## Pacuta_ATAC_TP10_1559    #f1a340      3      #d95f02
## Pacuta_ATAC_TP10_1641    #b35806      3      #d95f02
## Pacuta_ATAC_TP1_1043     #f1a340      3      #d95f02
## Pacuta_ATAC_TP11_1103    #542788      3      #d95f02
## Pacuta_ATAC_TP11_1777    #b35806      3      #d95f02
## Pacuta_ATAC_TP11_2306    #d8daeb      3      #d95f02
## Pacuta_ATAC_TP1_1775     #542788      2      #1b9e77
## Pacuta_ATAC_TP1_2363     #d8daeb      3      #d95f02
## Pacuta_ATAC_TP3_1041     #fee0b6      3      #d95f02
## Pacuta_ATAC_TP3_1471     #998ec3      3      #d95f02
## Pacuta_ATAC_TP3_1637     #fee0b6      3      #d95f02
## Pacuta_ATAC_TP4_1060     #d8daeb      3      #d95f02
## Pacuta_ATAC_TP4_1762     #542788      3      #d95f02
## Pacuta_ATAC_TP4_2002     #542788      3      #d95f02
## Pacuta_ATAC_TP5_1059     #fee0b6      2      #1b9e77
## Pacuta_ATAC_TP5_1563     #f1a340      3      #d95f02
## Pacuta_ATAC_TP5_1757     #542788      3      #d95f02
## Pacuta_ATAC_TP6_1050     #d8daeb      2      #1b9e77
## Pacuta_ATAC_TP6_1468     #998ec3      2      #1b9e77
## Pacuta_ATAC_TP6_1542     #f1a340      3      #d95f02
## Pacuta_ATAC_TP7_1047     #f1a340      2      #1b9e77
## Pacuta_ATAC_TP7_1445     #fee0b6      2      #1b9e77
## Pacuta_ATAC_TP7_2413     #d8daeb      3      #d95f02
## Pacuta_ATAC_TP8_1051     #fee0b6      3      #d95f02
## Pacuta_ATAC_TP8_1755     #542788      2      #1b9e77
## Pacuta_ATAC_TP8_2012     #542788      3      #d95f02
## Pacuta_ATAC_TP9_1141     #542788      2      #1b9e77
## Pacuta_ATAC_TP9_1594     #998ec3      3      #d95f02
## Pacuta_ATAC_TP9_2357     #d8daeb      3      #d95f02
## Pacuta_ATHC_TP10_1205    #998ec3      2      #1b9e77
## Pacuta_ATHC_TP10_2197    #998ec3      2      #1b9e77
## Pacuta_ATHC_TP10_2550    #d8daeb      2      #1b9e77
## Pacuta_ATHC_TP11_1147    #b35806      3      #d95f02
## Pacuta_ATHC_TP1_1207     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP11_2668    #f1a340      2      #1b9e77
## Pacuta_ATHC_TP11_2879    #b35806      2      #1b9e77
## Pacuta_ATHC_TP1_2743     #f1a340      3      #d95f02
## Pacuta_ATHC_TP1_2977     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP3_1219     #998ec3      2      #1b9e77
## Pacuta_ATHC_TP3_2534     #d8daeb      3      #d95f02
## Pacuta_ATHC_TP3_2750     #542788      3      #d95f02
## Pacuta_ATHC_TP4_1220     #f1a340      3      #d95f02
## Pacuta_ATHC_TP4_2733     #542788      3      #d95f02
## Pacuta_ATHC_TP4_2993     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP5_1296     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP5_2212     #998ec3      2      #1b9e77
## Pacuta_ATHC_TP5_2877     #b35806      3      #d95f02
## Pacuta_ATHC_TP6_1254     #542788      3      #d95f02
## Pacuta_ATHC_TP6_2870     #b35806      3      #d95f02
## Pacuta_ATHC_TP6_2999     #f1a340      2      #1b9e77
## Pacuta_ATHC_TP7_1281     #998ec3      2      #1b9e77
## Pacuta_ATHC_TP7_2409     #f1a340      3      #d95f02
## Pacuta_ATHC_TP7_2878     #b35806      3      #d95f02
## Pacuta_ATHC_TP8_1459     #d8daeb      3      #d95f02
## Pacuta_ATHC_TP8_2564     #d8daeb      3      #d95f02
## Pacuta_ATHC_TP8_2861     #b35806      2      #1b9e77
## Pacuta_ATHC_TP9_1451     #d8daeb      3      #d95f02
## Pacuta_ATHC_TP9_2873     #b35806      3      #d95f02
## Pacuta_ATHC_TP9_2979     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP10_1225    #b35806      2      #1b9e77
## Pacuta_HTAC_TP10_1536    #998ec3      3      #d95f02
## Pacuta_HTAC_TP10_2064    #b35806      3      #d95f02
## Pacuta_HTAC_TP11_1582    #f1a340      3      #d95f02
## Pacuta_HTAC_TP11_1596    #542788      3      #d95f02
## Pacuta_HTAC_TP11_1647    #b35806      3      #d95f02
## Pacuta_HTAC_TP1_1653     #b35806      2      #1b9e77
## Pacuta_HTAC_TP1_2005     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP1_2414     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP3_1617     #542788      3      #d95f02
## Pacuta_HTAC_TP3_1642     #b35806      3      #d95f02
## Pacuta_HTAC_TP3_2026     #542788      2      #1b9e77
## Pacuta_HTAC_TP4_1581     #f1a340      3      #d95f02
## Pacuta_HTAC_TP4_1701     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP4_1767     #542788      3      #d95f02
## Pacuta_HTAC_TP5_1303     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP5_1571     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP5_1707     #fee0b6      3      #d95f02
## Pacuta_HTAC_TP6_1330     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP6_1466     #f1a340      3      #d95f02
## Pacuta_HTAC_TP6_1744     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP7_1487     #f1a340      2      #1b9e77
## Pacuta_HTAC_TP7_1728     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP7_2072     #b35806      3      #d95f02
## Pacuta_HTAC_TP8_1329     #f1a340      2      #1b9e77
## Pacuta_HTAC_TP8_1765     #542788      3      #d95f02
## Pacuta_HTAC_TP8_2513     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP9_1302     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP9_1486     #fee0b6      2      #1b9e77
## Pacuta_HTAC_TP9_1696     #f1a340      3      #d95f02
## Pacuta_HTHC_TP10_1238    #542788      3      #d95f02
## Pacuta_HTHC_TP10_1732    #f1a340      3      #d95f02
## Pacuta_HTHC_TP10_2300    #d8daeb      3      #d95f02
## Pacuta_HTHC_TP11_1416    #fee0b6      2      #1b9e77
## Pacuta_HTHC_TP11_2185    #b35806      3      #d95f02
## Pacuta_HTHC_TP1_1239     #542788      3      #d95f02
## Pacuta_HTHC_TP1_1676     #542788      2      #1b9e77
## Pacuta_HTHC_TP1_2210     #542788      2      #1b9e77
## Pacuta_HTHC_TP3_1227     #542788      3      #d95f02
## Pacuta_HTHC_TP3_1418     #f1a340      2      #1b9e77
## Pacuta_HTHC_TP3_2527     #d8daeb      2      #1b9e77
## Pacuta_HTHC_TP4_1169     #998ec3      2      #1b9e77
## Pacuta_HTHC_TP4_1343     #f1a340      3      #d95f02
## Pacuta_HTHC_TP4_2195     #542788      3      #d95f02
## Pacuta_HTHC_TP5_1168     #998ec3      2      #1b9e77
## Pacuta_HTHC_TP5_1415     #fee0b6      2      #1b9e77
## Pacuta_HTHC_TP5_2087     #b35806      3      #d95f02
## Pacuta_HTHC_TP6_1138     #b35806      3      #d95f02
## Pacuta_HTHC_TP6_1595     #b35806      3      #d95f02
## Pacuta_HTHC_TP6_1721     #f1a340      2      #1b9e77
## Pacuta_HTHC_TP7_1090     #f1a340      3      #d95f02
## Pacuta_HTHC_TP7_1427     #f1a340      3      #d95f02
## Pacuta_HTHC_TP7_1820     #fee0b6      3      #d95f02
## Pacuta_HTHC_TP8_1184     #f1a340      3      #d95f02
## Pacuta_HTHC_TP8_1709     #f1a340      3      #d95f02
## Pacuta_HTHC_TP8_2304     #d8daeb      3      #d95f02
## Pacuta_HTHC_TP9_1131     #b35806      3      #d95f02
## Pacuta_HTHC_TP9_2202     #b35806      3      #d95f02
## Pacuta_HTHC_TP9_2305     #d8daeb      3      #d95f02
## SRR6914151               #9d529e      2      #01665e
## SRR6914609               #9d529e      2      #01665e
## SRR6914908               #9d529e      2      #01665e
## SRR6934388               #9d529e      2      #01665e
## SRR6934542               #9d529e      2      #01665e
## SRR6935629               #9d529e      2      #01665e
## SRR6942678               #9d529e      2      #01665e
## SRR6942729               #9d529e      2      #01665e
## SRR6951423               #bc80bd      2      #01665e
## SRR6951744               #bc80bd      2      #01665e
## SRR6952431               #bc80bd      2      #01665e
## SRR6963586               #bc80bd      2      #01665e
## SRR6963878               #bc80bd      2      #01665e
## SRR6963891               #bc80bd      2      #01665e
## SRR6964364               #bc80bd      2      #01665e
## SRR6986864               #bc80bd      2      #01665e
## SRR6987146               #bc80bd      2      #01665e
## SRR7039808               #bc80bd      2      #01665e
## SRR7040514               #bc80bd      2      #01665e
## SRR7041301               #bc80bd      2      #01665e
## SRR7042978               #e1cbe1      2      #01665e
## SRR7043013               #e1cbe1      2      #01665e
## SRR7043704               #e1cbe1      2      #01665e
## SRR7046161               #e1cbe1      2      #01665e
## SRR7055829               #e1cbe1      2      #01665e
## SRR7058378               #e1cbe1      2      #01665e
## SRR7058566               #e1cbe1      2      #01665e
## SRR7058606               #e1cbe1      2      #01665e
## SRR7058616               #e1cbe1      2      #01665e
## SRR7059699               #e1cbe1      2      #01665e
## SRR7062295               #e1cbe1      2      #01665e
## SRR7062350               #e1cbe1      2      #01665e
##                                                  group group_color
## Pacuta_ATAC_TP10_1159                           Group5     #ff7f00
## Pacuta_ATAC_TP10_1559                           Group3     #c51b7d
## Pacuta_ATAC_TP10_1641                           Group3     #c51b7d
## Pacuta_ATAC_TP1_1043                            Group4     #6a3d9a
## Pacuta_ATAC_TP11_1103                           Group2     #33a02c
## Pacuta_ATAC_TP11_1777                           Group1     #1f78b4
## Pacuta_ATAC_TP11_2306                           Group2     #33a02c
## Pacuta_ATAC_TP1_1775                            Group5     #ff7f00
## Pacuta_ATAC_TP1_2363                            Group3     #c51b7d
## Pacuta_ATAC_TP3_1041                            Group2     #33a02c
## Pacuta_ATAC_TP3_1471                            Group3     #c51b7d
## Pacuta_ATAC_TP3_1637                            Group2     #33a02c
## Pacuta_ATAC_TP4_1060                            Group3     #c51b7d
## Pacuta_ATAC_TP4_1762                            Group2     #33a02c
## Pacuta_ATAC_TP4_2002                            Group3     #c51b7d
## Pacuta_ATAC_TP5_1059                            Group7     #b15928
## Pacuta_ATAC_TP5_1563                            Group3     #c51b7d
## Pacuta_ATAC_TP5_1757                            Group3     #c51b7d
## Pacuta_ATAC_TP6_1050                            Group6     #e31a1c
## Pacuta_ATAC_TP6_1468                            Group5     #ff7f00
## Pacuta_ATAC_TP6_1542                            Group4     #6a3d9a
## Pacuta_ATAC_TP7_1047                            Group6     #e31a1c
## Pacuta_ATAC_TP7_1445                           Ungroup     #808080
## Pacuta_ATAC_TP7_2413                            Group2     #33a02c
## Pacuta_ATAC_TP8_1051                            Group1     #1f78b4
## Pacuta_ATAC_TP8_1755                            Group5     #ff7f00
## Pacuta_ATAC_TP8_2012                            Group2     #33a02c
## Pacuta_ATAC_TP9_1141                            Group5     #ff7f00
## Pacuta_ATAC_TP9_1594                            Group3     #c51b7d
## Pacuta_ATAC_TP9_2357                            Group2     #33a02c
## Pacuta_ATHC_TP10_1205                          Ungroup     #808080
## Pacuta_ATHC_TP10_2197                           Group6     #e31a1c
## Pacuta_ATHC_TP10_2550                          Ungroup     #808080
## Pacuta_ATHC_TP11_1147                           Group2     #33a02c
## Pacuta_ATHC_TP1_1207                            Group6     #e31a1c
## Pacuta_ATHC_TP11_2668                           Group6     #e31a1c
## Pacuta_ATHC_TP11_2879                           Group6     #e31a1c
## Pacuta_ATHC_TP1_2743                            Group4     #6a3d9a
## Pacuta_ATHC_TP1_2977                            Group6     #e31a1c
## Pacuta_ATHC_TP3_1219                            Group6     #e31a1c
## Pacuta_ATHC_TP3_2534                            Group2     #33a02c
## Pacuta_ATHC_TP3_2750                            Group4     #6a3d9a
## Pacuta_ATHC_TP4_1220                            Group2     #33a02c
## Pacuta_ATHC_TP4_2733                            Group1     #1f78b4
## Pacuta_ATHC_TP4_2993                            Group6     #e31a1c
## Pacuta_ATHC_TP5_1296                            Group6     #e31a1c
## Pacuta_ATHC_TP5_2212                            Group5     #ff7f00
## Pacuta_ATHC_TP5_2877                            Group2     #33a02c
## Pacuta_ATHC_TP6_1254                           Ungroup     #808080
## Pacuta_ATHC_TP6_2870                            Group2     #33a02c
## Pacuta_ATHC_TP6_2999                            Group6     #e31a1c
## Pacuta_ATHC_TP7_1281                            Group6     #e31a1c
## Pacuta_ATHC_TP7_2409                            Group3     #c51b7d
## Pacuta_ATHC_TP7_2878                            Group2     #33a02c
## Pacuta_ATHC_TP8_1459                            Group1     #1f78b4
## Pacuta_ATHC_TP8_2564                            Group3     #c51b7d
## Pacuta_ATHC_TP8_2861                            Group6     #e31a1c
## Pacuta_ATHC_TP9_1451                            Group2     #33a02c
## Pacuta_ATHC_TP9_2873                            Group2     #33a02c
## Pacuta_ATHC_TP9_2979                            Group6     #e31a1c
## Pacuta_HTAC_TP10_1225                           Group6     #e31a1c
## Pacuta_HTAC_TP10_1536                           Group3     #c51b7d
## Pacuta_HTAC_TP10_2064                           Group2     #33a02c
## Pacuta_HTAC_TP11_1582                           Group3     #c51b7d
## Pacuta_HTAC_TP11_1596                           Group2     #33a02c
## Pacuta_HTAC_TP11_1647                           Group3     #c51b7d
## Pacuta_HTAC_TP1_1653                            Group8     #000000
## Pacuta_HTAC_TP1_2005                            Group3     #c51b7d
## Pacuta_HTAC_TP1_2414                            Group2     #33a02c
## Pacuta_HTAC_TP3_1617                            Group3     #c51b7d
## Pacuta_HTAC_TP3_1642                            Group1     #1f78b4
## Pacuta_HTAC_TP3_2026                            Group6     #e31a1c
## Pacuta_HTAC_TP4_1581                            Group4     #6a3d9a
## Pacuta_HTAC_TP4_1701                            Group3     #c51b7d
## Pacuta_HTAC_TP4_1767                            Group1     #1f78b4
## Pacuta_HTAC_TP5_1303                            Group6     #e31a1c
## Pacuta_HTAC_TP5_1571                            Group6     #e31a1c
## Pacuta_HTAC_TP5_1707                            Group2     #33a02c
## Pacuta_HTAC_TP6_1330                            Group6     #e31a1c
## Pacuta_HTAC_TP6_1466                            Group2     #33a02c
## Pacuta_HTAC_TP6_1744                            Group5     #ff7f00
## Pacuta_HTAC_TP7_1487                            Group6     #e31a1c
## Pacuta_HTAC_TP7_1728                            Group2     #33a02c
## Pacuta_HTAC_TP7_2072                            Group1     #1f78b4
## Pacuta_HTAC_TP8_1329                            Group6     #e31a1c
## Pacuta_HTAC_TP8_1765                            Group3     #c51b7d
## Pacuta_HTAC_TP8_2513                            Group2     #33a02c
## Pacuta_HTAC_TP9_1302                            Group5     #ff7f00
## Pacuta_HTAC_TP9_1486                            Group5     #ff7f00
## Pacuta_HTAC_TP9_1696                            Group2     #33a02c
## Pacuta_HTHC_TP10_1238                           Group3     #c51b7d
## Pacuta_HTHC_TP10_1732                           Group3     #c51b7d
## Pacuta_HTHC_TP10_2300                           Group2     #33a02c
## Pacuta_HTHC_TP11_1416                           Group6     #e31a1c
## Pacuta_HTHC_TP11_2185                          Ungroup     #808080
## Pacuta_HTHC_TP1_1239                            Group3     #c51b7d
## Pacuta_HTHC_TP1_1676                            Group6     #e31a1c
## Pacuta_HTHC_TP1_2210                            Group7     #b15928
## Pacuta_HTHC_TP3_1227                            Group2     #33a02c
## Pacuta_HTHC_TP3_1418                            Group6     #e31a1c
## Pacuta_HTHC_TP3_2527                            Group6     #e31a1c
## Pacuta_HTHC_TP4_1169                            Group6     #e31a1c
## Pacuta_HTHC_TP4_1343                            Group4     #6a3d9a
## Pacuta_HTHC_TP4_2195                            Group3     #c51b7d
## Pacuta_HTHC_TP5_1168                            Group6     #e31a1c
## Pacuta_HTHC_TP5_1415                           Ungroup     #808080
## Pacuta_HTHC_TP5_2087                            Group2     #33a02c
## Pacuta_HTHC_TP6_1138                            Group2     #33a02c
## Pacuta_HTHC_TP6_1595                            Group2     #33a02c
## Pacuta_HTHC_TP6_1721                            Group8     #000000
## Pacuta_HTHC_TP7_1090                            Group3     #c51b7d
## Pacuta_HTHC_TP7_1427                            Group4     #6a3d9a
## Pacuta_HTHC_TP7_1820                            Group3     #c51b7d
## Pacuta_HTHC_TP8_1184                            Group4     #6a3d9a
## Pacuta_HTHC_TP8_1709                            Group2     #33a02c
## Pacuta_HTHC_TP8_2304                            Group2     #33a02c
## Pacuta_HTHC_TP9_1131                            Group2     #33a02c
## Pacuta_HTHC_TP9_2202                            Group2     #33a02c
## Pacuta_HTHC_TP9_2305                            Group2     #33a02c
## SRR6914151            SRA-Singapore-Raffles_Lighthouse     #9d529e
## SRR6914609            SRA-Singapore-Raffles_Lighthouse     #9d529e
## SRR6914908            SRA-Singapore-Raffles_Lighthouse     #9d529e
## SRR6934388            SRA-Singapore-Raffles_Lighthouse     #9d529e
## SRR6934542            SRA-Singapore-Raffles_Lighthouse     #9d529e
## SRR6935629            SRA-Singapore-Raffles_Lighthouse     #9d529e
## SRR6942678            SRA-Singapore-Raffles_Lighthouse     #9d529e
## SRR6942729            SRA-Singapore-Raffles_Lighthouse     #9d529e
## SRR6951423                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR6951744                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR6952431                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR6963586                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR6963878                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR6963891                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR6964364                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR6986864                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR6987146                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR7039808                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR7040514                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR7041301                   SRA-Singapore-Kusu_Island     #bc80bd
## SRR7042978               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7043013               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7043704               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7046161               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7055829               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7058378               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7058566               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7058606               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7058616               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7059699               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7062295               SRA-Singapore-St_Johns_Island     #e1cbe1
## SRR7062350               SRA-Singapore-St_Johns_Island     #e1cbe1
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
q  <- read.table("PCAngsd.angsd.beagle.gz.Admixture.admix.7.Q")
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
##                       Reef.42.43                  Lilipuna.Fringe 
##                        "#542788"                        "#f1a340" 
##                             HIMB                          Reef.18 
##                        "#b35806"                        "#d8daeb" 
##                       Reef.11.13                       Reef.35.36 
##                        "#fee0b6"                        "#998ec3" 
## SRA-Singapore-Raffles_Lighthouse        SRA-Singapore-Kusu_Island 
##                        "#9d529e"                        "#bc80bd" 
##    SRA-Singapore-St_Johns_Island                                2 
##                        "#e1cbe1"                        "#1b9e77" 
##                                3                                2 
##                        "#d95f02"                        "#01665e" 
##                           Group5                           Group3 
##                        "#ff7f00"                        "#c51b7d" 
##                           Group4                           Group2 
##                        "#6a3d9a"                        "#33a02c" 
##                           Group1                           Group7 
##                        "#1f78b4"                        "#b15928" 
##                           Group6                          Ungroup 
##                        "#e31a1c"                        "#808080" 
##                           Group8 SRA-Singapore-Raffles_Lighthouse 
##                        "#000000"                        "#9d529e" 
##        SRA-Singapore-Kusu_Island    SRA-Singapore-St_Johns_Island 
##                        "#bc80bd"                        "#e1cbe1"
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

```
## Warning: Duplicate levels detected

## Warning: Duplicate levels detected
```

```{=html}
<div id="htmlwidget-76dc36a5017027c0a93b" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-76dc36a5017027c0a93b">{"x":{"visdat":{"14f536e524b0":["function () ","plotlyVisDat"]},"cur_data":"14f536e524b0","attrs":{"14f536e524b0":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305","SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729","SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301","SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#9d529e","#bc80bd","#e1cbe1","#1b9e77","#d95f02","#01665e","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000","#9d529e","#bc80bd","#e1cbe1"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Group","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.087909776151017,-0.0914000446598861,-0.0924680310078454,-0.0907406066449548,-0.0883683173656369,-0.0913173498012151,-0.0878304958163144],"y":[0.0398843854978978,0.043310300762778,0.0433853613585064,0.0435593051756594,0.0404960693662204,0.043320610938873,0.0402791294894221],"text":["Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP8_1051","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP8_1459","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP7_2072"],"mode":"markers","marker":{"color":"rgba(215,59,45,1)","size":11,"line":{"color":"rgba(215,59,45,1)"}},"type":"scatter","name":"Group1","textfont":{"color":"rgba(215,59,45,1)"},"error_y":{"color":"rgba(215,59,45,1)"},"error_x":{"color":"rgba(215,59,45,1)"},"line":{"color":"rgba(215,59,45,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.102931740088072,-0.106187237350478,-0.106185203085442,-0.104187840337855,-0.107842222051049,-0.105137952779799,-0.104607038290398,-0.10379450703178,-0.105986562733763,-0.107308438177123,-0.103301573525756,-0.105934440468215,-0.104767734697234,-0.105978528275292,-0.104303934550798,-0.103756527848323,-0.105841199096512,-0.10705706568918,-0.1066772864339,-0.108151181773894,-0.105074285160288,-0.106633546737049,-0.103895187267956,-0.10976540915027,-0.105562793194793,-0.105858609501087,-0.107352039518693,-0.106201894712033,-0.105402597385353,-0.107828141187385,-0.105810413876808,-0.104957010950411,-0.106864503651448,-0.104831356670049],"y":[0.0583114037291299,0.0591338139086443,0.059439324134166,0.0588554955194258,0.0593865399930813,0.0590772302288572,0.0585410685859713,0.0583602260739379,0.0591559574369448,0.0600078226059689,0.0577667218035774,0.0589743542811655,0.0591237754043505,0.059621251500395,0.0588136783989166,0.0582347754157852,0.0593937192234727,0.0592932713056906,0.0595319493354748,0.0598297648308543,0.0591559025646922,0.0598405935082927,0.0587420992920614,0.0599992588392903,0.0590714082537799,0.0595093510358138,0.0593970312861219,0.0594788171403164,0.0593624284293784,0.0597737087864577,0.0593140578879109,0.0589096180027599,0.0596102260854338,0.059057455117443],"text":["Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"color":"rgba(174,91,48,1)","size":11,"line":{"color":"rgba(174,91,48,1)"}},"type":"scatter","name":"Group2","textfont":{"color":"rgba(174,91,48,1)"},"error_y":{"color":"rgba(174,91,48,1)"},"error_x":{"color":"rgba(174,91,48,1)"},"line":{"color":"rgba(174,91,48,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0343828979334148,-0.0334183458632226,-0.0340005388918923,-0.0336200180833008,-0.0342253906080032,-0.0357189274312273,-0.0358114089398999,-0.034171903543001,-0.0339310174923889,-0.0348147968354916,-0.035857124789644,-0.0341932458154212,-0.0356356608675336,-0.0354267739391556,-0.0341794575528002,-0.0338762356712997,-0.0358602044119605,-0.0332761634267843,-0.0343280688808464,-0.0348294776145679,-0.0343789576573794,-0.0355174542961062,-0.0345033148000463,-0.0339776886546019],"y":[-0.0698024674140147,-0.0682093267782368,-0.0693218584637933,-0.068977123678469,-0.0695490957872131,-0.0713667951076149,-0.0713432767171898,-0.0695406805868687,-0.0693928403510257,-0.0701514182038005,-0.0714264330025069,-0.0695427972946983,-0.071358184958669,-0.0709295353292799,-0.0694587977815881,-0.0692892975633885,-0.0712093331583534,-0.0687240085303204,-0.0699986600914017,-0.0694356439685457,-0.069597699016445,-0.0712288947357569,-0.069920442680824,-0.0695084809585879],"text":["Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP9_1594","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP8_2564","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP8_1765","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1820"],"mode":"markers","marker":{"color":"rgba(99,117,103,1)","size":11,"line":{"color":"rgba(99,117,103,1)"}},"type":"scatter","name":"Group3","textfont":{"color":"rgba(99,117,103,1)"},"error_y":{"color":"rgba(99,117,103,1)"},"error_x":{"color":"rgba(99,117,103,1)"},"line":{"color":"rgba(99,117,103,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.00477958827065582,-0.00479260315074789,-0.00481782622456821,-0.00505227022641038,-0.00477717203855318,-0.00452204591203607,-0.00481800683413358,-0.00483883645441126],"y":[-0.0771068929326349,-0.0777587572117368,-0.0772890939851324,-0.0797063965305195,-0.080354501550278,-0.0779676396996429,-0.0778657014473013,-0.0764748327898629],"text":["Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP6_1542","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP3_2750","Pacuta_HTAC_TP4_1581","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP8_1184"],"mode":"markers","marker":{"color":"rgba(61,130,147,1)","size":11,"line":{"color":"rgba(61,130,147,1)"}},"type":"scatter","name":"Group4","textfont":{"color":"rgba(61,130,147,1)"},"error_y":{"color":"rgba(61,130,147,1)"},"error_x":{"color":"rgba(61,130,147,1)"},"line":{"color":"rgba(61,130,147,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.00169484858710476,-0.00178748246462859,-0.0019816952312961,-0.00179917755411026,-0.001691218333494,-0.00179039969233185,-0.00181895733519894,-0.00141340548294186,-0.00170913522588395],"y":[-0.0485168343000597,-0.0489522400541657,-0.0494550825577194,-0.0487782666228017,-0.0485871090329391,-0.0495546483599106,-0.0485092018271027,-0.048937477818861,-0.0485269045078775],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP9_1141","Pacuta_ATHC_TP5_2212","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486"],"mode":"markers","marker":{"color":"rgba(171,44,134,1)","size":11,"line":{"color":"rgba(171,44,134,1)"}},"type":"scatter","name":"Group5","textfont":{"color":"rgba(171,44,134,1)"},"error_y":{"color":"rgba(171,44,134,1)"},"error_x":{"color":"rgba(171,44,134,1)"},"line":{"color":"rgba(171,44,134,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0638817942846192,0.0641298473888903,0.0638866881129803,0.0648487908173522,0.0636631574644534,0.0643254744612316,0.0641245609584253,0.063736129346863,0.0633920165001564,0.0634745699282088,0.0635648495327394,0.0643894930195544,0.0639150508997213,0.0635435517184501,0.0635360280834863,0.0637396993125961,0.0637711767947223,0.0645599077495936,0.0630968914742519,0.0638517977136142,0.0648159861792405,0.0642917662615296,0.0636154171348093,0.0637934497063493,0.0639506052585641,0.0646375520744721,0.063987644124737],"y":[-0.108586861886645,-0.109351572775111,-0.108837491294561,-0.110741562310039,-0.106522759915366,-0.109647972867507,-0.109703766746795,-0.109163191641733,-0.108363007808047,-0.108680992985497,-0.108527029031323,-0.109677068987259,-0.109807614829258,-0.108810065357337,-0.107973157083784,-0.109012278329329,-0.109061629793114,-0.110056214206267,-0.107649833286902,-0.108783354618302,-0.110733300895891,-0.109227819438082,-0.108898854501322,-0.108629931381199,-0.108868770720931,-0.110512882824837,-0.109397382064606],"text":["Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP7_1047","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP8_1329","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP5_1168"],"mode":"markers","marker":{"color":"rgba(89,50,89,1)","size":11,"line":{"color":"rgba(89,50,89,1)"}},"type":"scatter","name":"Group6","textfont":{"color":"rgba(89,50,89,1)"},"error_y":{"color":"rgba(89,50,89,1)"},"error_x":{"color":"rgba(89,50,89,1)"},"line":{"color":"rgba(89,50,89,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0393778589823457,0.0408739497214032],"y":[-0.0299824503655067,-0.0435376700074661],"text":["Pacuta_ATAC_TP5_1059","Pacuta_HTHC_TP1_2210"],"mode":"markers","marker":{"color":"rgba(80,80,80,1)","size":11,"line":{"color":"rgba(80,80,80,1)"}},"type":"scatter","name":"Group7","textfont":{"color":"rgba(80,80,80,1)"},"error_y":{"color":"rgba(80,80,80,1)"},"error_x":{"color":"rgba(80,80,80,1)"},"line":{"color":"rgba(80,80,80,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0552945321305707,0.0552982105742001],"y":[-0.0130829460627016,-0.0107659213878848],"text":["Pacuta_HTAC_TP1_1653","Pacuta_HTHC_TP6_1721"],"mode":"markers","marker":{"color":"rgba(225,203,225,1)","size":11,"line":{"color":"rgba(225,203,225,1)"}},"type":"scatter","name":"Group8","textfont":{"color":"rgba(225,203,225,1)"},"error_y":{"color":"rgba(225,203,225,1)"},"error_x":{"color":"rgba(225,203,225,1)"},"line":{"color":"rgba(225,203,225,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.11844319003458,0.118564624451305,0.118465504236608,0.119086502541586,0.119392393954628,0.118709326852248,0.11964958995854,0.119795092420036,0.119062613829404,0.118320289481152,0.119018651345233,0.118980639082607],"y":[0.114729831772802,0.115025388351994,0.114758498420943,0.115505960273821,0.115610091480121,0.115022070439929,0.115988277397334,0.116039673334351,0.11530425775305,0.11470305809261,0.115318187072981,0.115293792425329],"text":["SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301"],"mode":"markers","marker":{"color":"rgba(147,183,174,1)","size":11,"line":{"color":"rgba(147,183,174,1)"}},"type":"scatter","name":"SRA-Singapore-Kusu_Island","textfont":{"color":"rgba(147,183,174,1)"},"error_y":{"color":"rgba(147,183,174,1)"},"error_x":{"color":"rgba(147,183,174,1)"},"line":{"color":"rgba(147,183,174,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.100985628893805,0.101620413665608,0.102258720215135,0.10163719935419,0.101738888368613,0.101343137892189,0.101607040727096,0.102021945008195],"y":[0.0808842099049849,0.0816525623361383,0.0823754978204067,0.0815820601417709,0.0816996210813305,0.0811079900231838,0.0814535623365954,0.0820433341337974],"text":["SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729"],"mode":"markers","marker":{"color":"rgba(198,148,198,1)","size":11,"line":{"color":"rgba(198,148,198,1)"}},"type":"scatter","name":"SRA-Singapore-Raffles_Lighthouse","textfont":{"color":"rgba(198,148,198,1)"},"error_y":{"color":"rgba(198,148,198,1)"},"error_x":{"color":"rgba(198,148,198,1)"},"line":{"color":"rgba(198,148,198,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.111737971811951,0.109657420020779,0.110548962437426,0.110006695254346,0.111248954267185,0.111330167987098,0.110822988020866,0.111442142286358,0.110877149701611,0.111476930431115,0.110581132168019,0.11072502771064],"y":[0.108507746318127,0.106562511580627,0.10745108442712,0.106925354595738,0.108134494655654,0.108138913926835,0.107692974243287,0.108257976367335,0.107690038719373,0.108118853532052,0.107357989306043,0.107472123388406],"text":["SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"color":"rgba(181,122,50,1)","size":11,"line":{"color":"rgba(181,122,50,1)"}},"type":"scatter","name":"SRA-Singapore-St_Johns_Island","textfont":{"color":"rgba(181,122,50,1)"},"error_y":{"color":"rgba(181,122,50,1)"},"error_x":{"color":"rgba(181,122,50,1)"},"line":{"color":"rgba(181,122,50,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.00860652917958696,0.031073303286445,-0.0289521592857298,-0.0776338915679618,-0.0353248816206353,-0.0711002952011141],"y":[-0.0315839936553536,-0.0439052883914302,0.0086939526525521,0.0357197992029017,-0.0125408584784229,0.00601249859220779],"text":["Pacuta_ATAC_TP7_1445","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP6_1254","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP5_1415"],"mode":"markers","marker":{"color":"rgba(182,118,182,1)","size":11,"line":{"color":"rgba(182,118,182,1)"}},"type":"scatter","name":"Ungroup","textfont":{"color":"rgba(182,118,182,1)"},"error_y":{"color":"rgba(182,118,182,1)"},"error_x":{"color":"rgba(182,118,182,1)"},"line":{"color":"rgba(182,118,182,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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

```
## Warning: Duplicate levels detected

## Warning: Duplicate levels detected
```

```{=html}
<div id="htmlwidget-32bc199eb36295564175" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-32bc199eb36295564175">{"x":{"visdat":{"14f57dcbbaf7":["function () ","plotlyVisDat"]},"cur_data":"14f57dcbbaf7","attrs":{"14f57dcbbaf7":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305","SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729","SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301","SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#9d529e","#bc80bd","#e1cbe1","#1b9e77","#d95f02","#01665e","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000","#9d529e","#bc80bd","#e1cbe1"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Reef","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0334183458632226,-0.087909776151017,-0.105986562733763,0.0643254744612316,-0.105934440468215,-0.104767734697234,-0.105978528275292,0.0639150508997213,-0.103756527848323,0.0635360280834863,-0.105841199096512,-0.0354267739391556,0.0552945321305707,-0.0883683173656369,-0.0878304958163144,-0.0353248816206353,-0.107352039518693,-0.106201894712033,-0.105402597385353,-0.104957010950411,-0.106864503651448],"y":[-0.0682093267782368,0.0398843854978978,0.0591559574369448,-0.109647972867507,0.0589743542811655,0.0591237754043505,0.059621251500395,-0.109807614829258,0.0582347754157852,-0.107973157083784,0.0593937192234727,-0.0709295353292799,-0.0130829460627016,0.0404960693662204,0.0402791294894221,-0.0125408584784229,0.0593970312861219,0.0594788171403164,0.0593624284293784,0.0589096180027599,0.0596102260854338],"text":["Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP11_1777","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_2873","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP7_2072","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202"],"mode":"markers","marker":{"color":"rgba(205,142,106,1)","size":11,"line":{"color":"rgba(205,142,106,1)"}},"type":"scatter","name":"HIMB","textfont":{"color":"rgba(205,142,106,1)"},"error_y":{"color":"rgba(205,142,106,1)"},"error_x":{"color":"rgba(205,142,106,1)"},"line":{"color":"rgba(205,142,106,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0343828979334148,-0.00477958827065582,-0.0358114089398999,-0.00479260315074789,0.0641298473888903,0.0636631574644534,-0.00481782622456821,-0.103301573525756,0.0635648495327394,-0.0348147968354916,-0.0356356608675336,-0.00477717203855318,-0.105074285160288,0.0638517977136142,0.0648159861792405,-0.10976540915027,-0.0348294776145679,0.0637934497063493,-0.00452204591203607,0.0552982105742001,-0.0345033148000463,-0.00481800683413358,-0.00483883645441126,-0.107828141187385],"y":[-0.0698024674140147,-0.0771068929326349,-0.0713432767171898,-0.0777587572117368,-0.109351572775111,-0.106522759915366,-0.0772890939851324,0.0577667218035774,-0.108527029031323,-0.0701514182038005,-0.071358184958669,-0.080354501550278,0.0591559025646922,-0.108783354618302,-0.110733300895891,0.0599992588392903,-0.0694356439685457,-0.108629931381199,-0.0779676396996429,-0.0107659213878848,-0.069920442680824,-0.0778657014473013,-0.0764748327898629,0.0597737087864577],"text":["Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_2409","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709"],"mode":"markers","marker":{"color":"rgba(228,147,53,1)","size":11,"line":{"color":"rgba(228,147,53,1)"}},"type":"scatter","name":"Lilipuna.Fringe","textfont":{"color":"rgba(228,147,53,1)"},"error_y":{"color":"rgba(228,147,53,1)"},"error_x":{"color":"rgba(228,147,53,1)"},"line":{"color":"rgba(228,147,53,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.106185203085442,-0.104187840337855,0.0393778589823457,-0.00860652917958696,-0.0914000446598861,0.0648487908173522,0.0641245609584253,0.0633920165001564,0.0634745699282088,-0.108151181773894,-0.00170913522588395,0.0642917662615296,-0.0711002952011141,-0.0339776886546019],"y":[0.059439324134166,0.0588554955194258,-0.0299824503655067,-0.0315839936553536,0.043310300762778,-0.110741562310039,-0.109703766746795,-0.108363007808047,-0.108680992985497,0.0598297648308543,-0.0485269045078775,-0.109227819438082,0.00601249859220779,-0.0695084809585879],"text":["Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP8_1051","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP9_1486","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP7_1820"],"mode":"markers","marker":{"color":"rgba(170,154,194,1)","size":11,"line":{"color":"rgba(170,154,194,1)"}},"type":"scatter","name":"Reef.11.13","textfont":{"color":"rgba(170,154,194,1)"},"error_y":{"color":"rgba(170,154,194,1)"},"error_x":{"color":"rgba(170,154,194,1)"},"line":{"color":"rgba(170,154,194,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.106187237350478,-0.0340005388918923,-0.0342253906080032,0.0638817942846192,-0.105137952779799,-0.10379450703178,-0.0289521592857298,-0.107308438177123,-0.0907406066449548,-0.035857124789644,-0.104303934550798,-0.0341794575528002,-0.1066772864339,-0.0358602044119605,-0.106633546737049,-0.103895187267956,-0.105562793194793,0.0639506052585641,-0.105810413876808,-0.104831356670049],"y":[0.0591338139086443,-0.0693218584637933,-0.0695490957872131,-0.108586861886645,0.0590772302288572,0.0583602260739379,0.0086939526525521,0.0600078226059689,0.0435593051756594,-0.0714264330025069,0.0588136783989166,-0.0694587977815881,0.0595319493354748,-0.0712093331583534,0.0598405935082927,0.0587420992920614,0.0590714082537799,-0.108868770720931,0.0593140578879109,0.059057455117443],"text":["Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP9_1451","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP8_2513","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"color":"rgba(241,222,202,1)","size":11,"line":{"color":"rgba(241,222,202,1)"}},"type":"scatter","name":"Reef.18","textfont":{"color":"rgba(241,222,202,1)"},"error_y":{"color":"rgba(241,222,202,1)"},"error_x":{"color":"rgba(241,222,202,1)"},"line":{"color":"rgba(241,222,202,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0336200180833008,-0.0019816952312961,-0.0339310174923889,0.031073303286445,0.0638866881129803,0.063736129346863,-0.00179039969233185,0.0643894930195544,0.0635435517184501,-0.0341932458154212,0.0637711767947223,0.0645599077495936,0.0630968914742519,-0.00181895733519894,-0.00141340548294186,0.0646375520744721,0.063987644124737],"y":[-0.068977123678469,-0.0494550825577194,-0.0693928403510257,-0.0439052883914302,-0.108837491294561,-0.109163191641733,-0.0495546483599106,-0.109677068987259,-0.108810065357337,-0.0695427972946983,-0.109061629793114,-0.110056214206267,-0.107649833286902,-0.0485092018271027,-0.048937477818861,-0.110512882824837,-0.109397382064606],"text":["Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP9_1594","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP9_1302","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP5_1168"],"mode":"markers","marker":{"color":"rgba(159,84,160,1)","size":11,"line":{"color":"rgba(159,84,160,1)"}},"type":"scatter","name":"Reef.35.36","textfont":{"color":"rgba(159,84,160,1)"},"error_y":{"color":"rgba(159,84,160,1)"},"error_x":{"color":"rgba(159,84,160,1)"},"line":{"color":"rgba(159,84,160,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.00169484858710476,-0.102931740088072,-0.00178748246462859,-0.107842222051049,-0.0357189274312273,-0.034171903543001,-0.00179917755411026,-0.104607038290398,-0.001691218333494,-0.00505227022641038,-0.0924680310078454,-0.0776338915679618,-0.10705706568918,-0.0338762356712997,0.0637396993125961,-0.0913173498012151,-0.0332761634267843,-0.0343280688808464,-0.0343789576573794,0.0636154171348093,0.0408739497214032,-0.105858609501087,-0.0355174542961062],"y":[-0.0485168343000597,0.0583114037291299,-0.0489522400541657,0.0593865399930813,-0.0713667951076149,-0.0695406805868687,-0.0487782666228017,0.0585410685859713,-0.0485871090329391,-0.0797063965305195,0.0433853613585064,0.0357197992029017,0.0592932713056906,-0.0692892975633885,-0.109012278329329,0.043320610938873,-0.0687240085303204,-0.0699986600914017,-0.069597699016445,-0.108898854501322,-0.0435376700074661,0.0595093510358138,-0.0712288947357569],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP6_1254","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP8_1765","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP4_2195"],"mode":"markers","marker":{"color":"rgba(84,39,136,1)","size":11,"line":{"color":"rgba(84,39,136,1)"}},"type":"scatter","name":"Reef.42.43","textfont":{"color":"rgba(84,39,136,1)"},"error_y":{"color":"rgba(84,39,136,1)"},"error_x":{"color":"rgba(84,39,136,1)"},"line":{"color":"rgba(84,39,136,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.11844319003458,0.118564624451305,0.118465504236608,0.119086502541586,0.119392393954628,0.118709326852248,0.11964958995854,0.119795092420036,0.119062613829404,0.118320289481152,0.119018651345233,0.118980639082607],"y":[0.114729831772802,0.115025388351994,0.114758498420943,0.115505960273821,0.115610091480121,0.115022070439929,0.115988277397334,0.116039673334351,0.11530425775305,0.11470305809261,0.115318187072981,0.115293792425329],"text":["SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301"],"mode":"markers","marker":{"color":"rgba(147,183,174,1)","size":11,"line":{"color":"rgba(147,183,174,1)"}},"type":"scatter","name":"SRA-Singapore-Kusu_Island","textfont":{"color":"rgba(147,183,174,1)"},"error_y":{"color":"rgba(147,183,174,1)"},"error_x":{"color":"rgba(147,183,174,1)"},"line":{"color":"rgba(147,183,174,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.100985628893805,0.101620413665608,0.102258720215135,0.10163719935419,0.101738888368613,0.101343137892189,0.101607040727096,0.102021945008195],"y":[0.0808842099049849,0.0816525623361383,0.0823754978204067,0.0815820601417709,0.0816996210813305,0.0811079900231838,0.0814535623365954,0.0820433341337974],"text":["SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729"],"mode":"markers","marker":{"color":"rgba(198,148,198,1)","size":11,"line":{"color":"rgba(198,148,198,1)"}},"type":"scatter","name":"SRA-Singapore-Raffles_Lighthouse","textfont":{"color":"rgba(198,148,198,1)"},"error_y":{"color":"rgba(198,148,198,1)"},"error_x":{"color":"rgba(198,148,198,1)"},"line":{"color":"rgba(198,148,198,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.111737971811951,0.109657420020779,0.110548962437426,0.110006695254346,0.111248954267185,0.111330167987098,0.110822988020866,0.111442142286358,0.110877149701611,0.111476930431115,0.110581132168019,0.11072502771064],"y":[0.108507746318127,0.106562511580627,0.10745108442712,0.106925354595738,0.108134494655654,0.108138913926835,0.107692974243287,0.108257976367335,0.107690038719373,0.108118853532052,0.107357989306043,0.107472123388406],"text":["SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"color":"rgba(181,122,50,1)","size":11,"line":{"color":"rgba(181,122,50,1)"}},"type":"scatter","name":"SRA-Singapore-St_Johns_Island","textfont":{"color":"rgba(181,122,50,1)"},"error_y":{"color":"rgba(181,122,50,1)"},"error_x":{"color":"rgba(181,122,50,1)"},"line":{"color":"rgba(181,122,50,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-56c1a38cbdc670bb8240" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-56c1a38cbdc670bb8240">{"x":{"visdat":{"14f54a5feea6":["function () ","plotlyVisDat"]},"cur_data":"14f54a5feea6","attrs":{"14f54a5feea6":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305","SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729","SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301","SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#9d529e","#bc80bd","#e1cbe1","#1b9e77","#d95f02","#01665e","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000","#9d529e","#bc80bd","#e1cbe1"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Ploidy","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.00169484858710476,-0.0343828979334148,-0.0334183458632226,-0.00477958827065582,-0.102931740088072,-0.087909776151017,-0.106187237350478,-0.00178748246462859,-0.0340005388918923,-0.106185203085442,-0.0336200180833008,-0.104187840337855,-0.0342253906080032,-0.107842222051049,-0.0357189274312273,0.0393778589823457,-0.0358114089398999,-0.034171903543001,0.0638817942846192,-0.0019816952312961,-0.00479260315074789,0.0641298473888903,-0.00860652917958696,-0.105137952779799,-0.0914000446598861,-0.00179917755411026,-0.104607038290398,-0.001691218333494,-0.0339310174923889,-0.10379450703178,0.031073303286445,0.0638866881129803,-0.0289521592857298,-0.105986562733763,0.0648487908173522,0.0636631574644534,0.0643254744612316,-0.00481782622456821,0.0641245609584253,0.063736129346863,-0.107308438177123,-0.00505227022641038,-0.103301573525756,-0.0924680310078454,0.0633920165001564,0.0634745699282088,-0.00179039969233185,-0.105934440468215,-0.0776338915679618,-0.104767734697234,0.0635648495327394,0.0643894930195544,-0.0348147968354916,-0.105978528275292,-0.0907406066449548,-0.035857124789644,0.0639150508997213,-0.104303934550798,-0.103756527848323,0.0635435517184501,0.0635360280834863,-0.0341932458154212,-0.105841199096512,-0.0356356608675336,-0.10705706568918,-0.0354267739391556,0.0552945321305707,-0.0341794575528002,-0.1066772864339,-0.0338762356712997,-0.0883683173656369,0.0637396993125961,-0.00477717203855318,-0.0358602044119605,-0.0913173498012151,0.0637711767947223,0.0645599077495936,-0.108151181773894,0.0630968914742519,-0.105074285160288,-0.00181895733519894,0.0638517977136142,-0.106633546737049,-0.0878304958163144,0.0648159861792405,-0.0332761634267843,-0.103895187267956,-0.00141340548294186,-0.00170913522588395,-0.10976540915027,-0.0343280688808464,-0.0348294776145679,-0.105562793194793,0.0642917662615296,-0.0353248816206353,-0.0343789576573794,0.0636154171348093,0.0408739497214032,-0.105858609501087,0.0637934497063493,0.0639506052585641,0.0646375520744721,-0.00452204591203607,-0.0355174542961062,0.063987644124737,-0.0711002952011141,-0.107352039518693,-0.106201894712033,-0.105402597385353,0.0552982105742001,-0.0345033148000463,-0.00481800683413358,-0.0339776886546019,-0.00483883645441126,-0.107828141187385,-0.105810413876808,-0.104957010950411,-0.106864503651448,-0.104831356670049,0.100985628893805,0.101620413665608,0.102258720215135,0.10163719935419,0.101738888368613,0.101343137892189,0.101607040727096,0.102021945008195,0.11844319003458,0.118564624451305,0.118465504236608,0.119086502541586,0.119392393954628,0.118709326852248,0.11964958995854,0.119795092420036,0.119062613829404,0.118320289481152,0.119018651345233,0.118980639082607,0.111737971811951,0.109657420020779,0.110548962437426,0.110006695254346,0.111248954267185,0.111330167987098,0.110822988020866,0.111442142286358,0.110877149701611,0.111476930431115,0.110581132168019,0.11072502771064],"y":[-0.0485168343000597,-0.0698024674140147,-0.0682093267782368,-0.0771068929326349,0.0583114037291299,0.0398843854978978,0.0591338139086443,-0.0489522400541657,-0.0693218584637933,0.059439324134166,-0.068977123678469,0.0588554955194258,-0.0695490957872131,0.0593865399930813,-0.0713667951076149,-0.0299824503655067,-0.0713432767171898,-0.0695406805868687,-0.108586861886645,-0.0494550825577194,-0.0777587572117368,-0.109351572775111,-0.0315839936553536,0.0590772302288572,0.043310300762778,-0.0487782666228017,0.0585410685859713,-0.0485871090329391,-0.0693928403510257,0.0583602260739379,-0.0439052883914302,-0.108837491294561,0.0086939526525521,0.0591559574369448,-0.110741562310039,-0.106522759915366,-0.109647972867507,-0.0772890939851324,-0.109703766746795,-0.109163191641733,0.0600078226059689,-0.0797063965305195,0.0577667218035774,0.0433853613585064,-0.108363007808047,-0.108680992985497,-0.0495546483599106,0.0589743542811655,0.0357197992029017,0.0591237754043505,-0.108527029031323,-0.109677068987259,-0.0701514182038005,0.059621251500395,0.0435593051756594,-0.0714264330025069,-0.109807614829258,0.0588136783989166,0.0582347754157852,-0.108810065357337,-0.107973157083784,-0.0695427972946983,0.0593937192234727,-0.071358184958669,0.0592932713056906,-0.0709295353292799,-0.0130829460627016,-0.0694587977815881,0.0595319493354748,-0.0692892975633885,0.0404960693662204,-0.109012278329329,-0.080354501550278,-0.0712093331583534,0.043320610938873,-0.109061629793114,-0.110056214206267,0.0598297648308543,-0.107649833286902,0.0591559025646922,-0.0485092018271027,-0.108783354618302,0.0598405935082927,0.0402791294894221,-0.110733300895891,-0.0687240085303204,0.0587420992920614,-0.048937477818861,-0.0485269045078775,0.0599992588392903,-0.0699986600914017,-0.0694356439685457,0.0590714082537799,-0.109227819438082,-0.0125408584784229,-0.069597699016445,-0.108898854501322,-0.0435376700074661,0.0595093510358138,-0.108629931381199,-0.108868770720931,-0.110512882824837,-0.0779676396996429,-0.0712288947357569,-0.109397382064606,0.00601249859220779,0.0593970312861219,0.0594788171403164,0.0593624284293784,-0.0107659213878848,-0.069920442680824,-0.0778657014473013,-0.0695084809585879,-0.0764748327898629,0.0597737087864577,0.0593140578879109,0.0589096180027599,0.0596102260854338,0.059057455117443,0.0808842099049849,0.0816525623361383,0.0823754978204067,0.0815820601417709,0.0816996210813305,0.0811079900231838,0.0814535623365954,0.0820433341337974,0.114729831772802,0.115025388351994,0.114758498420943,0.115505960273821,0.115610091480121,0.115022070439929,0.115988277397334,0.116039673334351,0.11530425775305,0.11470305809261,0.115318187072981,0.115293792425329,0.108507746318127,0.106562511580627,0.10745108442712,0.106925354595738,0.108134494655654,0.108138913926835,0.107692974243287,0.108257976367335,0.107690038719373,0.108118853532052,0.107357989306043,0.107472123388406],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305","SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729","SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301","SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"colorbar":{"title":"ploidy","ticklen":2},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(235,157,69,1)"],["0.0833333333333335","rgba(184,94,12,1)"],["0.125","rgba(217,201,206,1)"],["0.166666666666667","rgba(249,223,191,1)"],["0.208333333333333","rgba(176,158,193,1)"],["0.25","rgba(157,98,167,1)"],["0.291666666666667","rgba(179,115,180,1)"],["0.333333333333333","rgba(213,178,213,1)"],["0.375","rgba(120,177,158,1)"],["0.416666666666667","rgba(168,128,61,1)"],["0.458333333333333","rgba(132,104,62,1)"],["0.5","rgba(157,119,64,1)"],["0.541666666666667","rgba(231,88,78,1)"],["0.583333333333333","rgba(163,48,137,1)"],["0.625","rgba(104,102,119,1)"],["0.666666666666667","rgba(67,146,96,1)"],["0.708333333333333","rgba(113,112,139,1)"],["0.75","rgba(190,80,37,1)"],["0.791666666666667","rgba(212,66,50,1)"],["0.833333333333333","rgba(106,106,106,1)"],["0.875","rgba(24,16,24,1)"],["0.916666666666667","rgba(160,86,161,1)"],["0.958333333333333","rgba(190,131,190,1)"],["1","rgba(225,203,225,1)"]],"showscale":false,"color":[2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,2,3,3,2,2,3,2,2,3,3,2,3,2,3,3,2,2,2,3,2,2,2,3,2,2,3,3,3,3,2,2,2,3,3,3,2,2,3,3,3,3,2,3,3,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,3,2,2,3,2,3,2,2,3,3,2,3,3,2,2,3,3,3,3,2,3,3,2,2,3,2,2,2,3,3,2,2,3,3,3,2,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"size":11,"line":{"colorbar":{"title":"","ticklen":2},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(235,157,69,1)"],["0.0833333333333335","rgba(184,94,12,1)"],["0.125","rgba(217,201,206,1)"],["0.166666666666667","rgba(249,223,191,1)"],["0.208333333333333","rgba(176,158,193,1)"],["0.25","rgba(157,98,167,1)"],["0.291666666666667","rgba(179,115,180,1)"],["0.333333333333333","rgba(213,178,213,1)"],["0.375","rgba(120,177,158,1)"],["0.416666666666667","rgba(168,128,61,1)"],["0.458333333333333","rgba(132,104,62,1)"],["0.5","rgba(157,119,64,1)"],["0.541666666666667","rgba(231,88,78,1)"],["0.583333333333333","rgba(163,48,137,1)"],["0.625","rgba(104,102,119,1)"],["0.666666666666667","rgba(67,146,96,1)"],["0.708333333333333","rgba(113,112,139,1)"],["0.75","rgba(190,80,37,1)"],["0.791666666666667","rgba(212,66,50,1)"],["0.833333333333333","rgba(106,106,106,1)"],["0.875","rgba(24,16,24,1)"],["0.916666666666667","rgba(160,86,161,1)"],["0.958333333333333","rgba(190,131,190,1)"],["1","rgba(225,203,225,1)"]],"showscale":false,"color":[2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,2,3,3,2,2,3,2,2,3,3,2,3,2,3,3,2,2,2,3,2,2,2,3,2,2,3,3,3,3,2,2,2,3,3,3,2,2,3,3,3,3,2,3,3,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,3,2,2,3,2,3,2,2,3,3,2,3,3,2,2,3,3,3,3,2,3,3,2,2,3,2,2,2,3,3,2,2,3,3,3,2,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}},"type":"scatter","xaxis":"x","yaxis":"y","frame":null},{"x":[-0.10976540915027,0.119795092420036],"y":[-0.110741562310039,0.116039673334351],"type":"scatter","mode":"markers","opacity":0,"hoverinfo":"none","showlegend":false,"marker":{"colorbar":{"title":"ploidy","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(235,157,69,1)"],["0.0833333333333335","rgba(184,94,12,1)"],["0.125","rgba(217,201,206,1)"],["0.166666666666667","rgba(249,223,191,1)"],["0.208333333333333","rgba(176,158,193,1)"],["0.25","rgba(157,98,167,1)"],["0.291666666666667","rgba(179,115,180,1)"],["0.333333333333333","rgba(213,178,213,1)"],["0.375","rgba(120,177,158,1)"],["0.416666666666667","rgba(168,128,61,1)"],["0.458333333333333","rgba(132,104,62,1)"],["0.5","rgba(157,119,64,1)"],["0.541666666666667","rgba(231,88,78,1)"],["0.583333333333333","rgba(163,48,137,1)"],["0.625","rgba(104,102,119,1)"],["0.666666666666667","rgba(67,146,96,1)"],["0.708333333333333","rgba(113,112,139,1)"],["0.75","rgba(190,80,37,1)"],["0.791666666666667","rgba(212,66,50,1)"],["0.833333333333333","rgba(106,106,106,1)"],["0.875","rgba(24,16,24,1)"],["0.916666666666667","rgba(160,86,161,1)"],["0.958333333333333","rgba(190,131,190,1)"],["1","rgba(225,203,225,1)"]],"showscale":true,"color":[2,3],"line":{"color":"rgba(255,127,14,1)"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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

```
## Warning: Duplicate levels detected

## Warning: Duplicate levels detected
```

```{=html}
<div id="htmlwidget-9bede1a7a1f58d47c48b" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-9bede1a7a1f58d47c48b">{"x":{"visdat":{"14f560b24a78":["function () ","plotlyVisDat"]},"cur_data":"14f560b24a78","attrs":{"14f560b24a78":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305","SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729","SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301","SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#9d529e","#bc80bd","#e1cbe1","#1b9e77","#d95f02","#01665e","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000","#9d529e","#bc80bd","#e1cbe1"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Group","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0538333961763806,-0.0515835767663872,-0.0611219706757994,-0.0494836920981077,-0.0492011872391827,-0.0530080003490253,-0.0486896121303684],"y":[-0.0873118374741984,-0.086337688274275,-0.102383402709373,-0.083815992999526,-0.0802508105936718,-0.0891606333044414,-0.0797514299002811],"text":["Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP8_1051","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP8_1459","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP7_2072"],"mode":"markers","marker":{"color":"rgba(215,59,45,1)","size":11,"line":{"color":"rgba(215,59,45,1)"}},"type":"scatter","name":"Group1","textfont":{"color":"rgba(215,59,45,1)"},"error_y":{"color":"rgba(215,59,45,1)"},"error_x":{"color":"rgba(215,59,45,1)"},"line":{"color":"rgba(215,59,45,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0489450962004277,-0.059013649486611,-0.0573236031218694,-0.0509358823936688,-0.0632956521039723,-0.0553564465445297,-0.0507270406143958,-0.0505389556753465,-0.0609230079724621,-0.0622248029799049,-0.0540851552531239,-0.0551353549725535,-0.0526178381174508,-0.0577909828311013,-0.0514635109927413,-0.0498387256099919,-0.0557205653119484,-0.0625987047830899,-0.060109103379336,-0.0636132014351463,-0.0531368267822988,-0.0608057931386368,-0.0506396577851495,-0.0693887345846079,-0.0547479723576256,-0.0555465740849365,-0.0649283315075732,-0.0577912527954553,-0.0552832363622027,-0.0643050275114951,-0.0557120353660638,-0.0521838737667409,-0.0589156035430733,-0.0524243967376363],"y":[-0.0969941362589877,-0.115412262806254,-0.113241364175914,-0.100660992690893,-0.123706890147924,-0.109261012068346,-0.0991740693413415,-0.0992654111619984,-0.119875579377021,-0.123132070705555,-0.105571273752684,-0.108392763386094,-0.104354620550549,-0.114482523820953,-0.101529294398627,-0.097747252794264,-0.1094417660821,-0.122431189919647,-0.119040080078323,-0.125009952400077,-0.105032403033994,-0.120336263801797,-0.0994384469639176,-0.135196923550645,-0.107091606485858,-0.110231462529279,-0.126933139530872,-0.114163388756836,-0.109494120805786,-0.12594205247791,-0.109774013430729,-0.102511337341296,-0.115866359564381,-0.102708465114943],"text":["Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"color":"rgba(174,91,48,1)","size":11,"line":{"color":"rgba(174,91,48,1)"}},"type":"scatter","name":"Group2","textfont":{"color":"rgba(174,91,48,1)"},"error_y":{"color":"rgba(174,91,48,1)"},"error_x":{"color":"rgba(174,91,48,1)"},"line":{"color":"rgba(174,91,48,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.055223595654752,-0.0509810994843474,-0.0529772900204963,-0.0508079539001352,-0.0542798132796589,-0.0666169445343891,-0.065788859731392,-0.0548955083873108,-0.0532026536257134,-0.0605374970597841,-0.0669508368130695,-0.0549328813472718,-0.0656714480499856,-0.063359880194833,-0.0538267540411748,-0.0524904025957669,-0.0667344062685995,-0.0494405152665275,-0.057027972296337,-0.0560062003081574,-0.0566683468784877,-0.0647633871466509,-0.0580479852196485,-0.0555992239993313],"y":[0.038891673454119,0.0357984593168777,0.0376933181870344,0.0363435400287314,0.0382013219091017,0.0450807546017575,0.0444998464308033,0.038727608142556,0.0377465972652728,0.0417218111667468,0.0450298276125442,0.0386430071381795,0.0446432622701469,0.0433510098756955,0.0378838304501787,0.0370323514602239,0.0449185400082439,0.0354997330018897,0.0400113046102533,0.0376584169677684,0.0393042565794872,0.0439521940602111,0.0407474710133996,0.0390747649943793],"text":["Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP9_1594","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP8_2564","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP8_1765","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1820"],"mode":"markers","marker":{"color":"rgba(99,117,103,1)","size":11,"line":{"color":"rgba(99,117,103,1)"}},"type":"scatter","name":"Group3","textfont":{"color":"rgba(99,117,103,1)"},"error_y":{"color":"rgba(99,117,103,1)"},"error_x":{"color":"rgba(99,117,103,1)"},"line":{"color":"rgba(99,117,103,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0434551193461075,-0.0440119096217065,-0.0399228764628456,-0.0489393244675798,-0.0502864259851041,-0.0432656515097961,-0.0439361020406452,-0.0380838068658365],"y":[0.0656574296298447,0.0666361772862748,0.0616148942264521,0.0736715015871354,0.0758371033620305,0.066342355973377,0.066769515617457,0.05843827729274],"text":["Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP6_1542","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP3_2750","Pacuta_HTAC_TP4_1581","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP8_1184"],"mode":"markers","marker":{"color":"rgba(61,130,147,1)","size":11,"line":{"color":"rgba(61,130,147,1)"}},"type":"scatter","name":"Group4","textfont":{"color":"rgba(61,130,147,1)"},"error_y":{"color":"rgba(61,130,147,1)"},"error_x":{"color":"rgba(61,130,147,1)"},"line":{"color":"rgba(61,130,147,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0286497253181373,-0.0313480591605319,-0.0353122275587051,-0.0296475088815271,-0.0292667754932113,-0.034906885846489,-0.0317708367532473,-0.0305946110757402,-0.0301183870569937],"y":[0.0334348851609058,0.0363276861368901,0.0398587145676375,0.0345427533837624,0.0341622646790375,0.0400449816337064,0.0361162644449065,0.0356456397522978,0.034652324116181],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP9_1141","Pacuta_ATHC_TP5_2212","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486"],"mode":"markers","marker":{"color":"rgba(171,44,134,1)","size":11,"line":{"color":"rgba(171,44,134,1)"}},"type":"scatter","name":"Group5","textfont":{"color":"rgba(171,44,134,1)"},"error_y":{"color":"rgba(171,44,134,1)"},"error_x":{"color":"rgba(171,44,134,1)"},"line":{"color":"rgba(171,44,134,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0141721983410679,-0.0151714379070437,-0.0139519932200089,-0.016592815347768,-0.0149478616218807,-0.0155820722003208,-0.0150702673132995,-0.0139914618091777,-0.0130039834849391,-0.0133459406237415,-0.0135906657687692,-0.0154339080103487,-0.015037908195721,-0.0138220742642933,-0.0125063351701625,-0.0138832854326729,-0.014005066225258,-0.0158242716653922,-0.0120979003176351,-0.0143714790497787,-0.016970040020823,-0.0148482807932271,-0.0130764844655136,-0.013175418685665,-0.0135233172500207,-0.0161506084134536,-0.0145942781351581],"y":[0.119857949067155,0.126396446741563,0.117848344084191,0.139138216975891,0.128934984023181,0.130014403463006,0.126931839716419,0.117330413703393,0.109067432275582,0.111649741024461,0.114733882804252,0.128028352404739,0.122996202826022,0.11392920902865,0.105871094365514,0.115690553894148,0.117020770402239,0.131450105991684,0.10151615784482,0.120741803381205,0.138619382533828,0.123158001840636,0.110978149179797,0.111771631077098,0.115608658627416,0.1338581901647,0.122469690993247],"text":["Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP7_1047","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP8_1329","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP5_1168"],"mode":"markers","marker":{"color":"rgba(89,50,89,1)","size":11,"line":{"color":"rgba(89,50,89,1)"}},"type":"scatter","name":"Group6","textfont":{"color":"rgba(89,50,89,1)"},"error_y":{"color":"rgba(89,50,89,1)"},"error_x":{"color":"rgba(89,50,89,1)"},"line":{"color":"rgba(89,50,89,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.00503019000597857,-0.00317282439360773],"y":[0.0438526441761239,0.0617470546541416],"text":["Pacuta_ATAC_TP5_1059","Pacuta_HTHC_TP1_2210"],"mode":"markers","marker":{"color":"rgba(80,80,80,1)","size":11,"line":{"color":"rgba(80,80,80,1)"}},"type":"scatter","name":"Group7","textfont":{"color":"rgba(80,80,80,1)"},"error_y":{"color":"rgba(80,80,80,1)"},"error_x":{"color":"rgba(80,80,80,1)"},"line":{"color":"rgba(80,80,80,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0233961522618102,0.023006134220501],"y":[0.0450378986431911,0.0380658669829281],"text":["Pacuta_HTAC_TP1_1653","Pacuta_HTHC_TP6_1721"],"mode":"markers","marker":{"color":"rgba(225,203,225,1)","size":11,"line":{"color":"rgba(225,203,225,1)"}},"type":"scatter","name":"Group8","textfont":{"color":"rgba(225,203,225,1)"},"error_y":{"color":"rgba(225,203,225,1)"},"error_x":{"color":"rgba(225,203,225,1)"},"line":{"color":"rgba(225,203,225,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.154894980634615,0.161053367851533,0.152144974532457,0.171715648562726,0.175244033300999,0.157517536503936,0.183223822885767,0.182778233394646,0.16955059539825,0.150680735162812,0.1693036033091,0.170153185483835],"y":[-0.0253168473894569,-0.0268033546552437,-0.0248612283848494,-0.0288369875858683,-0.0291969559762794,-0.0258821985725571,-0.0310620564088931,-0.0309166043081994,-0.0281593681506104,-0.0248070886199878,-0.0282825607559514,-0.0284205330230914],"text":["SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301"],"mode":"markers","marker":{"color":"rgba(147,183,174,1)","size":11,"line":{"color":"rgba(147,183,174,1)"}},"type":"scatter","name":"SRA-Singapore-Kusu_Island","textfont":{"color":"rgba(147,183,174,1)"},"error_y":{"color":"rgba(147,183,174,1)"},"error_x":{"color":"rgba(147,183,174,1)"},"line":{"color":"rgba(147,183,174,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.124297238215673,0.139090679154586,0.148070045033372,0.138637472398615,0.141137021025991,0.129582672435747,0.136474617842345,0.142118722308991],"y":[-0.0168533578805978,-0.0194068505725478,-0.020879225470935,-0.0192348733177437,-0.0195358550356349,-0.0175623856721748,-0.0187452249554802,-0.0196976035577411],"text":["SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729"],"mode":"markers","marker":{"color":"rgba(198,148,198,1)","size":11,"line":{"color":"rgba(198,148,198,1)"}},"type":"scatter","name":"SRA-Singapore-Raffles_Lighthouse","textfont":{"color":"rgba(198,148,198,1)"},"error_y":{"color":"rgba(198,148,198,1)"},"error_x":{"color":"rgba(198,148,198,1)"},"line":{"color":"rgba(198,148,198,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.161888831821276,0.118098904122466,0.142469351443509,0.127946670119986,0.157242608766759,0.158031300807886,0.149052590377245,0.160577516053076,0.150051028568593,0.159535523896547,0.144628975677042,0.146569422546426],"y":[-0.026574951462903,-0.0186009035741342,-0.02287882910178,-0.0205858463179394,-0.0256525064662516,-0.0257290284438211,-0.0241818449354243,-0.0261922429240915,-0.0241529853757617,-0.0259599704707734,-0.0233252245314268,-0.023680041314706],"text":["SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"color":"rgba(181,122,50,1)","size":11,"line":{"color":"rgba(181,122,50,1)"}},"type":"scatter","name":"SRA-Singapore-St_Johns_Island","textfont":{"color":"rgba(181,122,50,1)"},"error_y":{"color":"rgba(181,122,50,1)"},"error_x":{"color":"rgba(181,122,50,1)"},"line":{"color":"rgba(181,122,50,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0312741071016752,-0.00236082011067008,-0.0247156670280202,-0.0465076056632526,-0.0397385664719065,-0.0502114329043837],"y":[0.0190685776447495,0.0472475287788579,-0.0308594474938219,-0.0714225158658498,-0.00980969864080465,-0.0394084112097653],"text":["Pacuta_ATAC_TP7_1445","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP6_1254","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP5_1415"],"mode":"markers","marker":{"color":"rgba(182,118,182,1)","size":11,"line":{"color":"rgba(182,118,182,1)"}},"type":"scatter","name":"Ungroup","textfont":{"color":"rgba(182,118,182,1)"},"error_y":{"color":"rgba(182,118,182,1)"},"error_x":{"color":"rgba(182,118,182,1)"},"line":{"color":"rgba(182,118,182,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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

```
## Warning: Duplicate levels detected

## Warning: Duplicate levels detected
```

```{=html}
<div id="htmlwidget-09d0726f2dadf2a49101" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-09d0726f2dadf2a49101">{"x":{"visdat":{"14f559373de0":["function () ","plotlyVisDat"]},"cur_data":"14f559373de0","attrs":{"14f559373de0":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305","SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729","SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301","SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#9d529e","#bc80bd","#e1cbe1","#1b9e77","#d95f02","#01665e","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000","#9d529e","#bc80bd","#e1cbe1"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Reef","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0509810994843474,-0.0538333961763806,-0.0609230079724621,-0.0155820722003208,-0.0551353549725535,-0.0526178381174508,-0.0577909828311013,-0.015037908195721,-0.0498387256099919,-0.0125063351701625,-0.0557205653119484,-0.063359880194833,0.0233961522618102,-0.0492011872391827,-0.0486896121303684,-0.0397385664719065,-0.0649283315075732,-0.0577912527954553,-0.0552832363622027,-0.0521838737667409,-0.0589156035430733],"y":[0.0357984593168777,-0.0873118374741984,-0.119875579377021,0.130014403463006,-0.108392763386094,-0.104354620550549,-0.114482523820953,0.122996202826022,-0.097747252794264,0.105871094365514,-0.1094417660821,0.0433510098756955,0.0450378986431911,-0.0802508105936718,-0.0797514299002811,-0.00980969864080465,-0.126933139530872,-0.114163388756836,-0.109494120805786,-0.102511337341296,-0.115866359564381],"text":["Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP11_1777","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_2873","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP7_2072","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202"],"mode":"markers","marker":{"color":"rgba(205,142,106,1)","size":11,"line":{"color":"rgba(205,142,106,1)"}},"type":"scatter","name":"HIMB","textfont":{"color":"rgba(205,142,106,1)"},"error_y":{"color":"rgba(205,142,106,1)"},"error_x":{"color":"rgba(205,142,106,1)"},"line":{"color":"rgba(205,142,106,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.055223595654752,-0.0434551193461075,-0.065788859731392,-0.0440119096217065,-0.0151714379070437,-0.0149478616218807,-0.0399228764628456,-0.0540851552531239,-0.0135906657687692,-0.0605374970597841,-0.0656714480499856,-0.0502864259851041,-0.0531368267822988,-0.0143714790497787,-0.016970040020823,-0.0693887345846079,-0.0560062003081574,-0.013175418685665,-0.0432656515097961,0.023006134220501,-0.0580479852196485,-0.0439361020406452,-0.0380838068658365,-0.0643050275114951],"y":[0.038891673454119,0.0656574296298447,0.0444998464308033,0.0666361772862748,0.126396446741563,0.128934984023181,0.0616148942264521,-0.105571273752684,0.114733882804252,0.0417218111667468,0.0446432622701469,0.0758371033620305,-0.105032403033994,0.120741803381205,0.138619382533828,-0.135196923550645,0.0376584169677684,0.111771631077098,0.066342355973377,0.0380658669829281,0.0407474710133996,0.066769515617457,0.05843827729274,-0.12594205247791],"text":["Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_2409","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709"],"mode":"markers","marker":{"color":"rgba(228,147,53,1)","size":11,"line":{"color":"rgba(228,147,53,1)"}},"type":"scatter","name":"Lilipuna.Fringe","textfont":{"color":"rgba(228,147,53,1)"},"error_y":{"color":"rgba(228,147,53,1)"},"error_x":{"color":"rgba(228,147,53,1)"},"line":{"color":"rgba(228,147,53,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0573236031218694,-0.0509358823936688,0.00503019000597857,-0.0312741071016752,-0.0515835767663872,-0.016592815347768,-0.0150702673132995,-0.0130039834849391,-0.0133459406237415,-0.0636132014351463,-0.0301183870569937,-0.0148482807932271,-0.0502114329043837,-0.0555992239993313],"y":[-0.113241364175914,-0.100660992690893,0.0438526441761239,0.0190685776447495,-0.086337688274275,0.139138216975891,0.126931839716419,0.109067432275582,0.111649741024461,-0.125009952400077,0.034652324116181,0.123158001840636,-0.0394084112097653,0.0390747649943793],"text":["Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP8_1051","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP9_1486","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP7_1820"],"mode":"markers","marker":{"color":"rgba(170,154,194,1)","size":11,"line":{"color":"rgba(170,154,194,1)"}},"type":"scatter","name":"Reef.11.13","textfont":{"color":"rgba(170,154,194,1)"},"error_y":{"color":"rgba(170,154,194,1)"},"error_x":{"color":"rgba(170,154,194,1)"},"line":{"color":"rgba(170,154,194,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.059013649486611,-0.0529772900204963,-0.0542798132796589,-0.0141721983410679,-0.0553564465445297,-0.0505389556753465,-0.0247156670280202,-0.0622248029799049,-0.0494836920981077,-0.0669508368130695,-0.0514635109927413,-0.0538267540411748,-0.060109103379336,-0.0667344062685995,-0.0608057931386368,-0.0506396577851495,-0.0547479723576256,-0.0135233172500207,-0.0557120353660638,-0.0524243967376363],"y":[-0.115412262806254,0.0376933181870344,0.0382013219091017,0.119857949067155,-0.109261012068346,-0.0992654111619984,-0.0308594474938219,-0.123132070705555,-0.083815992999526,0.0450298276125442,-0.101529294398627,0.0378838304501787,-0.119040080078323,0.0449185400082439,-0.120336263801797,-0.0994384469639176,-0.107091606485858,0.115608658627416,-0.109774013430729,-0.102708465114943],"text":["Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP9_1451","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP8_2513","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"color":"rgba(241,222,202,1)","size":11,"line":{"color":"rgba(241,222,202,1)"}},"type":"scatter","name":"Reef.18","textfont":{"color":"rgba(241,222,202,1)"},"error_y":{"color":"rgba(241,222,202,1)"},"error_x":{"color":"rgba(241,222,202,1)"},"line":{"color":"rgba(241,222,202,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0508079539001352,-0.0353122275587051,-0.0532026536257134,-0.00236082011067008,-0.0139519932200089,-0.0139914618091777,-0.034906885846489,-0.0154339080103487,-0.0138220742642933,-0.0549328813472718,-0.014005066225258,-0.0158242716653922,-0.0120979003176351,-0.0317708367532473,-0.0305946110757402,-0.0161506084134536,-0.0145942781351581],"y":[0.0363435400287314,0.0398587145676375,0.0377465972652728,0.0472475287788579,0.117848344084191,0.117330413703393,0.0400449816337064,0.128028352404739,0.11392920902865,0.0386430071381795,0.117020770402239,0.131450105991684,0.10151615784482,0.0361162644449065,0.0356456397522978,0.1338581901647,0.122469690993247],"text":["Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP9_1594","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP9_1302","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP5_1168"],"mode":"markers","marker":{"color":"rgba(159,84,160,1)","size":11,"line":{"color":"rgba(159,84,160,1)"}},"type":"scatter","name":"Reef.35.36","textfont":{"color":"rgba(159,84,160,1)"},"error_y":{"color":"rgba(159,84,160,1)"},"error_x":{"color":"rgba(159,84,160,1)"},"line":{"color":"rgba(159,84,160,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0286497253181373,-0.0489450962004277,-0.0313480591605319,-0.0632956521039723,-0.0666169445343891,-0.0548955083873108,-0.0296475088815271,-0.0507270406143958,-0.0292667754932113,-0.0489393244675798,-0.0611219706757994,-0.0465076056632526,-0.0625987047830899,-0.0524904025957669,-0.0138832854326729,-0.0530080003490253,-0.0494405152665275,-0.057027972296337,-0.0566683468784877,-0.0130764844655136,-0.00317282439360773,-0.0555465740849365,-0.0647633871466509],"y":[0.0334348851609058,-0.0969941362589877,0.0363276861368901,-0.123706890147924,0.0450807546017575,0.038727608142556,0.0345427533837624,-0.0991740693413415,0.0341622646790375,0.0736715015871354,-0.102383402709373,-0.0714225158658498,-0.122431189919647,0.0370323514602239,0.115690553894148,-0.0891606333044414,0.0354997330018897,0.0400113046102533,0.0393042565794872,0.110978149179797,0.0617470546541416,-0.110231462529279,0.0439521940602111],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP6_1254","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP8_1765","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP4_2195"],"mode":"markers","marker":{"color":"rgba(84,39,136,1)","size":11,"line":{"color":"rgba(84,39,136,1)"}},"type":"scatter","name":"Reef.42.43","textfont":{"color":"rgba(84,39,136,1)"},"error_y":{"color":"rgba(84,39,136,1)"},"error_x":{"color":"rgba(84,39,136,1)"},"line":{"color":"rgba(84,39,136,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.154894980634615,0.161053367851533,0.152144974532457,0.171715648562726,0.175244033300999,0.157517536503936,0.183223822885767,0.182778233394646,0.16955059539825,0.150680735162812,0.1693036033091,0.170153185483835],"y":[-0.0253168473894569,-0.0268033546552437,-0.0248612283848494,-0.0288369875858683,-0.0291969559762794,-0.0258821985725571,-0.0310620564088931,-0.0309166043081994,-0.0281593681506104,-0.0248070886199878,-0.0282825607559514,-0.0284205330230914],"text":["SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301"],"mode":"markers","marker":{"color":"rgba(147,183,174,1)","size":11,"line":{"color":"rgba(147,183,174,1)"}},"type":"scatter","name":"SRA-Singapore-Kusu_Island","textfont":{"color":"rgba(147,183,174,1)"},"error_y":{"color":"rgba(147,183,174,1)"},"error_x":{"color":"rgba(147,183,174,1)"},"line":{"color":"rgba(147,183,174,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.124297238215673,0.139090679154586,0.148070045033372,0.138637472398615,0.141137021025991,0.129582672435747,0.136474617842345,0.142118722308991],"y":[-0.0168533578805978,-0.0194068505725478,-0.020879225470935,-0.0192348733177437,-0.0195358550356349,-0.0175623856721748,-0.0187452249554802,-0.0196976035577411],"text":["SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729"],"mode":"markers","marker":{"color":"rgba(198,148,198,1)","size":11,"line":{"color":"rgba(198,148,198,1)"}},"type":"scatter","name":"SRA-Singapore-Raffles_Lighthouse","textfont":{"color":"rgba(198,148,198,1)"},"error_y":{"color":"rgba(198,148,198,1)"},"error_x":{"color":"rgba(198,148,198,1)"},"line":{"color":"rgba(198,148,198,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.161888831821276,0.118098904122466,0.142469351443509,0.127946670119986,0.157242608766759,0.158031300807886,0.149052590377245,0.160577516053076,0.150051028568593,0.159535523896547,0.144628975677042,0.146569422546426],"y":[-0.026574951462903,-0.0186009035741342,-0.02287882910178,-0.0205858463179394,-0.0256525064662516,-0.0257290284438211,-0.0241818449354243,-0.0261922429240915,-0.0241529853757617,-0.0259599704707734,-0.0233252245314268,-0.023680041314706],"text":["SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"color":"rgba(181,122,50,1)","size":11,"line":{"color":"rgba(181,122,50,1)"}},"type":"scatter","name":"SRA-Singapore-St_Johns_Island","textfont":{"color":"rgba(181,122,50,1)"},"error_y":{"color":"rgba(181,122,50,1)"},"error_x":{"color":"rgba(181,122,50,1)"},"line":{"color":"rgba(181,122,50,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-3be607a16a9303dd4a81" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-3be607a16a9303dd4a81">{"x":{"visdat":{"14f539bf6ae2":["function () ","plotlyVisDat"]},"cur_data":"14f539bf6ae2","attrs":{"14f539bf6ae2":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305","SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729","SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301","SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#9d529e","#bc80bd","#e1cbe1","#1b9e77","#d95f02","#01665e","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000","#9d529e","#bc80bd","#e1cbe1"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Ploidy","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0286497253181373,-0.055223595654752,-0.0509810994843474,-0.0434551193461075,-0.0489450962004277,-0.0538333961763806,-0.059013649486611,-0.0313480591605319,-0.0529772900204963,-0.0573236031218694,-0.0508079539001352,-0.0509358823936688,-0.0542798132796589,-0.0632956521039723,-0.0666169445343891,0.00503019000597857,-0.065788859731392,-0.0548955083873108,-0.0141721983410679,-0.0353122275587051,-0.0440119096217065,-0.0151714379070437,-0.0312741071016752,-0.0553564465445297,-0.0515835767663872,-0.0296475088815271,-0.0507270406143958,-0.0292667754932113,-0.0532026536257134,-0.0505389556753465,-0.00236082011067008,-0.0139519932200089,-0.0247156670280202,-0.0609230079724621,-0.016592815347768,-0.0149478616218807,-0.0155820722003208,-0.0399228764628456,-0.0150702673132995,-0.0139914618091777,-0.0622248029799049,-0.0489393244675798,-0.0540851552531239,-0.0611219706757994,-0.0130039834849391,-0.0133459406237415,-0.034906885846489,-0.0551353549725535,-0.0465076056632526,-0.0526178381174508,-0.0135906657687692,-0.0154339080103487,-0.0605374970597841,-0.0577909828311013,-0.0494836920981077,-0.0669508368130695,-0.015037908195721,-0.0514635109927413,-0.0498387256099919,-0.0138220742642933,-0.0125063351701625,-0.0549328813472718,-0.0557205653119484,-0.0656714480499856,-0.0625987047830899,-0.063359880194833,0.0233961522618102,-0.0538267540411748,-0.060109103379336,-0.0524904025957669,-0.0492011872391827,-0.0138832854326729,-0.0502864259851041,-0.0667344062685995,-0.0530080003490253,-0.014005066225258,-0.0158242716653922,-0.0636132014351463,-0.0120979003176351,-0.0531368267822988,-0.0317708367532473,-0.0143714790497787,-0.0608057931386368,-0.0486896121303684,-0.016970040020823,-0.0494405152665275,-0.0506396577851495,-0.0305946110757402,-0.0301183870569937,-0.0693887345846079,-0.057027972296337,-0.0560062003081574,-0.0547479723576256,-0.0148482807932271,-0.0397385664719065,-0.0566683468784877,-0.0130764844655136,-0.00317282439360773,-0.0555465740849365,-0.013175418685665,-0.0135233172500207,-0.0161506084134536,-0.0432656515097961,-0.0647633871466509,-0.0145942781351581,-0.0502114329043837,-0.0649283315075732,-0.0577912527954553,-0.0552832363622027,0.023006134220501,-0.0580479852196485,-0.0439361020406452,-0.0555992239993313,-0.0380838068658365,-0.0643050275114951,-0.0557120353660638,-0.0521838737667409,-0.0589156035430733,-0.0524243967376363,0.124297238215673,0.139090679154586,0.148070045033372,0.138637472398615,0.141137021025991,0.129582672435747,0.136474617842345,0.142118722308991,0.154894980634615,0.161053367851533,0.152144974532457,0.171715648562726,0.175244033300999,0.157517536503936,0.183223822885767,0.182778233394646,0.16955059539825,0.150680735162812,0.1693036033091,0.170153185483835,0.161888831821276,0.118098904122466,0.142469351443509,0.127946670119986,0.157242608766759,0.158031300807886,0.149052590377245,0.160577516053076,0.150051028568593,0.159535523896547,0.144628975677042,0.146569422546426],"y":[0.0334348851609058,0.038891673454119,0.0357984593168777,0.0656574296298447,-0.0969941362589877,-0.0873118374741984,-0.115412262806254,0.0363276861368901,0.0376933181870344,-0.113241364175914,0.0363435400287314,-0.100660992690893,0.0382013219091017,-0.123706890147924,0.0450807546017575,0.0438526441761239,0.0444998464308033,0.038727608142556,0.119857949067155,0.0398587145676375,0.0666361772862748,0.126396446741563,0.0190685776447495,-0.109261012068346,-0.086337688274275,0.0345427533837624,-0.0991740693413415,0.0341622646790375,0.0377465972652728,-0.0992654111619984,0.0472475287788579,0.117848344084191,-0.0308594474938219,-0.119875579377021,0.139138216975891,0.128934984023181,0.130014403463006,0.0616148942264521,0.126931839716419,0.117330413703393,-0.123132070705555,0.0736715015871354,-0.105571273752684,-0.102383402709373,0.109067432275582,0.111649741024461,0.0400449816337064,-0.108392763386094,-0.0714225158658498,-0.104354620550549,0.114733882804252,0.128028352404739,0.0417218111667468,-0.114482523820953,-0.083815992999526,0.0450298276125442,0.122996202826022,-0.101529294398627,-0.097747252794264,0.11392920902865,0.105871094365514,0.0386430071381795,-0.1094417660821,0.0446432622701469,-0.122431189919647,0.0433510098756955,0.0450378986431911,0.0378838304501787,-0.119040080078323,0.0370323514602239,-0.0802508105936718,0.115690553894148,0.0758371033620305,0.0449185400082439,-0.0891606333044414,0.117020770402239,0.131450105991684,-0.125009952400077,0.10151615784482,-0.105032403033994,0.0361162644449065,0.120741803381205,-0.120336263801797,-0.0797514299002811,0.138619382533828,0.0354997330018897,-0.0994384469639176,0.0356456397522978,0.034652324116181,-0.135196923550645,0.0400113046102533,0.0376584169677684,-0.107091606485858,0.123158001840636,-0.00980969864080465,0.0393042565794872,0.110978149179797,0.0617470546541416,-0.110231462529279,0.111771631077098,0.115608658627416,0.1338581901647,0.066342355973377,0.0439521940602111,0.122469690993247,-0.0394084112097653,-0.126933139530872,-0.114163388756836,-0.109494120805786,0.0380658669829281,0.0407474710133996,0.066769515617457,0.0390747649943793,0.05843827729274,-0.12594205247791,-0.109774013430729,-0.102511337341296,-0.115866359564381,-0.102708465114943,-0.0168533578805978,-0.0194068505725478,-0.020879225470935,-0.0192348733177437,-0.0195358550356349,-0.0175623856721748,-0.0187452249554802,-0.0196976035577411,-0.0253168473894569,-0.0268033546552437,-0.0248612283848494,-0.0288369875858683,-0.0291969559762794,-0.0258821985725571,-0.0310620564088931,-0.0309166043081994,-0.0281593681506104,-0.0248070886199878,-0.0282825607559514,-0.0284205330230914,-0.026574951462903,-0.0186009035741342,-0.02287882910178,-0.0205858463179394,-0.0256525064662516,-0.0257290284438211,-0.0241818449354243,-0.0261922429240915,-0.0241529853757617,-0.0259599704707734,-0.0233252245314268,-0.023680041314706],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305","SRR6914151","SRR6914609","SRR6914908","SRR6934388","SRR6934542","SRR6935629","SRR6942678","SRR6942729","SRR6951423","SRR6951744","SRR6952431","SRR6963586","SRR6963878","SRR6963891","SRR6964364","SRR6986864","SRR6987146","SRR7039808","SRR7040514","SRR7041301","SRR7042978","SRR7043013","SRR7043704","SRR7046161","SRR7055829","SRR7058378","SRR7058566","SRR7058606","SRR7058616","SRR7059699","SRR7062295","SRR7062350"],"mode":"markers","marker":{"colorbar":{"title":"ploidy","ticklen":2},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(235,157,69,1)"],["0.0833333333333335","rgba(184,94,12,1)"],["0.125","rgba(217,201,206,1)"],["0.166666666666667","rgba(249,223,191,1)"],["0.208333333333333","rgba(176,158,193,1)"],["0.25","rgba(157,98,167,1)"],["0.291666666666667","rgba(179,115,180,1)"],["0.333333333333333","rgba(213,178,213,1)"],["0.375","rgba(120,177,158,1)"],["0.416666666666667","rgba(168,128,61,1)"],["0.458333333333333","rgba(132,104,62,1)"],["0.5","rgba(157,119,64,1)"],["0.541666666666667","rgba(231,88,78,1)"],["0.583333333333333","rgba(163,48,137,1)"],["0.625","rgba(104,102,119,1)"],["0.666666666666667","rgba(67,146,96,1)"],["0.708333333333333","rgba(113,112,139,1)"],["0.75","rgba(190,80,37,1)"],["0.791666666666667","rgba(212,66,50,1)"],["0.833333333333333","rgba(106,106,106,1)"],["0.875","rgba(24,16,24,1)"],["0.916666666666667","rgba(160,86,161,1)"],["0.958333333333333","rgba(190,131,190,1)"],["1","rgba(225,203,225,1)"]],"showscale":false,"color":[2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,2,3,3,2,2,3,2,2,3,3,2,3,2,3,3,2,2,2,3,2,2,2,3,2,2,3,3,3,3,2,2,2,3,3,3,2,2,3,3,3,3,2,3,3,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,3,2,2,3,2,3,2,2,3,3,2,3,3,2,2,3,3,3,3,2,3,3,2,2,3,2,2,2,3,3,2,2,3,3,3,2,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2],"size":11,"line":{"colorbar":{"title":"","ticklen":2},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(235,157,69,1)"],["0.0833333333333335","rgba(184,94,12,1)"],["0.125","rgba(217,201,206,1)"],["0.166666666666667","rgba(249,223,191,1)"],["0.208333333333333","rgba(176,158,193,1)"],["0.25","rgba(157,98,167,1)"],["0.291666666666667","rgba(179,115,180,1)"],["0.333333333333333","rgba(213,178,213,1)"],["0.375","rgba(120,177,158,1)"],["0.416666666666667","rgba(168,128,61,1)"],["0.458333333333333","rgba(132,104,62,1)"],["0.5","rgba(157,119,64,1)"],["0.541666666666667","rgba(231,88,78,1)"],["0.583333333333333","rgba(163,48,137,1)"],["0.625","rgba(104,102,119,1)"],["0.666666666666667","rgba(67,146,96,1)"],["0.708333333333333","rgba(113,112,139,1)"],["0.75","rgba(190,80,37,1)"],["0.791666666666667","rgba(212,66,50,1)"],["0.833333333333333","rgba(106,106,106,1)"],["0.875","rgba(24,16,24,1)"],["0.916666666666667","rgba(160,86,161,1)"],["0.958333333333333","rgba(190,131,190,1)"],["1","rgba(225,203,225,1)"]],"showscale":false,"color":[2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,2,3,3,2,2,3,2,2,3,3,2,3,2,3,3,2,2,2,3,2,2,2,3,2,2,3,3,3,3,2,2,2,3,3,3,2,2,3,3,3,3,2,3,3,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,3,2,2,3,2,3,2,2,3,3,2,3,3,2,2,3,3,3,3,2,3,3,2,2,3,2,2,2,3,3,2,2,3,3,3,2,3,3,3,3,3,3,3,3,3,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]}},"type":"scatter","xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0693887345846079,0.183223822885767],"y":[-0.135196923550645,0.139138216975891],"type":"scatter","mode":"markers","opacity":0,"hoverinfo":"none","showlegend":false,"marker":{"colorbar":{"title":"ploidy","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(235,157,69,1)"],["0.0833333333333335","rgba(184,94,12,1)"],["0.125","rgba(217,201,206,1)"],["0.166666666666667","rgba(249,223,191,1)"],["0.208333333333333","rgba(176,158,193,1)"],["0.25","rgba(157,98,167,1)"],["0.291666666666667","rgba(179,115,180,1)"],["0.333333333333333","rgba(213,178,213,1)"],["0.375","rgba(120,177,158,1)"],["0.416666666666667","rgba(168,128,61,1)"],["0.458333333333333","rgba(132,104,62,1)"],["0.5","rgba(157,119,64,1)"],["0.541666666666667","rgba(231,88,78,1)"],["0.583333333333333","rgba(163,48,137,1)"],["0.625","rgba(104,102,119,1)"],["0.666666666666667","rgba(67,146,96,1)"],["0.708333333333333","rgba(113,112,139,1)"],["0.75","rgba(190,80,37,1)"],["0.791666666666667","rgba(212,66,50,1)"],["0.833333333333333","rgba(106,106,106,1)"],["0.875","rgba(24,16,24,1)"],["0.916666666666667","rgba(160,86,161,1)"],["0.958333333333333","rgba(190,131,190,1)"],["1","rgba(225,203,225,1)"]],"showscale":true,"color":[2,3],"line":{"color":"rgba(255,127,14,1)"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
