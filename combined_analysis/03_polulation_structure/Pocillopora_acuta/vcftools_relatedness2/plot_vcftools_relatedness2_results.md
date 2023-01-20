---
title: "Plot `vcftools --relatedness2` results for *P. acuta* RNA-seq samples from this study and SRA"
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
samples.info <- read.table("../../../samples_Pacuta_ALL.annotations.txt", header=T, comment.char='')
rownames(samples.info) <- samples.info$sample
samples.info
```

```
##                                      sample species          treatment
## Pacuta_ATAC_TP11_1777 Pacuta_ATAC_TP11_1777  Pacuta               ATAC
## Pacuta_ATAC_TP8_1051   Pacuta_ATAC_TP8_1051  Pacuta               ATAC
## Pacuta_ATHC_TP4_2733   Pacuta_ATHC_TP4_2733  Pacuta               ATHC
## Pacuta_ATHC_TP8_1459   Pacuta_ATHC_TP8_1459  Pacuta               ATHC
## Pacuta_HTAC_TP3_1642   Pacuta_HTAC_TP3_1642  Pacuta               HTAC
## Pacuta_HTAC_TP4_1767   Pacuta_HTAC_TP4_1767  Pacuta               HTAC
## Pacuta_HTAC_TP7_2072   Pacuta_HTAC_TP7_2072  Pacuta               HTAC
## Pacuta_ATAC_TP11_1103 Pacuta_ATAC_TP11_1103  Pacuta               ATAC
## Pacuta_ATAC_TP11_2306 Pacuta_ATAC_TP11_2306  Pacuta               ATAC
## Pacuta_ATAC_TP3_1041   Pacuta_ATAC_TP3_1041  Pacuta               ATAC
## Pacuta_ATAC_TP3_1637   Pacuta_ATAC_TP3_1637  Pacuta               ATAC
## Pacuta_ATAC_TP4_1762   Pacuta_ATAC_TP4_1762  Pacuta               ATAC
## Pacuta_ATAC_TP7_2413   Pacuta_ATAC_TP7_2413  Pacuta               ATAC
## Pacuta_ATAC_TP8_2012   Pacuta_ATAC_TP8_2012  Pacuta               ATAC
## Pacuta_ATAC_TP9_2357   Pacuta_ATAC_TP9_2357  Pacuta               ATAC
## Pacuta_ATHC_TP11_1147 Pacuta_ATHC_TP11_1147  Pacuta               ATHC
## Pacuta_ATHC_TP3_2534   Pacuta_ATHC_TP3_2534  Pacuta               ATHC
## Pacuta_ATHC_TP4_1220   Pacuta_ATHC_TP4_1220  Pacuta               ATHC
## Pacuta_ATHC_TP5_2877   Pacuta_ATHC_TP5_2877  Pacuta               ATHC
## Pacuta_ATHC_TP6_2870   Pacuta_ATHC_TP6_2870  Pacuta               ATHC
## Pacuta_ATHC_TP7_2878   Pacuta_ATHC_TP7_2878  Pacuta               ATHC
## Pacuta_ATHC_TP9_1451   Pacuta_ATHC_TP9_1451  Pacuta               ATHC
## Pacuta_ATHC_TP9_2873   Pacuta_ATHC_TP9_2873  Pacuta               ATHC
## Pacuta_HTAC_TP10_2064 Pacuta_HTAC_TP10_2064  Pacuta               HTAC
## Pacuta_HTAC_TP11_1596 Pacuta_HTAC_TP11_1596  Pacuta               HTAC
## Pacuta_HTAC_TP1_2414   Pacuta_HTAC_TP1_2414  Pacuta               HTAC
## Pacuta_HTAC_TP5_1707   Pacuta_HTAC_TP5_1707  Pacuta               HTAC
## Pacuta_HTAC_TP6_1466   Pacuta_HTAC_TP6_1466  Pacuta               HTAC
## Pacuta_HTAC_TP7_1728   Pacuta_HTAC_TP7_1728  Pacuta               HTAC
## Pacuta_HTAC_TP8_2513   Pacuta_HTAC_TP8_2513  Pacuta               HTAC
## Pacuta_HTAC_TP9_1696   Pacuta_HTAC_TP9_1696  Pacuta               HTAC
## Pacuta_HTHC_TP10_2300 Pacuta_HTHC_TP10_2300  Pacuta               HTHC
## Pacuta_HTHC_TP3_1227   Pacuta_HTHC_TP3_1227  Pacuta               HTHC
## Pacuta_HTHC_TP5_2087   Pacuta_HTHC_TP5_2087  Pacuta               HTHC
## Pacuta_HTHC_TP6_1138   Pacuta_HTHC_TP6_1138  Pacuta               HTHC
## Pacuta_HTHC_TP6_1595   Pacuta_HTHC_TP6_1595  Pacuta               HTHC
## Pacuta_HTHC_TP8_1709   Pacuta_HTHC_TP8_1709  Pacuta               HTHC
## Pacuta_HTHC_TP8_2304   Pacuta_HTHC_TP8_2304  Pacuta               HTHC
## Pacuta_HTHC_TP9_1131   Pacuta_HTHC_TP9_1131  Pacuta               HTHC
## Pacuta_HTHC_TP9_2202   Pacuta_HTHC_TP9_2202  Pacuta               HTHC
## Pacuta_HTHC_TP9_2305   Pacuta_HTHC_TP9_2305  Pacuta               HTHC
## Pacuta_ATAC_TP10_1559 Pacuta_ATAC_TP10_1559  Pacuta               ATAC
## Pacuta_ATAC_TP10_1641 Pacuta_ATAC_TP10_1641  Pacuta               ATAC
## Pacuta_ATAC_TP1_2363   Pacuta_ATAC_TP1_2363  Pacuta               ATAC
## Pacuta_ATAC_TP3_1471   Pacuta_ATAC_TP3_1471  Pacuta               ATAC
## Pacuta_ATAC_TP4_1060   Pacuta_ATAC_TP4_1060  Pacuta               ATAC
## Pacuta_ATAC_TP4_2002   Pacuta_ATAC_TP4_2002  Pacuta               ATAC
## Pacuta_ATAC_TP5_1563   Pacuta_ATAC_TP5_1563  Pacuta               ATAC
## Pacuta_ATAC_TP5_1757   Pacuta_ATAC_TP5_1757  Pacuta               ATAC
## Pacuta_ATAC_TP9_1594   Pacuta_ATAC_TP9_1594  Pacuta               ATAC
## Pacuta_ATHC_TP7_2409   Pacuta_ATHC_TP7_2409  Pacuta               ATHC
## Pacuta_ATHC_TP8_2564   Pacuta_ATHC_TP8_2564  Pacuta               ATHC
## Pacuta_HTAC_TP10_1536 Pacuta_HTAC_TP10_1536  Pacuta               HTAC
## Pacuta_HTAC_TP11_1582 Pacuta_HTAC_TP11_1582  Pacuta               HTAC
## Pacuta_HTAC_TP11_1647 Pacuta_HTAC_TP11_1647  Pacuta               HTAC
## Pacuta_HTAC_TP1_2005   Pacuta_HTAC_TP1_2005  Pacuta               HTAC
## Pacuta_HTAC_TP3_1617   Pacuta_HTAC_TP3_1617  Pacuta               HTAC
## Pacuta_HTAC_TP4_1701   Pacuta_HTAC_TP4_1701  Pacuta               HTAC
## Pacuta_HTAC_TP8_1765   Pacuta_HTAC_TP8_1765  Pacuta               HTAC
## Pacuta_HTHC_TP10_1238 Pacuta_HTHC_TP10_1238  Pacuta               HTHC
## Pacuta_HTHC_TP10_1732 Pacuta_HTHC_TP10_1732  Pacuta               HTHC
## Pacuta_HTHC_TP1_1239   Pacuta_HTHC_TP1_1239  Pacuta               HTHC
## Pacuta_HTHC_TP4_2195   Pacuta_HTHC_TP4_2195  Pacuta               HTHC
## Pacuta_HTHC_TP7_1090   Pacuta_HTHC_TP7_1090  Pacuta               HTHC
## Pacuta_HTHC_TP7_1820   Pacuta_HTHC_TP7_1820  Pacuta               HTHC
## Pacuta_ATAC_TP1_1043   Pacuta_ATAC_TP1_1043  Pacuta               ATAC
## Pacuta_ATAC_TP6_1542   Pacuta_ATAC_TP6_1542  Pacuta               ATAC
## Pacuta_ATHC_TP1_2743   Pacuta_ATHC_TP1_2743  Pacuta               ATHC
## Pacuta_ATHC_TP3_2750   Pacuta_ATHC_TP3_2750  Pacuta               ATHC
## Pacuta_HTAC_TP4_1581   Pacuta_HTAC_TP4_1581  Pacuta               HTAC
## Pacuta_HTHC_TP4_1343   Pacuta_HTHC_TP4_1343  Pacuta               HTHC
## Pacuta_HTHC_TP7_1427   Pacuta_HTHC_TP7_1427  Pacuta               HTHC
## Pacuta_HTHC_TP8_1184   Pacuta_HTHC_TP8_1184  Pacuta               HTHC
## Pacuta_ATAC_TP10_1159 Pacuta_ATAC_TP10_1159  Pacuta               ATAC
## Pacuta_ATAC_TP1_1775   Pacuta_ATAC_TP1_1775  Pacuta               ATAC
## Pacuta_ATAC_TP6_1468   Pacuta_ATAC_TP6_1468  Pacuta               ATAC
## Pacuta_ATAC_TP8_1755   Pacuta_ATAC_TP8_1755  Pacuta               ATAC
## Pacuta_ATAC_TP9_1141   Pacuta_ATAC_TP9_1141  Pacuta               ATAC
## Pacuta_ATHC_TP5_2212   Pacuta_ATHC_TP5_2212  Pacuta               ATHC
## Pacuta_HTAC_TP6_1744   Pacuta_HTAC_TP6_1744  Pacuta               HTAC
## Pacuta_HTAC_TP9_1302   Pacuta_HTAC_TP9_1302  Pacuta               HTAC
## Pacuta_HTAC_TP9_1486   Pacuta_HTAC_TP9_1486  Pacuta               HTAC
## Pacuta_ATAC_TP6_1050   Pacuta_ATAC_TP6_1050  Pacuta               ATAC
## Pacuta_ATAC_TP7_1047   Pacuta_ATAC_TP7_1047  Pacuta               ATAC
## Pacuta_ATHC_TP10_2197 Pacuta_ATHC_TP10_2197  Pacuta               ATHC
## Pacuta_ATHC_TP1_1207   Pacuta_ATHC_TP1_1207  Pacuta               ATHC
## Pacuta_ATHC_TP11_2668 Pacuta_ATHC_TP11_2668  Pacuta               ATHC
## Pacuta_ATHC_TP11_2879 Pacuta_ATHC_TP11_2879  Pacuta               ATHC
## Pacuta_ATHC_TP1_2977   Pacuta_ATHC_TP1_2977  Pacuta               ATHC
## Pacuta_ATHC_TP3_1219   Pacuta_ATHC_TP3_1219  Pacuta               ATHC
## Pacuta_ATHC_TP4_2993   Pacuta_ATHC_TP4_2993  Pacuta               ATHC
## Pacuta_ATHC_TP5_1296   Pacuta_ATHC_TP5_1296  Pacuta               ATHC
## Pacuta_ATHC_TP6_2999   Pacuta_ATHC_TP6_2999  Pacuta               ATHC
## Pacuta_ATHC_TP7_1281   Pacuta_ATHC_TP7_1281  Pacuta               ATHC
## Pacuta_ATHC_TP8_2861   Pacuta_ATHC_TP8_2861  Pacuta               ATHC
## Pacuta_ATHC_TP9_2979   Pacuta_ATHC_TP9_2979  Pacuta               ATHC
## Pacuta_HTAC_TP10_1225 Pacuta_HTAC_TP10_1225  Pacuta               HTAC
## Pacuta_HTAC_TP3_2026   Pacuta_HTAC_TP3_2026  Pacuta               HTAC
## Pacuta_HTAC_TP5_1303   Pacuta_HTAC_TP5_1303  Pacuta               HTAC
## Pacuta_HTAC_TP5_1571   Pacuta_HTAC_TP5_1571  Pacuta               HTAC
## Pacuta_HTAC_TP6_1330   Pacuta_HTAC_TP6_1330  Pacuta               HTAC
## Pacuta_HTAC_TP7_1487   Pacuta_HTAC_TP7_1487  Pacuta               HTAC
## Pacuta_HTAC_TP8_1329   Pacuta_HTAC_TP8_1329  Pacuta               HTAC
## Pacuta_HTHC_TP11_1416 Pacuta_HTHC_TP11_1416  Pacuta               HTHC
## Pacuta_HTHC_TP1_1676   Pacuta_HTHC_TP1_1676  Pacuta               HTHC
## Pacuta_HTHC_TP3_1418   Pacuta_HTHC_TP3_1418  Pacuta               HTHC
## Pacuta_HTHC_TP3_2527   Pacuta_HTHC_TP3_2527  Pacuta               HTHC
## Pacuta_HTHC_TP4_1169   Pacuta_HTHC_TP4_1169  Pacuta               HTHC
## Pacuta_HTHC_TP5_1168   Pacuta_HTHC_TP5_1168  Pacuta               HTHC
## Pacuta_ATAC_TP5_1059   Pacuta_ATAC_TP5_1059  Pacuta               ATAC
## Pacuta_HTHC_TP1_2210   Pacuta_HTHC_TP1_2210  Pacuta               HTHC
## Pacuta_HTAC_TP1_1653   Pacuta_HTAC_TP1_1653  Pacuta               HTAC
## Pacuta_HTHC_TP6_1721   Pacuta_HTHC_TP6_1721  Pacuta               HTHC
## Pacuta_ATAC_TP7_1445   Pacuta_ATAC_TP7_1445  Pacuta               ATAC
## Pacuta_ATHC_TP10_1205 Pacuta_ATHC_TP10_1205  Pacuta               ATHC
## Pacuta_ATHC_TP10_2550 Pacuta_ATHC_TP10_2550  Pacuta               ATHC
## Pacuta_ATHC_TP6_1254   Pacuta_ATHC_TP6_1254  Pacuta               ATHC
## Pacuta_HTHC_TP11_2185 Pacuta_HTHC_TP11_2185  Pacuta               HTHC
## Pacuta_HTHC_TP5_1415   Pacuta_HTHC_TP5_1415  Pacuta               HTHC
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
## Pacuta_ATAC_TP11_1777      TP11      1777                             HIMB
## Pacuta_ATAC_TP8_1051        TP8      1051                       Reef.11.13
## Pacuta_ATHC_TP4_2733        TP4      2733                       Reef.42.43
## Pacuta_ATHC_TP8_1459        TP8      1459                          Reef.18
## Pacuta_HTAC_TP3_1642        TP3      1642                             HIMB
## Pacuta_HTAC_TP4_1767        TP4      1767                       Reef.42.43
## Pacuta_HTAC_TP7_2072        TP7      2072                             HIMB
## Pacuta_ATAC_TP11_1103      TP11      1103                       Reef.42.43
## Pacuta_ATAC_TP11_2306      TP11      2306                          Reef.18
## Pacuta_ATAC_TP3_1041        TP3      1041                       Reef.11.13
## Pacuta_ATAC_TP3_1637        TP3      1637                       Reef.11.13
## Pacuta_ATAC_TP4_1762        TP4      1762                       Reef.42.43
## Pacuta_ATAC_TP7_2413        TP7      2413                          Reef.18
## Pacuta_ATAC_TP8_2012        TP8      2012                       Reef.42.43
## Pacuta_ATAC_TP9_2357        TP9      2357                          Reef.18
## Pacuta_ATHC_TP11_1147      TP11      1147                             HIMB
## Pacuta_ATHC_TP3_2534        TP3      2534                          Reef.18
## Pacuta_ATHC_TP4_1220        TP4      1220                  Lilipuna.Fringe
## Pacuta_ATHC_TP5_2877        TP5      2877                             HIMB
## Pacuta_ATHC_TP6_2870        TP6      2870                             HIMB
## Pacuta_ATHC_TP7_2878        TP7      2878                             HIMB
## Pacuta_ATHC_TP9_1451        TP9      1451                          Reef.18
## Pacuta_ATHC_TP9_2873        TP9      2873                             HIMB
## Pacuta_HTAC_TP10_2064      TP10      2064                             HIMB
## Pacuta_HTAC_TP11_1596      TP11      1596                       Reef.42.43
## Pacuta_HTAC_TP1_2414        TP1      2414                          Reef.18
## Pacuta_HTAC_TP5_1707        TP5      1707                       Reef.11.13
## Pacuta_HTAC_TP6_1466        TP6      1466                  Lilipuna.Fringe
## Pacuta_HTAC_TP7_1728        TP7      1728                          Reef.18
## Pacuta_HTAC_TP8_2513        TP8      2513                          Reef.18
## Pacuta_HTAC_TP9_1696        TP9      1696                  Lilipuna.Fringe
## Pacuta_HTHC_TP10_2300      TP10      2300                          Reef.18
## Pacuta_HTHC_TP3_1227        TP3      1227                       Reef.42.43
## Pacuta_HTHC_TP5_2087        TP5      2087                             HIMB
## Pacuta_HTHC_TP6_1138        TP6      1138                             HIMB
## Pacuta_HTHC_TP6_1595        TP6      1595                             HIMB
## Pacuta_HTHC_TP8_1709        TP8      1709                  Lilipuna.Fringe
## Pacuta_HTHC_TP8_2304        TP8      2304                          Reef.18
## Pacuta_HTHC_TP9_1131        TP9      1131                             HIMB
## Pacuta_HTHC_TP9_2202        TP9      2202                             HIMB
## Pacuta_HTHC_TP9_2305        TP9      2305                          Reef.18
## Pacuta_ATAC_TP10_1559      TP10      1559                  Lilipuna.Fringe
## Pacuta_ATAC_TP10_1641      TP10      1641                             HIMB
## Pacuta_ATAC_TP1_2363        TP1      2363                          Reef.18
## Pacuta_ATAC_TP3_1471        TP3      1471                       Reef.35.36
## Pacuta_ATAC_TP4_1060        TP4      1060                          Reef.18
## Pacuta_ATAC_TP4_2002        TP4      2002                       Reef.42.43
## Pacuta_ATAC_TP5_1563        TP5      1563                  Lilipuna.Fringe
## Pacuta_ATAC_TP5_1757        TP5      1757                       Reef.42.43
## Pacuta_ATAC_TP9_1594        TP9      1594                       Reef.35.36
## Pacuta_ATHC_TP7_2409        TP7      2409                  Lilipuna.Fringe
## Pacuta_ATHC_TP8_2564        TP8      2564                          Reef.18
## Pacuta_HTAC_TP10_1536      TP10      1536                       Reef.35.36
## Pacuta_HTAC_TP11_1582      TP11      1582                  Lilipuna.Fringe
## Pacuta_HTAC_TP11_1647      TP11      1647                             HIMB
## Pacuta_HTAC_TP1_2005        TP1      2005                          Reef.18
## Pacuta_HTAC_TP3_1617        TP3      1617                       Reef.42.43
## Pacuta_HTAC_TP4_1701        TP4      1701                          Reef.18
## Pacuta_HTAC_TP8_1765        TP8      1765                       Reef.42.43
## Pacuta_HTHC_TP10_1238      TP10      1238                       Reef.42.43
## Pacuta_HTHC_TP10_1732      TP10      1732                  Lilipuna.Fringe
## Pacuta_HTHC_TP1_1239        TP1      1239                       Reef.42.43
## Pacuta_HTHC_TP4_2195        TP4      2195                       Reef.42.43
## Pacuta_HTHC_TP7_1090        TP7      1090                  Lilipuna.Fringe
## Pacuta_HTHC_TP7_1820        TP7      1820                       Reef.11.13
## Pacuta_ATAC_TP1_1043        TP1      1043                  Lilipuna.Fringe
## Pacuta_ATAC_TP6_1542        TP6      1542                  Lilipuna.Fringe
## Pacuta_ATHC_TP1_2743        TP1      2743                  Lilipuna.Fringe
## Pacuta_ATHC_TP3_2750        TP3      2750                       Reef.42.43
## Pacuta_HTAC_TP4_1581        TP4      1581                  Lilipuna.Fringe
## Pacuta_HTHC_TP4_1343        TP4      1343                  Lilipuna.Fringe
## Pacuta_HTHC_TP7_1427        TP7      1427                  Lilipuna.Fringe
## Pacuta_HTHC_TP8_1184        TP8      1184                  Lilipuna.Fringe
## Pacuta_ATAC_TP10_1159      TP10      1159                       Reef.42.43
## Pacuta_ATAC_TP1_1775        TP1      1775                       Reef.42.43
## Pacuta_ATAC_TP6_1468        TP6      1468                       Reef.35.36
## Pacuta_ATAC_TP8_1755        TP8      1755                       Reef.42.43
## Pacuta_ATAC_TP9_1141        TP9      1141                       Reef.42.43
## Pacuta_ATHC_TP5_2212        TP5      2212                       Reef.35.36
## Pacuta_HTAC_TP6_1744        TP6      1744                       Reef.35.36
## Pacuta_HTAC_TP9_1302        TP9      1302                       Reef.35.36
## Pacuta_HTAC_TP9_1486        TP9      1486                       Reef.11.13
## Pacuta_ATAC_TP6_1050        TP6      1050                          Reef.18
## Pacuta_ATAC_TP7_1047        TP7      1047                  Lilipuna.Fringe
## Pacuta_ATHC_TP10_2197      TP10      2197                       Reef.35.36
## Pacuta_ATHC_TP1_1207        TP1      1207                       Reef.11.13
## Pacuta_ATHC_TP11_2668      TP11      2668                  Lilipuna.Fringe
## Pacuta_ATHC_TP11_2879      TP11      2879                             HIMB
## Pacuta_ATHC_TP1_2977        TP1      2977                       Reef.11.13
## Pacuta_ATHC_TP3_1219        TP3      1219                       Reef.35.36
## Pacuta_ATHC_TP4_2993        TP4      2993                       Reef.11.13
## Pacuta_ATHC_TP5_1296        TP5      1296                       Reef.11.13
## Pacuta_ATHC_TP6_2999        TP6      2999                  Lilipuna.Fringe
## Pacuta_ATHC_TP7_1281        TP7      1281                       Reef.35.36
## Pacuta_ATHC_TP8_2861        TP8      2861                             HIMB
## Pacuta_ATHC_TP9_2979        TP9      2979                       Reef.35.36
## Pacuta_HTAC_TP10_1225      TP10      1225                             HIMB
## Pacuta_HTAC_TP3_2026        TP3      2026                       Reef.42.43
## Pacuta_HTAC_TP5_1303        TP5      1303                       Reef.35.36
## Pacuta_HTAC_TP5_1571        TP5      1571                       Reef.35.36
## Pacuta_HTAC_TP6_1330        TP6      1330                       Reef.35.36
## Pacuta_HTAC_TP7_1487        TP7      1487                  Lilipuna.Fringe
## Pacuta_HTAC_TP8_1329        TP8      1329                  Lilipuna.Fringe
## Pacuta_HTHC_TP11_1416      TP11      1416                       Reef.11.13
## Pacuta_HTHC_TP1_1676        TP1      1676                       Reef.42.43
## Pacuta_HTHC_TP3_1418        TP3      1418                  Lilipuna.Fringe
## Pacuta_HTHC_TP3_2527        TP3      2527                          Reef.18
## Pacuta_HTHC_TP4_1169        TP4      1169                       Reef.35.36
## Pacuta_HTHC_TP5_1168        TP5      1168                       Reef.35.36
## Pacuta_ATAC_TP5_1059        TP5      1059                       Reef.11.13
## Pacuta_HTHC_TP1_2210        TP1      2210                       Reef.42.43
## Pacuta_HTAC_TP1_1653        TP1      1653                             HIMB
## Pacuta_HTHC_TP6_1721        TP6      1721                  Lilipuna.Fringe
## Pacuta_ATAC_TP7_1445        TP7      1445                       Reef.11.13
## Pacuta_ATHC_TP10_1205      TP10      1205                       Reef.35.36
## Pacuta_ATHC_TP10_2550      TP10      2550                          Reef.18
## Pacuta_ATHC_TP6_1254        TP6      1254                       Reef.42.43
## Pacuta_HTHC_TP11_2185      TP11      2185                             HIMB
## Pacuta_HTHC_TP5_1415        TP5      1415                       Reef.11.13
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
## Pacuta_ATAC_TP11_1777    #b35806      3      #d95f02
## Pacuta_ATAC_TP8_1051     #fee0b6      3      #d95f02
## Pacuta_ATHC_TP4_2733     #542788      3      #d95f02
## Pacuta_ATHC_TP8_1459     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP3_1642     #b35806      3      #d95f02
## Pacuta_HTAC_TP4_1767     #542788      3      #d95f02
## Pacuta_HTAC_TP7_2072     #b35806      3      #d95f02
## Pacuta_ATAC_TP11_1103    #542788      3      #d95f02
## Pacuta_ATAC_TP11_2306    #d8daeb      3      #d95f02
## Pacuta_ATAC_TP3_1041     #fee0b6      3      #d95f02
## Pacuta_ATAC_TP3_1637     #fee0b6      3      #d95f02
## Pacuta_ATAC_TP4_1762     #542788      3      #d95f02
## Pacuta_ATAC_TP7_2413     #d8daeb      3      #d95f02
## Pacuta_ATAC_TP8_2012     #542788      3      #d95f02
## Pacuta_ATAC_TP9_2357     #d8daeb      3      #d95f02
## Pacuta_ATHC_TP11_1147    #b35806      3      #d95f02
## Pacuta_ATHC_TP3_2534     #d8daeb      3      #d95f02
## Pacuta_ATHC_TP4_1220     #f1a340      3      #d95f02
## Pacuta_ATHC_TP5_2877     #b35806      3      #d95f02
## Pacuta_ATHC_TP6_2870     #b35806      3      #d95f02
## Pacuta_ATHC_TP7_2878     #b35806      3      #d95f02
## Pacuta_ATHC_TP9_1451     #d8daeb      3      #d95f02
## Pacuta_ATHC_TP9_2873     #b35806      3      #d95f02
## Pacuta_HTAC_TP10_2064    #b35806      3      #d95f02
## Pacuta_HTAC_TP11_1596    #542788      3      #d95f02
## Pacuta_HTAC_TP1_2414     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP5_1707     #fee0b6      3      #d95f02
## Pacuta_HTAC_TP6_1466     #f1a340      3      #d95f02
## Pacuta_HTAC_TP7_1728     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP8_2513     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP9_1696     #f1a340      3      #d95f02
## Pacuta_HTHC_TP10_2300    #d8daeb      3      #d95f02
## Pacuta_HTHC_TP3_1227     #542788      3      #d95f02
## Pacuta_HTHC_TP5_2087     #b35806      3      #d95f02
## Pacuta_HTHC_TP6_1138     #b35806      3      #d95f02
## Pacuta_HTHC_TP6_1595     #b35806      3      #d95f02
## Pacuta_HTHC_TP8_1709     #f1a340      3      #d95f02
## Pacuta_HTHC_TP8_2304     #d8daeb      3      #d95f02
## Pacuta_HTHC_TP9_1131     #b35806      3      #d95f02
## Pacuta_HTHC_TP9_2202     #b35806      3      #d95f02
## Pacuta_HTHC_TP9_2305     #d8daeb      3      #d95f02
## Pacuta_ATAC_TP10_1559    #f1a340      3      #d95f02
## Pacuta_ATAC_TP10_1641    #b35806      3      #d95f02
## Pacuta_ATAC_TP1_2363     #d8daeb      3      #d95f02
## Pacuta_ATAC_TP3_1471     #998ec3      3      #d95f02
## Pacuta_ATAC_TP4_1060     #d8daeb      3      #d95f02
## Pacuta_ATAC_TP4_2002     #542788      3      #d95f02
## Pacuta_ATAC_TP5_1563     #f1a340      3      #d95f02
## Pacuta_ATAC_TP5_1757     #542788      3      #d95f02
## Pacuta_ATAC_TP9_1594     #998ec3      3      #d95f02
## Pacuta_ATHC_TP7_2409     #f1a340      3      #d95f02
## Pacuta_ATHC_TP8_2564     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP10_1536    #998ec3      3      #d95f02
## Pacuta_HTAC_TP11_1582    #f1a340      3      #d95f02
## Pacuta_HTAC_TP11_1647    #b35806      3      #d95f02
## Pacuta_HTAC_TP1_2005     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP3_1617     #542788      3      #d95f02
## Pacuta_HTAC_TP4_1701     #d8daeb      3      #d95f02
## Pacuta_HTAC_TP8_1765     #542788      3      #d95f02
## Pacuta_HTHC_TP10_1238    #542788      3      #d95f02
## Pacuta_HTHC_TP10_1732    #f1a340      3      #d95f02
## Pacuta_HTHC_TP1_1239     #542788      3      #d95f02
## Pacuta_HTHC_TP4_2195     #542788      3      #d95f02
## Pacuta_HTHC_TP7_1090     #f1a340      3      #d95f02
## Pacuta_HTHC_TP7_1820     #fee0b6      3      #d95f02
## Pacuta_ATAC_TP1_1043     #f1a340      3      #d95f02
## Pacuta_ATAC_TP6_1542     #f1a340      3      #d95f02
## Pacuta_ATHC_TP1_2743     #f1a340      3      #d95f02
## Pacuta_ATHC_TP3_2750     #542788      3      #d95f02
## Pacuta_HTAC_TP4_1581     #f1a340      3      #d95f02
## Pacuta_HTHC_TP4_1343     #f1a340      3      #d95f02
## Pacuta_HTHC_TP7_1427     #f1a340      3      #d95f02
## Pacuta_HTHC_TP8_1184     #f1a340      3      #d95f02
## Pacuta_ATAC_TP10_1159    #542788      2      #1b9e77
## Pacuta_ATAC_TP1_1775     #542788      2      #1b9e77
## Pacuta_ATAC_TP6_1468     #998ec3      2      #1b9e77
## Pacuta_ATAC_TP8_1755     #542788      2      #1b9e77
## Pacuta_ATAC_TP9_1141     #542788      2      #1b9e77
## Pacuta_ATHC_TP5_2212     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP6_1744     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP9_1302     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP9_1486     #fee0b6      2      #1b9e77
## Pacuta_ATAC_TP6_1050     #d8daeb      2      #1b9e77
## Pacuta_ATAC_TP7_1047     #f1a340      2      #1b9e77
## Pacuta_ATHC_TP10_2197    #998ec3      2      #1b9e77
## Pacuta_ATHC_TP1_1207     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP11_2668    #f1a340      2      #1b9e77
## Pacuta_ATHC_TP11_2879    #b35806      2      #1b9e77
## Pacuta_ATHC_TP1_2977     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP3_1219     #998ec3      2      #1b9e77
## Pacuta_ATHC_TP4_2993     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP5_1296     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP6_2999     #f1a340      2      #1b9e77
## Pacuta_ATHC_TP7_1281     #998ec3      2      #1b9e77
## Pacuta_ATHC_TP8_2861     #b35806      2      #1b9e77
## Pacuta_ATHC_TP9_2979     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP10_1225    #b35806      2      #1b9e77
## Pacuta_HTAC_TP3_2026     #542788      2      #1b9e77
## Pacuta_HTAC_TP5_1303     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP5_1571     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP6_1330     #998ec3      2      #1b9e77
## Pacuta_HTAC_TP7_1487     #f1a340      2      #1b9e77
## Pacuta_HTAC_TP8_1329     #f1a340      2      #1b9e77
## Pacuta_HTHC_TP11_1416    #fee0b6      2      #1b9e77
## Pacuta_HTHC_TP1_1676     #542788      2      #1b9e77
## Pacuta_HTHC_TP3_1418     #f1a340      2      #1b9e77
## Pacuta_HTHC_TP3_2527     #d8daeb      2      #1b9e77
## Pacuta_HTHC_TP4_1169     #998ec3      2      #1b9e77
## Pacuta_HTHC_TP5_1168     #998ec3      2      #1b9e77
## Pacuta_ATAC_TP5_1059     #fee0b6      2      #1b9e77
## Pacuta_HTHC_TP1_2210     #542788      2      #1b9e77
## Pacuta_HTAC_TP1_1653     #b35806      2      #1b9e77
## Pacuta_HTHC_TP6_1721     #f1a340      2      #1b9e77
## Pacuta_ATAC_TP7_1445     #fee0b6      2      #1b9e77
## Pacuta_ATHC_TP10_1205    #998ec3      2      #1b9e77
## Pacuta_ATHC_TP10_2550    #d8daeb      2      #1b9e77
## Pacuta_ATHC_TP6_1254     #542788      3      #d95f02
## Pacuta_HTHC_TP11_2185    #b35806      3      #d95f02
## Pacuta_HTHC_TP5_1415     #fee0b6      2      #1b9e77
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
## Pacuta_ATAC_TP11_1777                           Group1     #1f78b4
## Pacuta_ATAC_TP8_1051                            Group1     #1f78b4
## Pacuta_ATHC_TP4_2733                            Group1     #1f78b4
## Pacuta_ATHC_TP8_1459                            Group1     #1f78b4
## Pacuta_HTAC_TP3_1642                            Group1     #1f78b4
## Pacuta_HTAC_TP4_1767                            Group1     #1f78b4
## Pacuta_HTAC_TP7_2072                            Group1     #1f78b4
## Pacuta_ATAC_TP11_1103                           Group2     #33a02c
## Pacuta_ATAC_TP11_2306                           Group2     #33a02c
## Pacuta_ATAC_TP3_1041                            Group2     #33a02c
## Pacuta_ATAC_TP3_1637                            Group2     #33a02c
## Pacuta_ATAC_TP4_1762                            Group2     #33a02c
## Pacuta_ATAC_TP7_2413                            Group2     #33a02c
## Pacuta_ATAC_TP8_2012                            Group2     #33a02c
## Pacuta_ATAC_TP9_2357                            Group2     #33a02c
## Pacuta_ATHC_TP11_1147                           Group2     #33a02c
## Pacuta_ATHC_TP3_2534                            Group2     #33a02c
## Pacuta_ATHC_TP4_1220                            Group2     #33a02c
## Pacuta_ATHC_TP5_2877                            Group2     #33a02c
## Pacuta_ATHC_TP6_2870                            Group2     #33a02c
## Pacuta_ATHC_TP7_2878                            Group2     #33a02c
## Pacuta_ATHC_TP9_1451                            Group2     #33a02c
## Pacuta_ATHC_TP9_2873                            Group2     #33a02c
## Pacuta_HTAC_TP10_2064                           Group2     #33a02c
## Pacuta_HTAC_TP11_1596                           Group2     #33a02c
## Pacuta_HTAC_TP1_2414                            Group2     #33a02c
## Pacuta_HTAC_TP5_1707                            Group2     #33a02c
## Pacuta_HTAC_TP6_1466                            Group2     #33a02c
## Pacuta_HTAC_TP7_1728                            Group2     #33a02c
## Pacuta_HTAC_TP8_2513                            Group2     #33a02c
## Pacuta_HTAC_TP9_1696                            Group2     #33a02c
## Pacuta_HTHC_TP10_2300                           Group2     #33a02c
## Pacuta_HTHC_TP3_1227                            Group2     #33a02c
## Pacuta_HTHC_TP5_2087                            Group2     #33a02c
## Pacuta_HTHC_TP6_1138                            Group2     #33a02c
## Pacuta_HTHC_TP6_1595                            Group2     #33a02c
## Pacuta_HTHC_TP8_1709                            Group2     #33a02c
## Pacuta_HTHC_TP8_2304                            Group2     #33a02c
## Pacuta_HTHC_TP9_1131                            Group2     #33a02c
## Pacuta_HTHC_TP9_2202                            Group2     #33a02c
## Pacuta_HTHC_TP9_2305                            Group2     #33a02c
## Pacuta_ATAC_TP10_1559                           Group3     #c51b7d
## Pacuta_ATAC_TP10_1641                           Group3     #c51b7d
## Pacuta_ATAC_TP1_2363                            Group3     #c51b7d
## Pacuta_ATAC_TP3_1471                            Group3     #c51b7d
## Pacuta_ATAC_TP4_1060                            Group3     #c51b7d
## Pacuta_ATAC_TP4_2002                            Group3     #c51b7d
## Pacuta_ATAC_TP5_1563                            Group3     #c51b7d
## Pacuta_ATAC_TP5_1757                            Group3     #c51b7d
## Pacuta_ATAC_TP9_1594                            Group3     #c51b7d
## Pacuta_ATHC_TP7_2409                            Group3     #c51b7d
## Pacuta_ATHC_TP8_2564                            Group3     #c51b7d
## Pacuta_HTAC_TP10_1536                           Group3     #c51b7d
## Pacuta_HTAC_TP11_1582                           Group3     #c51b7d
## Pacuta_HTAC_TP11_1647                           Group3     #c51b7d
## Pacuta_HTAC_TP1_2005                            Group3     #c51b7d
## Pacuta_HTAC_TP3_1617                            Group3     #c51b7d
## Pacuta_HTAC_TP4_1701                            Group3     #c51b7d
## Pacuta_HTAC_TP8_1765                            Group3     #c51b7d
## Pacuta_HTHC_TP10_1238                           Group3     #c51b7d
## Pacuta_HTHC_TP10_1732                           Group3     #c51b7d
## Pacuta_HTHC_TP1_1239                            Group3     #c51b7d
## Pacuta_HTHC_TP4_2195                            Group3     #c51b7d
## Pacuta_HTHC_TP7_1090                            Group3     #c51b7d
## Pacuta_HTHC_TP7_1820                            Group3     #c51b7d
## Pacuta_ATAC_TP1_1043                            Group4     #6a3d9a
## Pacuta_ATAC_TP6_1542                            Group4     #6a3d9a
## Pacuta_ATHC_TP1_2743                            Group4     #6a3d9a
## Pacuta_ATHC_TP3_2750                            Group4     #6a3d9a
## Pacuta_HTAC_TP4_1581                            Group4     #6a3d9a
## Pacuta_HTHC_TP4_1343                            Group4     #6a3d9a
## Pacuta_HTHC_TP7_1427                            Group4     #6a3d9a
## Pacuta_HTHC_TP8_1184                            Group4     #6a3d9a
## Pacuta_ATAC_TP10_1159                           Group5     #ff7f00
## Pacuta_ATAC_TP1_1775                            Group5     #ff7f00
## Pacuta_ATAC_TP6_1468                            Group5     #ff7f00
## Pacuta_ATAC_TP8_1755                            Group5     #ff7f00
## Pacuta_ATAC_TP9_1141                            Group5     #ff7f00
## Pacuta_ATHC_TP5_2212                            Group5     #ff7f00
## Pacuta_HTAC_TP6_1744                            Group5     #ff7f00
## Pacuta_HTAC_TP9_1302                            Group5     #ff7f00
## Pacuta_HTAC_TP9_1486                            Group5     #ff7f00
## Pacuta_ATAC_TP6_1050                            Group6     #e31a1c
## Pacuta_ATAC_TP7_1047                            Group6     #e31a1c
## Pacuta_ATHC_TP10_2197                           Group6     #e31a1c
## Pacuta_ATHC_TP1_1207                            Group6     #e31a1c
## Pacuta_ATHC_TP11_2668                           Group6     #e31a1c
## Pacuta_ATHC_TP11_2879                           Group6     #e31a1c
## Pacuta_ATHC_TP1_2977                            Group6     #e31a1c
## Pacuta_ATHC_TP3_1219                            Group6     #e31a1c
## Pacuta_ATHC_TP4_2993                            Group6     #e31a1c
## Pacuta_ATHC_TP5_1296                            Group6     #e31a1c
## Pacuta_ATHC_TP6_2999                            Group6     #e31a1c
## Pacuta_ATHC_TP7_1281                            Group6     #e31a1c
## Pacuta_ATHC_TP8_2861                            Group6     #e31a1c
## Pacuta_ATHC_TP9_2979                            Group6     #e31a1c
## Pacuta_HTAC_TP10_1225                           Group6     #e31a1c
## Pacuta_HTAC_TP3_2026                            Group6     #e31a1c
## Pacuta_HTAC_TP5_1303                            Group6     #e31a1c
## Pacuta_HTAC_TP5_1571                            Group6     #e31a1c
## Pacuta_HTAC_TP6_1330                            Group6     #e31a1c
## Pacuta_HTAC_TP7_1487                            Group6     #e31a1c
## Pacuta_HTAC_TP8_1329                            Group6     #e31a1c
## Pacuta_HTHC_TP11_1416                           Group6     #e31a1c
## Pacuta_HTHC_TP1_1676                            Group6     #e31a1c
## Pacuta_HTHC_TP3_1418                            Group6     #e31a1c
## Pacuta_HTHC_TP3_2527                            Group6     #e31a1c
## Pacuta_HTHC_TP4_1169                            Group6     #e31a1c
## Pacuta_HTHC_TP5_1168                            Group6     #e31a1c
## Pacuta_ATAC_TP5_1059                            Group7     #b15928
## Pacuta_HTHC_TP1_2210                            Group7     #b15928
## Pacuta_HTAC_TP1_1653                            Group8     #000000
## Pacuta_HTHC_TP6_1721                            Group8     #000000
## Pacuta_ATAC_TP7_1445                           Ungroup     #808080
## Pacuta_ATHC_TP10_1205                          Ungroup     #808080
## Pacuta_ATHC_TP10_2550                          Ungroup     #808080
## Pacuta_ATHC_TP6_1254                           Ungroup     #808080
## Pacuta_HTHC_TP11_2185                          Ungroup     #808080
## Pacuta_HTHC_TP5_1415                           Ungroup     #808080
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
