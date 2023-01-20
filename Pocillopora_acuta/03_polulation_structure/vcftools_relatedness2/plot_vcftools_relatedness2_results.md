---
title: "Plot `vcftools --relatedness2` results for *P. acuta* RNA-seq samples"
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
samples.info <- read.table("../../samples_Pacuta.annotations.txt", header=T, comment.char='')
rownames(samples.info) <- samples.info$sample
samples.info
```

```
##                                      sample species treatment timepoint plugid
## Pacuta_ATAC_TP11_1777 Pacuta_ATAC_TP11_1777  Pacuta      ATAC      TP11   1777
## Pacuta_ATAC_TP8_1051   Pacuta_ATAC_TP8_1051  Pacuta      ATAC       TP8   1051
## Pacuta_ATHC_TP4_2733   Pacuta_ATHC_TP4_2733  Pacuta      ATHC       TP4   2733
## Pacuta_ATHC_TP8_1459   Pacuta_ATHC_TP8_1459  Pacuta      ATHC       TP8   1459
## Pacuta_HTAC_TP3_1642   Pacuta_HTAC_TP3_1642  Pacuta      HTAC       TP3   1642
## Pacuta_HTAC_TP4_1767   Pacuta_HTAC_TP4_1767  Pacuta      HTAC       TP4   1767
## Pacuta_HTAC_TP7_2072   Pacuta_HTAC_TP7_2072  Pacuta      HTAC       TP7   2072
## Pacuta_ATAC_TP11_1103 Pacuta_ATAC_TP11_1103  Pacuta      ATAC      TP11   1103
## Pacuta_ATAC_TP11_2306 Pacuta_ATAC_TP11_2306  Pacuta      ATAC      TP11   2306
## Pacuta_ATAC_TP3_1041   Pacuta_ATAC_TP3_1041  Pacuta      ATAC       TP3   1041
## Pacuta_ATAC_TP3_1637   Pacuta_ATAC_TP3_1637  Pacuta      ATAC       TP3   1637
## Pacuta_ATAC_TP4_1762   Pacuta_ATAC_TP4_1762  Pacuta      ATAC       TP4   1762
## Pacuta_ATAC_TP7_2413   Pacuta_ATAC_TP7_2413  Pacuta      ATAC       TP7   2413
## Pacuta_ATAC_TP8_2012   Pacuta_ATAC_TP8_2012  Pacuta      ATAC       TP8   2012
## Pacuta_ATAC_TP9_2357   Pacuta_ATAC_TP9_2357  Pacuta      ATAC       TP9   2357
## Pacuta_ATHC_TP11_1147 Pacuta_ATHC_TP11_1147  Pacuta      ATHC      TP11   1147
## Pacuta_ATHC_TP3_2534   Pacuta_ATHC_TP3_2534  Pacuta      ATHC       TP3   2534
## Pacuta_ATHC_TP4_1220   Pacuta_ATHC_TP4_1220  Pacuta      ATHC       TP4   1220
## Pacuta_ATHC_TP5_2877   Pacuta_ATHC_TP5_2877  Pacuta      ATHC       TP5   2877
## Pacuta_ATHC_TP6_2870   Pacuta_ATHC_TP6_2870  Pacuta      ATHC       TP6   2870
## Pacuta_ATHC_TP7_2878   Pacuta_ATHC_TP7_2878  Pacuta      ATHC       TP7   2878
## Pacuta_ATHC_TP9_1451   Pacuta_ATHC_TP9_1451  Pacuta      ATHC       TP9   1451
## Pacuta_ATHC_TP9_2873   Pacuta_ATHC_TP9_2873  Pacuta      ATHC       TP9   2873
## Pacuta_HTAC_TP10_2064 Pacuta_HTAC_TP10_2064  Pacuta      HTAC      TP10   2064
## Pacuta_HTAC_TP11_1596 Pacuta_HTAC_TP11_1596  Pacuta      HTAC      TP11   1596
## Pacuta_HTAC_TP1_2414   Pacuta_HTAC_TP1_2414  Pacuta      HTAC       TP1   2414
## Pacuta_HTAC_TP5_1707   Pacuta_HTAC_TP5_1707  Pacuta      HTAC       TP5   1707
## Pacuta_HTAC_TP6_1466   Pacuta_HTAC_TP6_1466  Pacuta      HTAC       TP6   1466
## Pacuta_HTAC_TP7_1728   Pacuta_HTAC_TP7_1728  Pacuta      HTAC       TP7   1728
## Pacuta_HTAC_TP8_2513   Pacuta_HTAC_TP8_2513  Pacuta      HTAC       TP8   2513
## Pacuta_HTAC_TP9_1696   Pacuta_HTAC_TP9_1696  Pacuta      HTAC       TP9   1696
## Pacuta_HTHC_TP10_2300 Pacuta_HTHC_TP10_2300  Pacuta      HTHC      TP10   2300
## Pacuta_HTHC_TP3_1227   Pacuta_HTHC_TP3_1227  Pacuta      HTHC       TP3   1227
## Pacuta_HTHC_TP5_2087   Pacuta_HTHC_TP5_2087  Pacuta      HTHC       TP5   2087
## Pacuta_HTHC_TP6_1138   Pacuta_HTHC_TP6_1138  Pacuta      HTHC       TP6   1138
## Pacuta_HTHC_TP6_1595   Pacuta_HTHC_TP6_1595  Pacuta      HTHC       TP6   1595
## Pacuta_HTHC_TP8_1709   Pacuta_HTHC_TP8_1709  Pacuta      HTHC       TP8   1709
## Pacuta_HTHC_TP8_2304   Pacuta_HTHC_TP8_2304  Pacuta      HTHC       TP8   2304
## Pacuta_HTHC_TP9_1131   Pacuta_HTHC_TP9_1131  Pacuta      HTHC       TP9   1131
## Pacuta_HTHC_TP9_2202   Pacuta_HTHC_TP9_2202  Pacuta      HTHC       TP9   2202
## Pacuta_HTHC_TP9_2305   Pacuta_HTHC_TP9_2305  Pacuta      HTHC       TP9   2305
## Pacuta_ATAC_TP10_1559 Pacuta_ATAC_TP10_1559  Pacuta      ATAC      TP10   1559
## Pacuta_ATAC_TP10_1641 Pacuta_ATAC_TP10_1641  Pacuta      ATAC      TP10   1641
## Pacuta_ATAC_TP1_2363   Pacuta_ATAC_TP1_2363  Pacuta      ATAC       TP1   2363
## Pacuta_ATAC_TP3_1471   Pacuta_ATAC_TP3_1471  Pacuta      ATAC       TP3   1471
## Pacuta_ATAC_TP4_1060   Pacuta_ATAC_TP4_1060  Pacuta      ATAC       TP4   1060
## Pacuta_ATAC_TP4_2002   Pacuta_ATAC_TP4_2002  Pacuta      ATAC       TP4   2002
## Pacuta_ATAC_TP5_1563   Pacuta_ATAC_TP5_1563  Pacuta      ATAC       TP5   1563
## Pacuta_ATAC_TP5_1757   Pacuta_ATAC_TP5_1757  Pacuta      ATAC       TP5   1757
## Pacuta_ATAC_TP9_1594   Pacuta_ATAC_TP9_1594  Pacuta      ATAC       TP9   1594
## Pacuta_ATHC_TP7_2409   Pacuta_ATHC_TP7_2409  Pacuta      ATHC       TP7   2409
## Pacuta_ATHC_TP8_2564   Pacuta_ATHC_TP8_2564  Pacuta      ATHC       TP8   2564
## Pacuta_HTAC_TP10_1536 Pacuta_HTAC_TP10_1536  Pacuta      HTAC      TP10   1536
## Pacuta_HTAC_TP11_1582 Pacuta_HTAC_TP11_1582  Pacuta      HTAC      TP11   1582
## Pacuta_HTAC_TP11_1647 Pacuta_HTAC_TP11_1647  Pacuta      HTAC      TP11   1647
## Pacuta_HTAC_TP1_2005   Pacuta_HTAC_TP1_2005  Pacuta      HTAC       TP1   2005
## Pacuta_HTAC_TP3_1617   Pacuta_HTAC_TP3_1617  Pacuta      HTAC       TP3   1617
## Pacuta_HTAC_TP4_1701   Pacuta_HTAC_TP4_1701  Pacuta      HTAC       TP4   1701
## Pacuta_HTAC_TP8_1765   Pacuta_HTAC_TP8_1765  Pacuta      HTAC       TP8   1765
## Pacuta_HTHC_TP10_1238 Pacuta_HTHC_TP10_1238  Pacuta      HTHC      TP10   1238
## Pacuta_HTHC_TP10_1732 Pacuta_HTHC_TP10_1732  Pacuta      HTHC      TP10   1732
## Pacuta_HTHC_TP1_1239   Pacuta_HTHC_TP1_1239  Pacuta      HTHC       TP1   1239
## Pacuta_HTHC_TP4_2195   Pacuta_HTHC_TP4_2195  Pacuta      HTHC       TP4   2195
## Pacuta_HTHC_TP7_1090   Pacuta_HTHC_TP7_1090  Pacuta      HTHC       TP7   1090
## Pacuta_HTHC_TP7_1820   Pacuta_HTHC_TP7_1820  Pacuta      HTHC       TP7   1820
## Pacuta_ATAC_TP1_1043   Pacuta_ATAC_TP1_1043  Pacuta      ATAC       TP1   1043
## Pacuta_ATAC_TP6_1542   Pacuta_ATAC_TP6_1542  Pacuta      ATAC       TP6   1542
## Pacuta_ATHC_TP1_2743   Pacuta_ATHC_TP1_2743  Pacuta      ATHC       TP1   2743
## Pacuta_ATHC_TP3_2750   Pacuta_ATHC_TP3_2750  Pacuta      ATHC       TP3   2750
## Pacuta_HTAC_TP4_1581   Pacuta_HTAC_TP4_1581  Pacuta      HTAC       TP4   1581
## Pacuta_HTHC_TP4_1343   Pacuta_HTHC_TP4_1343  Pacuta      HTHC       TP4   1343
## Pacuta_HTHC_TP7_1427   Pacuta_HTHC_TP7_1427  Pacuta      HTHC       TP7   1427
## Pacuta_HTHC_TP8_1184   Pacuta_HTHC_TP8_1184  Pacuta      HTHC       TP8   1184
## Pacuta_ATAC_TP10_1159 Pacuta_ATAC_TP10_1159  Pacuta      ATAC      TP10   1159
## Pacuta_ATAC_TP1_1775   Pacuta_ATAC_TP1_1775  Pacuta      ATAC       TP1   1775
## Pacuta_ATAC_TP6_1468   Pacuta_ATAC_TP6_1468  Pacuta      ATAC       TP6   1468
## Pacuta_ATAC_TP8_1755   Pacuta_ATAC_TP8_1755  Pacuta      ATAC       TP8   1755
## Pacuta_ATAC_TP9_1141   Pacuta_ATAC_TP9_1141  Pacuta      ATAC       TP9   1141
## Pacuta_ATHC_TP5_2212   Pacuta_ATHC_TP5_2212  Pacuta      ATHC       TP5   2212
## Pacuta_HTAC_TP6_1744   Pacuta_HTAC_TP6_1744  Pacuta      HTAC       TP6   1744
## Pacuta_HTAC_TP9_1302   Pacuta_HTAC_TP9_1302  Pacuta      HTAC       TP9   1302
## Pacuta_HTAC_TP9_1486   Pacuta_HTAC_TP9_1486  Pacuta      HTAC       TP9   1486
## Pacuta_ATAC_TP6_1050   Pacuta_ATAC_TP6_1050  Pacuta      ATAC       TP6   1050
## Pacuta_ATAC_TP7_1047   Pacuta_ATAC_TP7_1047  Pacuta      ATAC       TP7   1047
## Pacuta_ATHC_TP10_2197 Pacuta_ATHC_TP10_2197  Pacuta      ATHC      TP10   2197
## Pacuta_ATHC_TP1_1207   Pacuta_ATHC_TP1_1207  Pacuta      ATHC       TP1   1207
## Pacuta_ATHC_TP11_2668 Pacuta_ATHC_TP11_2668  Pacuta      ATHC      TP11   2668
## Pacuta_ATHC_TP11_2879 Pacuta_ATHC_TP11_2879  Pacuta      ATHC      TP11   2879
## Pacuta_ATHC_TP1_2977   Pacuta_ATHC_TP1_2977  Pacuta      ATHC       TP1   2977
## Pacuta_ATHC_TP3_1219   Pacuta_ATHC_TP3_1219  Pacuta      ATHC       TP3   1219
## Pacuta_ATHC_TP4_2993   Pacuta_ATHC_TP4_2993  Pacuta      ATHC       TP4   2993
## Pacuta_ATHC_TP5_1296   Pacuta_ATHC_TP5_1296  Pacuta      ATHC       TP5   1296
## Pacuta_ATHC_TP6_2999   Pacuta_ATHC_TP6_2999  Pacuta      ATHC       TP6   2999
## Pacuta_ATHC_TP7_1281   Pacuta_ATHC_TP7_1281  Pacuta      ATHC       TP7   1281
## Pacuta_ATHC_TP8_2861   Pacuta_ATHC_TP8_2861  Pacuta      ATHC       TP8   2861
## Pacuta_ATHC_TP9_2979   Pacuta_ATHC_TP9_2979  Pacuta      ATHC       TP9   2979
## Pacuta_HTAC_TP10_1225 Pacuta_HTAC_TP10_1225  Pacuta      HTAC      TP10   1225
## Pacuta_HTAC_TP3_2026   Pacuta_HTAC_TP3_2026  Pacuta      HTAC       TP3   2026
## Pacuta_HTAC_TP5_1303   Pacuta_HTAC_TP5_1303  Pacuta      HTAC       TP5   1303
## Pacuta_HTAC_TP5_1571   Pacuta_HTAC_TP5_1571  Pacuta      HTAC       TP5   1571
## Pacuta_HTAC_TP6_1330   Pacuta_HTAC_TP6_1330  Pacuta      HTAC       TP6   1330
## Pacuta_HTAC_TP7_1487   Pacuta_HTAC_TP7_1487  Pacuta      HTAC       TP7   1487
## Pacuta_HTAC_TP8_1329   Pacuta_HTAC_TP8_1329  Pacuta      HTAC       TP8   1329
## Pacuta_HTHC_TP11_1416 Pacuta_HTHC_TP11_1416  Pacuta      HTHC      TP11   1416
## Pacuta_HTHC_TP1_1676   Pacuta_HTHC_TP1_1676  Pacuta      HTHC       TP1   1676
## Pacuta_HTHC_TP3_1418   Pacuta_HTHC_TP3_1418  Pacuta      HTHC       TP3   1418
## Pacuta_HTHC_TP3_2527   Pacuta_HTHC_TP3_2527  Pacuta      HTHC       TP3   2527
## Pacuta_HTHC_TP4_1169   Pacuta_HTHC_TP4_1169  Pacuta      HTHC       TP4   1169
## Pacuta_HTHC_TP5_1168   Pacuta_HTHC_TP5_1168  Pacuta      HTHC       TP5   1168
## Pacuta_ATAC_TP5_1059   Pacuta_ATAC_TP5_1059  Pacuta      ATAC       TP5   1059
## Pacuta_HTHC_TP1_2210   Pacuta_HTHC_TP1_2210  Pacuta      HTHC       TP1   2210
## Pacuta_HTAC_TP1_1653   Pacuta_HTAC_TP1_1653  Pacuta      HTAC       TP1   1653
## Pacuta_HTHC_TP6_1721   Pacuta_HTHC_TP6_1721  Pacuta      HTHC       TP6   1721
## Pacuta_ATAC_TP7_1445   Pacuta_ATAC_TP7_1445  Pacuta      ATAC       TP7   1445
## Pacuta_ATHC_TP10_1205 Pacuta_ATHC_TP10_1205  Pacuta      ATHC      TP10   1205
## Pacuta_ATHC_TP10_2550 Pacuta_ATHC_TP10_2550  Pacuta      ATHC      TP10   2550
## Pacuta_ATHC_TP6_1254   Pacuta_ATHC_TP6_1254  Pacuta      ATHC       TP6   1254
## Pacuta_HTHC_TP11_2185 Pacuta_HTHC_TP11_2185  Pacuta      HTHC      TP11   2185
## Pacuta_HTHC_TP5_1415   Pacuta_HTHC_TP5_1415  Pacuta      HTHC       TP5   1415
##                                  reef reef_color ploidy ploidy_color   group
## Pacuta_ATAC_TP11_1777            HIMB    #b35806      3      #d95f02  Group1
## Pacuta_ATAC_TP8_1051       Reef.11.13    #fee0b6      3      #d95f02  Group1
## Pacuta_ATHC_TP4_2733       Reef.42.43    #542788      3      #d95f02  Group1
## Pacuta_ATHC_TP8_1459          Reef.18    #d8daeb      3      #d95f02  Group1
## Pacuta_HTAC_TP3_1642             HIMB    #b35806      3      #d95f02  Group1
## Pacuta_HTAC_TP4_1767       Reef.42.43    #542788      3      #d95f02  Group1
## Pacuta_HTAC_TP7_2072             HIMB    #b35806      3      #d95f02  Group1
## Pacuta_ATAC_TP11_1103      Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_ATAC_TP11_2306         Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATAC_TP3_1041       Reef.11.13    #fee0b6      3      #d95f02  Group2
## Pacuta_ATAC_TP3_1637       Reef.11.13    #fee0b6      3      #d95f02  Group2
## Pacuta_ATAC_TP4_1762       Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_ATAC_TP7_2413          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATAC_TP8_2012       Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_ATAC_TP9_2357          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATHC_TP11_1147            HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP3_2534          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATHC_TP4_1220  Lilipuna.Fringe    #f1a340      3      #d95f02  Group2
## Pacuta_ATHC_TP5_2877             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP6_2870             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP7_2878             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP9_1451          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATHC_TP9_2873             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTAC_TP10_2064            HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTAC_TP11_1596      Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_HTAC_TP1_2414          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTAC_TP5_1707       Reef.11.13    #fee0b6      3      #d95f02  Group2
## Pacuta_HTAC_TP6_1466  Lilipuna.Fringe    #f1a340      3      #d95f02  Group2
## Pacuta_HTAC_TP7_1728          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTAC_TP8_2513          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTAC_TP9_1696  Lilipuna.Fringe    #f1a340      3      #d95f02  Group2
## Pacuta_HTHC_TP10_2300         Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTHC_TP3_1227       Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_HTHC_TP5_2087             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP6_1138             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP6_1595             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP8_1709  Lilipuna.Fringe    #f1a340      3      #d95f02  Group2
## Pacuta_HTHC_TP8_2304          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTHC_TP9_1131             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP9_2202             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP9_2305          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATAC_TP10_1559 Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_ATAC_TP10_1641            HIMB    #b35806      3      #d95f02  Group3
## Pacuta_ATAC_TP1_2363          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_ATAC_TP3_1471       Reef.35.36    #998ec3      3      #d95f02  Group3
## Pacuta_ATAC_TP4_1060          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_ATAC_TP4_2002       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_ATAC_TP5_1563  Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_ATAC_TP5_1757       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_ATAC_TP9_1594       Reef.35.36    #998ec3      3      #d95f02  Group3
## Pacuta_ATHC_TP7_2409  Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_ATHC_TP8_2564          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_HTAC_TP10_1536      Reef.35.36    #998ec3      3      #d95f02  Group3
## Pacuta_HTAC_TP11_1582 Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_HTAC_TP11_1647            HIMB    #b35806      3      #d95f02  Group3
## Pacuta_HTAC_TP1_2005          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_HTAC_TP3_1617       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTAC_TP4_1701          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_HTAC_TP8_1765       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTHC_TP10_1238      Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTHC_TP10_1732 Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_HTHC_TP1_1239       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTHC_TP4_2195       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTHC_TP7_1090  Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_HTHC_TP7_1820       Reef.11.13    #fee0b6      3      #d95f02  Group3
## Pacuta_ATAC_TP1_1043  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_ATAC_TP6_1542  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_ATHC_TP1_2743  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_ATHC_TP3_2750       Reef.42.43    #542788      3      #d95f02  Group4
## Pacuta_HTAC_TP4_1581  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_HTHC_TP4_1343  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_HTHC_TP7_1427  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_HTHC_TP8_1184  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_ATAC_TP10_1159      Reef.42.43    #542788      2      #1b9e77  Group5
## Pacuta_ATAC_TP1_1775       Reef.42.43    #542788      2      #1b9e77  Group5
## Pacuta_ATAC_TP6_1468       Reef.35.36    #998ec3      2      #1b9e77  Group5
## Pacuta_ATAC_TP8_1755       Reef.42.43    #542788      2      #1b9e77  Group5
## Pacuta_ATAC_TP9_1141       Reef.42.43    #542788      2      #1b9e77  Group5
## Pacuta_ATHC_TP5_2212       Reef.35.36    #998ec3      2      #1b9e77  Group5
## Pacuta_HTAC_TP6_1744       Reef.35.36    #998ec3      2      #1b9e77  Group5
## Pacuta_HTAC_TP9_1302       Reef.35.36    #998ec3      2      #1b9e77  Group5
## Pacuta_HTAC_TP9_1486       Reef.11.13    #fee0b6      2      #1b9e77  Group5
## Pacuta_ATAC_TP6_1050          Reef.18    #d8daeb      2      #1b9e77  Group6
## Pacuta_ATAC_TP7_1047  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_ATHC_TP10_2197      Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_ATHC_TP1_1207       Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_ATHC_TP11_2668 Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_ATHC_TP11_2879            HIMB    #b35806      2      #1b9e77  Group6
## Pacuta_ATHC_TP1_2977       Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_ATHC_TP3_1219       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_ATHC_TP4_2993       Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_ATHC_TP5_1296       Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_ATHC_TP6_2999  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_ATHC_TP7_1281       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_ATHC_TP8_2861             HIMB    #b35806      2      #1b9e77  Group6
## Pacuta_ATHC_TP9_2979       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTAC_TP10_1225            HIMB    #b35806      2      #1b9e77  Group6
## Pacuta_HTAC_TP3_2026       Reef.42.43    #542788      2      #1b9e77  Group6
## Pacuta_HTAC_TP5_1303       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTAC_TP5_1571       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTAC_TP6_1330       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTAC_TP7_1487  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_HTAC_TP8_1329  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_HTHC_TP11_1416      Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_HTHC_TP1_1676       Reef.42.43    #542788      2      #1b9e77  Group6
## Pacuta_HTHC_TP3_1418  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_HTHC_TP3_2527          Reef.18    #d8daeb      2      #1b9e77  Group6
## Pacuta_HTHC_TP4_1169       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTHC_TP5_1168       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_ATAC_TP5_1059       Reef.11.13    #fee0b6      2      #1b9e77  Group7
## Pacuta_HTHC_TP1_2210       Reef.42.43    #542788      2      #1b9e77  Group7
## Pacuta_HTAC_TP1_1653             HIMB    #b35806      2      #1b9e77  Group8
## Pacuta_HTHC_TP6_1721  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group8
## Pacuta_ATAC_TP7_1445       Reef.11.13    #fee0b6      2      #1b9e77 Ungroup
## Pacuta_ATHC_TP10_1205      Reef.35.36    #998ec3      2      #1b9e77 Ungroup
## Pacuta_ATHC_TP10_2550         Reef.18    #d8daeb      2      #1b9e77 Ungroup
## Pacuta_ATHC_TP6_1254       Reef.42.43    #542788      3      #d95f02 Ungroup
## Pacuta_HTHC_TP11_2185            HIMB    #b35806      3      #d95f02 Ungroup
## Pacuta_HTHC_TP5_1415       Reef.11.13    #fee0b6      2      #1b9e77 Ungroup
##                       group_color
## Pacuta_ATAC_TP11_1777     #1f78b4
## Pacuta_ATAC_TP8_1051      #1f78b4
## Pacuta_ATHC_TP4_2733      #1f78b4
## Pacuta_ATHC_TP8_1459      #1f78b4
## Pacuta_HTAC_TP3_1642      #1f78b4
## Pacuta_HTAC_TP4_1767      #1f78b4
## Pacuta_HTAC_TP7_2072      #1f78b4
## Pacuta_ATAC_TP11_1103     #33a02c
## Pacuta_ATAC_TP11_2306     #33a02c
## Pacuta_ATAC_TP3_1041      #33a02c
## Pacuta_ATAC_TP3_1637      #33a02c
## Pacuta_ATAC_TP4_1762      #33a02c
## Pacuta_ATAC_TP7_2413      #33a02c
## Pacuta_ATAC_TP8_2012      #33a02c
## Pacuta_ATAC_TP9_2357      #33a02c
## Pacuta_ATHC_TP11_1147     #33a02c
## Pacuta_ATHC_TP3_2534      #33a02c
## Pacuta_ATHC_TP4_1220      #33a02c
## Pacuta_ATHC_TP5_2877      #33a02c
## Pacuta_ATHC_TP6_2870      #33a02c
## Pacuta_ATHC_TP7_2878      #33a02c
## Pacuta_ATHC_TP9_1451      #33a02c
## Pacuta_ATHC_TP9_2873      #33a02c
## Pacuta_HTAC_TP10_2064     #33a02c
## Pacuta_HTAC_TP11_1596     #33a02c
## Pacuta_HTAC_TP1_2414      #33a02c
## Pacuta_HTAC_TP5_1707      #33a02c
## Pacuta_HTAC_TP6_1466      #33a02c
## Pacuta_HTAC_TP7_1728      #33a02c
## Pacuta_HTAC_TP8_2513      #33a02c
## Pacuta_HTAC_TP9_1696      #33a02c
## Pacuta_HTHC_TP10_2300     #33a02c
## Pacuta_HTHC_TP3_1227      #33a02c
## Pacuta_HTHC_TP5_2087      #33a02c
## Pacuta_HTHC_TP6_1138      #33a02c
## Pacuta_HTHC_TP6_1595      #33a02c
## Pacuta_HTHC_TP8_1709      #33a02c
## Pacuta_HTHC_TP8_2304      #33a02c
## Pacuta_HTHC_TP9_1131      #33a02c
## Pacuta_HTHC_TP9_2202      #33a02c
## Pacuta_HTHC_TP9_2305      #33a02c
## Pacuta_ATAC_TP10_1559     #c51b7d
## Pacuta_ATAC_TP10_1641     #c51b7d
## Pacuta_ATAC_TP1_2363      #c51b7d
## Pacuta_ATAC_TP3_1471      #c51b7d
## Pacuta_ATAC_TP4_1060      #c51b7d
## Pacuta_ATAC_TP4_2002      #c51b7d
## Pacuta_ATAC_TP5_1563      #c51b7d
## Pacuta_ATAC_TP5_1757      #c51b7d
## Pacuta_ATAC_TP9_1594      #c51b7d
## Pacuta_ATHC_TP7_2409      #c51b7d
## Pacuta_ATHC_TP8_2564      #c51b7d
## Pacuta_HTAC_TP10_1536     #c51b7d
## Pacuta_HTAC_TP11_1582     #c51b7d
## Pacuta_HTAC_TP11_1647     #c51b7d
## Pacuta_HTAC_TP1_2005      #c51b7d
## Pacuta_HTAC_TP3_1617      #c51b7d
## Pacuta_HTAC_TP4_1701      #c51b7d
## Pacuta_HTAC_TP8_1765      #c51b7d
## Pacuta_HTHC_TP10_1238     #c51b7d
## Pacuta_HTHC_TP10_1732     #c51b7d
## Pacuta_HTHC_TP1_1239      #c51b7d
## Pacuta_HTHC_TP4_2195      #c51b7d
## Pacuta_HTHC_TP7_1090      #c51b7d
## Pacuta_HTHC_TP7_1820      #c51b7d
## Pacuta_ATAC_TP1_1043      #6a3d9a
## Pacuta_ATAC_TP6_1542      #6a3d9a
## Pacuta_ATHC_TP1_2743      #6a3d9a
## Pacuta_ATHC_TP3_2750      #6a3d9a
## Pacuta_HTAC_TP4_1581      #6a3d9a
## Pacuta_HTHC_TP4_1343      #6a3d9a
## Pacuta_HTHC_TP7_1427      #6a3d9a
## Pacuta_HTHC_TP8_1184      #6a3d9a
## Pacuta_ATAC_TP10_1159     #ff7f00
## Pacuta_ATAC_TP1_1775      #ff7f00
## Pacuta_ATAC_TP6_1468      #ff7f00
## Pacuta_ATAC_TP8_1755      #ff7f00
## Pacuta_ATAC_TP9_1141      #ff7f00
## Pacuta_ATHC_TP5_2212      #ff7f00
## Pacuta_HTAC_TP6_1744      #ff7f00
## Pacuta_HTAC_TP9_1302      #ff7f00
## Pacuta_HTAC_TP9_1486      #ff7f00
## Pacuta_ATAC_TP6_1050      #e31a1c
## Pacuta_ATAC_TP7_1047      #e31a1c
## Pacuta_ATHC_TP10_2197     #e31a1c
## Pacuta_ATHC_TP1_1207      #e31a1c
## Pacuta_ATHC_TP11_2668     #e31a1c
## Pacuta_ATHC_TP11_2879     #e31a1c
## Pacuta_ATHC_TP1_2977      #e31a1c
## Pacuta_ATHC_TP3_1219      #e31a1c
## Pacuta_ATHC_TP4_2993      #e31a1c
## Pacuta_ATHC_TP5_1296      #e31a1c
## Pacuta_ATHC_TP6_2999      #e31a1c
## Pacuta_ATHC_TP7_1281      #e31a1c
## Pacuta_ATHC_TP8_2861      #e31a1c
## Pacuta_ATHC_TP9_2979      #e31a1c
## Pacuta_HTAC_TP10_1225     #e31a1c
## Pacuta_HTAC_TP3_2026      #e31a1c
## Pacuta_HTAC_TP5_1303      #e31a1c
## Pacuta_HTAC_TP5_1571      #e31a1c
## Pacuta_HTAC_TP6_1330      #e31a1c
## Pacuta_HTAC_TP7_1487      #e31a1c
## Pacuta_HTAC_TP8_1329      #e31a1c
## Pacuta_HTHC_TP11_1416     #e31a1c
## Pacuta_HTHC_TP1_1676      #e31a1c
## Pacuta_HTHC_TP3_1418      #e31a1c
## Pacuta_HTHC_TP3_2527      #e31a1c
## Pacuta_HTHC_TP4_1169      #e31a1c
## Pacuta_HTHC_TP5_1168      #e31a1c
## Pacuta_ATAC_TP5_1059      #b15928
## Pacuta_HTHC_TP1_2210      #b15928
## Pacuta_HTAC_TP1_1653      #000000
## Pacuta_HTHC_TP6_1721      #000000
## Pacuta_ATAC_TP7_1445      #808080
## Pacuta_ATHC_TP10_1205     #808080
## Pacuta_ATHC_TP10_2550     #808080
## Pacuta_ATHC_TP6_1254      #808080
## Pacuta_HTHC_TP11_2185     #808080
## Pacuta_HTHC_TP5_1415      #808080
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
