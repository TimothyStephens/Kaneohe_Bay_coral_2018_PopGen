---
title: "Plot `ANGSD` results for *P. acuta* RNA-seq samples"
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
samples.info <- read.table("../../samples_Pacuta.annotations.txt", header=T, comment.char='')
rownames(samples.info) <- samples.info$sample

samples.order <- read.table("bam.filelist.labels", header=F)
samples.info <- samples.info[samples.order$V1,]
samples.info
```

```
##                                      sample species treatment timepoint plugid
## Pacuta_ATAC_TP10_1159 Pacuta_ATAC_TP10_1159  Pacuta      ATAC      TP10   1159
## Pacuta_ATAC_TP10_1559 Pacuta_ATAC_TP10_1559  Pacuta      ATAC      TP10   1559
## Pacuta_ATAC_TP10_1641 Pacuta_ATAC_TP10_1641  Pacuta      ATAC      TP10   1641
## Pacuta_ATAC_TP1_1043   Pacuta_ATAC_TP1_1043  Pacuta      ATAC       TP1   1043
## Pacuta_ATAC_TP11_1103 Pacuta_ATAC_TP11_1103  Pacuta      ATAC      TP11   1103
## Pacuta_ATAC_TP11_1777 Pacuta_ATAC_TP11_1777  Pacuta      ATAC      TP11   1777
## Pacuta_ATAC_TP11_2306 Pacuta_ATAC_TP11_2306  Pacuta      ATAC      TP11   2306
## Pacuta_ATAC_TP1_1775   Pacuta_ATAC_TP1_1775  Pacuta      ATAC       TP1   1775
## Pacuta_ATAC_TP1_2363   Pacuta_ATAC_TP1_2363  Pacuta      ATAC       TP1   2363
## Pacuta_ATAC_TP3_1041   Pacuta_ATAC_TP3_1041  Pacuta      ATAC       TP3   1041
## Pacuta_ATAC_TP3_1471   Pacuta_ATAC_TP3_1471  Pacuta      ATAC       TP3   1471
## Pacuta_ATAC_TP3_1637   Pacuta_ATAC_TP3_1637  Pacuta      ATAC       TP3   1637
## Pacuta_ATAC_TP4_1060   Pacuta_ATAC_TP4_1060  Pacuta      ATAC       TP4   1060
## Pacuta_ATAC_TP4_1762   Pacuta_ATAC_TP4_1762  Pacuta      ATAC       TP4   1762
## Pacuta_ATAC_TP4_2002   Pacuta_ATAC_TP4_2002  Pacuta      ATAC       TP4   2002
## Pacuta_ATAC_TP5_1059   Pacuta_ATAC_TP5_1059  Pacuta      ATAC       TP5   1059
## Pacuta_ATAC_TP5_1563   Pacuta_ATAC_TP5_1563  Pacuta      ATAC       TP5   1563
## Pacuta_ATAC_TP5_1757   Pacuta_ATAC_TP5_1757  Pacuta      ATAC       TP5   1757
## Pacuta_ATAC_TP6_1050   Pacuta_ATAC_TP6_1050  Pacuta      ATAC       TP6   1050
## Pacuta_ATAC_TP6_1468   Pacuta_ATAC_TP6_1468  Pacuta      ATAC       TP6   1468
## Pacuta_ATAC_TP6_1542   Pacuta_ATAC_TP6_1542  Pacuta      ATAC       TP6   1542
## Pacuta_ATAC_TP7_1047   Pacuta_ATAC_TP7_1047  Pacuta      ATAC       TP7   1047
## Pacuta_ATAC_TP7_1445   Pacuta_ATAC_TP7_1445  Pacuta      ATAC       TP7   1445
## Pacuta_ATAC_TP7_2413   Pacuta_ATAC_TP7_2413  Pacuta      ATAC       TP7   2413
## Pacuta_ATAC_TP8_1051   Pacuta_ATAC_TP8_1051  Pacuta      ATAC       TP8   1051
## Pacuta_ATAC_TP8_1755   Pacuta_ATAC_TP8_1755  Pacuta      ATAC       TP8   1755
## Pacuta_ATAC_TP8_2012   Pacuta_ATAC_TP8_2012  Pacuta      ATAC       TP8   2012
## Pacuta_ATAC_TP9_1141   Pacuta_ATAC_TP9_1141  Pacuta      ATAC       TP9   1141
## Pacuta_ATAC_TP9_1594   Pacuta_ATAC_TP9_1594  Pacuta      ATAC       TP9   1594
## Pacuta_ATAC_TP9_2357   Pacuta_ATAC_TP9_2357  Pacuta      ATAC       TP9   2357
## Pacuta_ATHC_TP10_1205 Pacuta_ATHC_TP10_1205  Pacuta      ATHC      TP10   1205
## Pacuta_ATHC_TP10_2197 Pacuta_ATHC_TP10_2197  Pacuta      ATHC      TP10   2197
## Pacuta_ATHC_TP10_2550 Pacuta_ATHC_TP10_2550  Pacuta      ATHC      TP10   2550
## Pacuta_ATHC_TP11_1147 Pacuta_ATHC_TP11_1147  Pacuta      ATHC      TP11   1147
## Pacuta_ATHC_TP1_1207   Pacuta_ATHC_TP1_1207  Pacuta      ATHC       TP1   1207
## Pacuta_ATHC_TP11_2668 Pacuta_ATHC_TP11_2668  Pacuta      ATHC      TP11   2668
## Pacuta_ATHC_TP11_2879 Pacuta_ATHC_TP11_2879  Pacuta      ATHC      TP11   2879
## Pacuta_ATHC_TP1_2743   Pacuta_ATHC_TP1_2743  Pacuta      ATHC       TP1   2743
## Pacuta_ATHC_TP1_2977   Pacuta_ATHC_TP1_2977  Pacuta      ATHC       TP1   2977
## Pacuta_ATHC_TP3_1219   Pacuta_ATHC_TP3_1219  Pacuta      ATHC       TP3   1219
## Pacuta_ATHC_TP3_2534   Pacuta_ATHC_TP3_2534  Pacuta      ATHC       TP3   2534
## Pacuta_ATHC_TP3_2750   Pacuta_ATHC_TP3_2750  Pacuta      ATHC       TP3   2750
## Pacuta_ATHC_TP4_1220   Pacuta_ATHC_TP4_1220  Pacuta      ATHC       TP4   1220
## Pacuta_ATHC_TP4_2733   Pacuta_ATHC_TP4_2733  Pacuta      ATHC       TP4   2733
## Pacuta_ATHC_TP4_2993   Pacuta_ATHC_TP4_2993  Pacuta      ATHC       TP4   2993
## Pacuta_ATHC_TP5_1296   Pacuta_ATHC_TP5_1296  Pacuta      ATHC       TP5   1296
## Pacuta_ATHC_TP5_2212   Pacuta_ATHC_TP5_2212  Pacuta      ATHC       TP5   2212
## Pacuta_ATHC_TP5_2877   Pacuta_ATHC_TP5_2877  Pacuta      ATHC       TP5   2877
## Pacuta_ATHC_TP6_1254   Pacuta_ATHC_TP6_1254  Pacuta      ATHC       TP6   1254
## Pacuta_ATHC_TP6_2870   Pacuta_ATHC_TP6_2870  Pacuta      ATHC       TP6   2870
## Pacuta_ATHC_TP6_2999   Pacuta_ATHC_TP6_2999  Pacuta      ATHC       TP6   2999
## Pacuta_ATHC_TP7_1281   Pacuta_ATHC_TP7_1281  Pacuta      ATHC       TP7   1281
## Pacuta_ATHC_TP7_2409   Pacuta_ATHC_TP7_2409  Pacuta      ATHC       TP7   2409
## Pacuta_ATHC_TP7_2878   Pacuta_ATHC_TP7_2878  Pacuta      ATHC       TP7   2878
## Pacuta_ATHC_TP8_1459   Pacuta_ATHC_TP8_1459  Pacuta      ATHC       TP8   1459
## Pacuta_ATHC_TP8_2564   Pacuta_ATHC_TP8_2564  Pacuta      ATHC       TP8   2564
## Pacuta_ATHC_TP8_2861   Pacuta_ATHC_TP8_2861  Pacuta      ATHC       TP8   2861
## Pacuta_ATHC_TP9_1451   Pacuta_ATHC_TP9_1451  Pacuta      ATHC       TP9   1451
## Pacuta_ATHC_TP9_2873   Pacuta_ATHC_TP9_2873  Pacuta      ATHC       TP9   2873
## Pacuta_ATHC_TP9_2979   Pacuta_ATHC_TP9_2979  Pacuta      ATHC       TP9   2979
## Pacuta_HTAC_TP10_1225 Pacuta_HTAC_TP10_1225  Pacuta      HTAC      TP10   1225
## Pacuta_HTAC_TP10_1536 Pacuta_HTAC_TP10_1536  Pacuta      HTAC      TP10   1536
## Pacuta_HTAC_TP10_2064 Pacuta_HTAC_TP10_2064  Pacuta      HTAC      TP10   2064
## Pacuta_HTAC_TP11_1582 Pacuta_HTAC_TP11_1582  Pacuta      HTAC      TP11   1582
## Pacuta_HTAC_TP11_1596 Pacuta_HTAC_TP11_1596  Pacuta      HTAC      TP11   1596
## Pacuta_HTAC_TP11_1647 Pacuta_HTAC_TP11_1647  Pacuta      HTAC      TP11   1647
## Pacuta_HTAC_TP1_1653   Pacuta_HTAC_TP1_1653  Pacuta      HTAC       TP1   1653
## Pacuta_HTAC_TP1_2005   Pacuta_HTAC_TP1_2005  Pacuta      HTAC       TP1   2005
## Pacuta_HTAC_TP1_2414   Pacuta_HTAC_TP1_2414  Pacuta      HTAC       TP1   2414
## Pacuta_HTAC_TP3_1617   Pacuta_HTAC_TP3_1617  Pacuta      HTAC       TP3   1617
## Pacuta_HTAC_TP3_1642   Pacuta_HTAC_TP3_1642  Pacuta      HTAC       TP3   1642
## Pacuta_HTAC_TP3_2026   Pacuta_HTAC_TP3_2026  Pacuta      HTAC       TP3   2026
## Pacuta_HTAC_TP4_1581   Pacuta_HTAC_TP4_1581  Pacuta      HTAC       TP4   1581
## Pacuta_HTAC_TP4_1701   Pacuta_HTAC_TP4_1701  Pacuta      HTAC       TP4   1701
## Pacuta_HTAC_TP4_1767   Pacuta_HTAC_TP4_1767  Pacuta      HTAC       TP4   1767
## Pacuta_HTAC_TP5_1303   Pacuta_HTAC_TP5_1303  Pacuta      HTAC       TP5   1303
## Pacuta_HTAC_TP5_1571   Pacuta_HTAC_TP5_1571  Pacuta      HTAC       TP5   1571
## Pacuta_HTAC_TP5_1707   Pacuta_HTAC_TP5_1707  Pacuta      HTAC       TP5   1707
## Pacuta_HTAC_TP6_1330   Pacuta_HTAC_TP6_1330  Pacuta      HTAC       TP6   1330
## Pacuta_HTAC_TP6_1466   Pacuta_HTAC_TP6_1466  Pacuta      HTAC       TP6   1466
## Pacuta_HTAC_TP6_1744   Pacuta_HTAC_TP6_1744  Pacuta      HTAC       TP6   1744
## Pacuta_HTAC_TP7_1487   Pacuta_HTAC_TP7_1487  Pacuta      HTAC       TP7   1487
## Pacuta_HTAC_TP7_1728   Pacuta_HTAC_TP7_1728  Pacuta      HTAC       TP7   1728
## Pacuta_HTAC_TP7_2072   Pacuta_HTAC_TP7_2072  Pacuta      HTAC       TP7   2072
## Pacuta_HTAC_TP8_1329   Pacuta_HTAC_TP8_1329  Pacuta      HTAC       TP8   1329
## Pacuta_HTAC_TP8_1765   Pacuta_HTAC_TP8_1765  Pacuta      HTAC       TP8   1765
## Pacuta_HTAC_TP8_2513   Pacuta_HTAC_TP8_2513  Pacuta      HTAC       TP8   2513
## Pacuta_HTAC_TP9_1302   Pacuta_HTAC_TP9_1302  Pacuta      HTAC       TP9   1302
## Pacuta_HTAC_TP9_1486   Pacuta_HTAC_TP9_1486  Pacuta      HTAC       TP9   1486
## Pacuta_HTAC_TP9_1696   Pacuta_HTAC_TP9_1696  Pacuta      HTAC       TP9   1696
## Pacuta_HTHC_TP10_1238 Pacuta_HTHC_TP10_1238  Pacuta      HTHC      TP10   1238
## Pacuta_HTHC_TP10_1732 Pacuta_HTHC_TP10_1732  Pacuta      HTHC      TP10   1732
## Pacuta_HTHC_TP10_2300 Pacuta_HTHC_TP10_2300  Pacuta      HTHC      TP10   2300
## Pacuta_HTHC_TP11_1416 Pacuta_HTHC_TP11_1416  Pacuta      HTHC      TP11   1416
## Pacuta_HTHC_TP11_2185 Pacuta_HTHC_TP11_2185  Pacuta      HTHC      TP11   2185
## Pacuta_HTHC_TP1_1239   Pacuta_HTHC_TP1_1239  Pacuta      HTHC       TP1   1239
## Pacuta_HTHC_TP1_1676   Pacuta_HTHC_TP1_1676  Pacuta      HTHC       TP1   1676
## Pacuta_HTHC_TP1_2210   Pacuta_HTHC_TP1_2210  Pacuta      HTHC       TP1   2210
## Pacuta_HTHC_TP3_1227   Pacuta_HTHC_TP3_1227  Pacuta      HTHC       TP3   1227
## Pacuta_HTHC_TP3_1418   Pacuta_HTHC_TP3_1418  Pacuta      HTHC       TP3   1418
## Pacuta_HTHC_TP3_2527   Pacuta_HTHC_TP3_2527  Pacuta      HTHC       TP3   2527
## Pacuta_HTHC_TP4_1169   Pacuta_HTHC_TP4_1169  Pacuta      HTHC       TP4   1169
## Pacuta_HTHC_TP4_1343   Pacuta_HTHC_TP4_1343  Pacuta      HTHC       TP4   1343
## Pacuta_HTHC_TP4_2195   Pacuta_HTHC_TP4_2195  Pacuta      HTHC       TP4   2195
## Pacuta_HTHC_TP5_1168   Pacuta_HTHC_TP5_1168  Pacuta      HTHC       TP5   1168
## Pacuta_HTHC_TP5_1415   Pacuta_HTHC_TP5_1415  Pacuta      HTHC       TP5   1415
## Pacuta_HTHC_TP5_2087   Pacuta_HTHC_TP5_2087  Pacuta      HTHC       TP5   2087
## Pacuta_HTHC_TP6_1138   Pacuta_HTHC_TP6_1138  Pacuta      HTHC       TP6   1138
## Pacuta_HTHC_TP6_1595   Pacuta_HTHC_TP6_1595  Pacuta      HTHC       TP6   1595
## Pacuta_HTHC_TP6_1721   Pacuta_HTHC_TP6_1721  Pacuta      HTHC       TP6   1721
## Pacuta_HTHC_TP7_1090   Pacuta_HTHC_TP7_1090  Pacuta      HTHC       TP7   1090
## Pacuta_HTHC_TP7_1427   Pacuta_HTHC_TP7_1427  Pacuta      HTHC       TP7   1427
## Pacuta_HTHC_TP7_1820   Pacuta_HTHC_TP7_1820  Pacuta      HTHC       TP7   1820
## Pacuta_HTHC_TP8_1184   Pacuta_HTHC_TP8_1184  Pacuta      HTHC       TP8   1184
## Pacuta_HTHC_TP8_1709   Pacuta_HTHC_TP8_1709  Pacuta      HTHC       TP8   1709
## Pacuta_HTHC_TP8_2304   Pacuta_HTHC_TP8_2304  Pacuta      HTHC       TP8   2304
## Pacuta_HTHC_TP9_1131   Pacuta_HTHC_TP9_1131  Pacuta      HTHC       TP9   1131
## Pacuta_HTHC_TP9_2202   Pacuta_HTHC_TP9_2202  Pacuta      HTHC       TP9   2202
## Pacuta_HTHC_TP9_2305   Pacuta_HTHC_TP9_2305  Pacuta      HTHC       TP9   2305
##                                  reef reef_color ploidy ploidy_color   group
## Pacuta_ATAC_TP10_1159      Reef.42.43    #542788      2      #1b9e77  Group5
## Pacuta_ATAC_TP10_1559 Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_ATAC_TP10_1641            HIMB    #b35806      3      #d95f02  Group3
## Pacuta_ATAC_TP1_1043  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_ATAC_TP11_1103      Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_ATAC_TP11_1777            HIMB    #b35806      3      #d95f02  Group1
## Pacuta_ATAC_TP11_2306         Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATAC_TP1_1775       Reef.42.43    #542788      2      #1b9e77  Group5
## Pacuta_ATAC_TP1_2363          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_ATAC_TP3_1041       Reef.11.13    #fee0b6      3      #d95f02  Group2
## Pacuta_ATAC_TP3_1471       Reef.35.36    #998ec3      3      #d95f02  Group3
## Pacuta_ATAC_TP3_1637       Reef.11.13    #fee0b6      3      #d95f02  Group2
## Pacuta_ATAC_TP4_1060          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_ATAC_TP4_1762       Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_ATAC_TP4_2002       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_ATAC_TP5_1059       Reef.11.13    #fee0b6      2      #1b9e77  Group7
## Pacuta_ATAC_TP5_1563  Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_ATAC_TP5_1757       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_ATAC_TP6_1050          Reef.18    #d8daeb      2      #1b9e77  Group6
## Pacuta_ATAC_TP6_1468       Reef.35.36    #998ec3      2      #1b9e77  Group5
## Pacuta_ATAC_TP6_1542  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_ATAC_TP7_1047  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_ATAC_TP7_1445       Reef.11.13    #fee0b6      2      #1b9e77 Ungroup
## Pacuta_ATAC_TP7_2413          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATAC_TP8_1051       Reef.11.13    #fee0b6      3      #d95f02  Group1
## Pacuta_ATAC_TP8_1755       Reef.42.43    #542788      2      #1b9e77  Group5
## Pacuta_ATAC_TP8_2012       Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_ATAC_TP9_1141       Reef.42.43    #542788      2      #1b9e77  Group5
## Pacuta_ATAC_TP9_1594       Reef.35.36    #998ec3      3      #d95f02  Group3
## Pacuta_ATAC_TP9_2357          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATHC_TP10_1205      Reef.35.36    #998ec3      2      #1b9e77 Ungroup
## Pacuta_ATHC_TP10_2197      Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_ATHC_TP10_2550         Reef.18    #d8daeb      2      #1b9e77 Ungroup
## Pacuta_ATHC_TP11_1147            HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP1_1207       Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_ATHC_TP11_2668 Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_ATHC_TP11_2879            HIMB    #b35806      2      #1b9e77  Group6
## Pacuta_ATHC_TP1_2743  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_ATHC_TP1_2977       Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_ATHC_TP3_1219       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_ATHC_TP3_2534          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATHC_TP3_2750       Reef.42.43    #542788      3      #d95f02  Group4
## Pacuta_ATHC_TP4_1220  Lilipuna.Fringe    #f1a340      3      #d95f02  Group2
## Pacuta_ATHC_TP4_2733       Reef.42.43    #542788      3      #d95f02  Group1
## Pacuta_ATHC_TP4_2993       Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_ATHC_TP5_1296       Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_ATHC_TP5_2212       Reef.35.36    #998ec3      2      #1b9e77  Group5
## Pacuta_ATHC_TP5_2877             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP6_1254       Reef.42.43    #542788      3      #d95f02 Ungroup
## Pacuta_ATHC_TP6_2870             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP6_2999  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_ATHC_TP7_1281       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_ATHC_TP7_2409  Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_ATHC_TP7_2878             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP8_1459          Reef.18    #d8daeb      3      #d95f02  Group1
## Pacuta_ATHC_TP8_2564          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_ATHC_TP8_2861             HIMB    #b35806      2      #1b9e77  Group6
## Pacuta_ATHC_TP9_1451          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_ATHC_TP9_2873             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_ATHC_TP9_2979       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTAC_TP10_1225            HIMB    #b35806      2      #1b9e77  Group6
## Pacuta_HTAC_TP10_1536      Reef.35.36    #998ec3      3      #d95f02  Group3
## Pacuta_HTAC_TP10_2064            HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTAC_TP11_1582 Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_HTAC_TP11_1596      Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_HTAC_TP11_1647            HIMB    #b35806      3      #d95f02  Group3
## Pacuta_HTAC_TP1_1653             HIMB    #b35806      2      #1b9e77  Group8
## Pacuta_HTAC_TP1_2005          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_HTAC_TP1_2414          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTAC_TP3_1617       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTAC_TP3_1642             HIMB    #b35806      3      #d95f02  Group1
## Pacuta_HTAC_TP3_2026       Reef.42.43    #542788      2      #1b9e77  Group6
## Pacuta_HTAC_TP4_1581  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_HTAC_TP4_1701          Reef.18    #d8daeb      3      #d95f02  Group3
## Pacuta_HTAC_TP4_1767       Reef.42.43    #542788      3      #d95f02  Group1
## Pacuta_HTAC_TP5_1303       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTAC_TP5_1571       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTAC_TP5_1707       Reef.11.13    #fee0b6      3      #d95f02  Group2
## Pacuta_HTAC_TP6_1330       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTAC_TP6_1466  Lilipuna.Fringe    #f1a340      3      #d95f02  Group2
## Pacuta_HTAC_TP6_1744       Reef.35.36    #998ec3      2      #1b9e77  Group5
## Pacuta_HTAC_TP7_1487  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_HTAC_TP7_1728          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTAC_TP7_2072             HIMB    #b35806      3      #d95f02  Group1
## Pacuta_HTAC_TP8_1329  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_HTAC_TP8_1765       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTAC_TP8_2513          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTAC_TP9_1302       Reef.35.36    #998ec3      2      #1b9e77  Group5
## Pacuta_HTAC_TP9_1486       Reef.11.13    #fee0b6      2      #1b9e77  Group5
## Pacuta_HTAC_TP9_1696  Lilipuna.Fringe    #f1a340      3      #d95f02  Group2
## Pacuta_HTHC_TP10_1238      Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTHC_TP10_1732 Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_HTHC_TP10_2300         Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTHC_TP11_1416      Reef.11.13    #fee0b6      2      #1b9e77  Group6
## Pacuta_HTHC_TP11_2185            HIMB    #b35806      3      #d95f02 Ungroup
## Pacuta_HTHC_TP1_1239       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTHC_TP1_1676       Reef.42.43    #542788      2      #1b9e77  Group6
## Pacuta_HTHC_TP1_2210       Reef.42.43    #542788      2      #1b9e77  Group7
## Pacuta_HTHC_TP3_1227       Reef.42.43    #542788      3      #d95f02  Group2
## Pacuta_HTHC_TP3_1418  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group6
## Pacuta_HTHC_TP3_2527          Reef.18    #d8daeb      2      #1b9e77  Group6
## Pacuta_HTHC_TP4_1169       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTHC_TP4_1343  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_HTHC_TP4_2195       Reef.42.43    #542788      3      #d95f02  Group3
## Pacuta_HTHC_TP5_1168       Reef.35.36    #998ec3      2      #1b9e77  Group6
## Pacuta_HTHC_TP5_1415       Reef.11.13    #fee0b6      2      #1b9e77 Ungroup
## Pacuta_HTHC_TP5_2087             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP6_1138             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP6_1595             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP6_1721  Lilipuna.Fringe    #f1a340      2      #1b9e77  Group8
## Pacuta_HTHC_TP7_1090  Lilipuna.Fringe    #f1a340      3      #d95f02  Group3
## Pacuta_HTHC_TP7_1427  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_HTHC_TP7_1820       Reef.11.13    #fee0b6      3      #d95f02  Group3
## Pacuta_HTHC_TP8_1184  Lilipuna.Fringe    #f1a340      3      #d95f02  Group4
## Pacuta_HTHC_TP8_1709  Lilipuna.Fringe    #f1a340      3      #d95f02  Group2
## Pacuta_HTHC_TP8_2304          Reef.18    #d8daeb      3      #d95f02  Group2
## Pacuta_HTHC_TP9_1131             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP9_2202             HIMB    #b35806      3      #d95f02  Group2
## Pacuta_HTHC_TP9_2305          Reef.18    #d8daeb      3      #d95f02  Group2
##                       group_color
## Pacuta_ATAC_TP10_1159     #ff7f00
## Pacuta_ATAC_TP10_1559     #c51b7d
## Pacuta_ATAC_TP10_1641     #c51b7d
## Pacuta_ATAC_TP1_1043      #6a3d9a
## Pacuta_ATAC_TP11_1103     #33a02c
## Pacuta_ATAC_TP11_1777     #1f78b4
## Pacuta_ATAC_TP11_2306     #33a02c
## Pacuta_ATAC_TP1_1775      #ff7f00
## Pacuta_ATAC_TP1_2363      #c51b7d
## Pacuta_ATAC_TP3_1041      #33a02c
## Pacuta_ATAC_TP3_1471      #c51b7d
## Pacuta_ATAC_TP3_1637      #33a02c
## Pacuta_ATAC_TP4_1060      #c51b7d
## Pacuta_ATAC_TP4_1762      #33a02c
## Pacuta_ATAC_TP4_2002      #c51b7d
## Pacuta_ATAC_TP5_1059      #b15928
## Pacuta_ATAC_TP5_1563      #c51b7d
## Pacuta_ATAC_TP5_1757      #c51b7d
## Pacuta_ATAC_TP6_1050      #e31a1c
## Pacuta_ATAC_TP6_1468      #ff7f00
## Pacuta_ATAC_TP6_1542      #6a3d9a
## Pacuta_ATAC_TP7_1047      #e31a1c
## Pacuta_ATAC_TP7_1445      #808080
## Pacuta_ATAC_TP7_2413      #33a02c
## Pacuta_ATAC_TP8_1051      #1f78b4
## Pacuta_ATAC_TP8_1755      #ff7f00
## Pacuta_ATAC_TP8_2012      #33a02c
## Pacuta_ATAC_TP9_1141      #ff7f00
## Pacuta_ATAC_TP9_1594      #c51b7d
## Pacuta_ATAC_TP9_2357      #33a02c
## Pacuta_ATHC_TP10_1205     #808080
## Pacuta_ATHC_TP10_2197     #e31a1c
## Pacuta_ATHC_TP10_2550     #808080
## Pacuta_ATHC_TP11_1147     #33a02c
## Pacuta_ATHC_TP1_1207      #e31a1c
## Pacuta_ATHC_TP11_2668     #e31a1c
## Pacuta_ATHC_TP11_2879     #e31a1c
## Pacuta_ATHC_TP1_2743      #6a3d9a
## Pacuta_ATHC_TP1_2977      #e31a1c
## Pacuta_ATHC_TP3_1219      #e31a1c
## Pacuta_ATHC_TP3_2534      #33a02c
## Pacuta_ATHC_TP3_2750      #6a3d9a
## Pacuta_ATHC_TP4_1220      #33a02c
## Pacuta_ATHC_TP4_2733      #1f78b4
## Pacuta_ATHC_TP4_2993      #e31a1c
## Pacuta_ATHC_TP5_1296      #e31a1c
## Pacuta_ATHC_TP5_2212      #ff7f00
## Pacuta_ATHC_TP5_2877      #33a02c
## Pacuta_ATHC_TP6_1254      #808080
## Pacuta_ATHC_TP6_2870      #33a02c
## Pacuta_ATHC_TP6_2999      #e31a1c
## Pacuta_ATHC_TP7_1281      #e31a1c
## Pacuta_ATHC_TP7_2409      #c51b7d
## Pacuta_ATHC_TP7_2878      #33a02c
## Pacuta_ATHC_TP8_1459      #1f78b4
## Pacuta_ATHC_TP8_2564      #c51b7d
## Pacuta_ATHC_TP8_2861      #e31a1c
## Pacuta_ATHC_TP9_1451      #33a02c
## Pacuta_ATHC_TP9_2873      #33a02c
## Pacuta_ATHC_TP9_2979      #e31a1c
## Pacuta_HTAC_TP10_1225     #e31a1c
## Pacuta_HTAC_TP10_1536     #c51b7d
## Pacuta_HTAC_TP10_2064     #33a02c
## Pacuta_HTAC_TP11_1582     #c51b7d
## Pacuta_HTAC_TP11_1596     #33a02c
## Pacuta_HTAC_TP11_1647     #c51b7d
## Pacuta_HTAC_TP1_1653      #000000
## Pacuta_HTAC_TP1_2005      #c51b7d
## Pacuta_HTAC_TP1_2414      #33a02c
## Pacuta_HTAC_TP3_1617      #c51b7d
## Pacuta_HTAC_TP3_1642      #1f78b4
## Pacuta_HTAC_TP3_2026      #e31a1c
## Pacuta_HTAC_TP4_1581      #6a3d9a
## Pacuta_HTAC_TP4_1701      #c51b7d
## Pacuta_HTAC_TP4_1767      #1f78b4
## Pacuta_HTAC_TP5_1303      #e31a1c
## Pacuta_HTAC_TP5_1571      #e31a1c
## Pacuta_HTAC_TP5_1707      #33a02c
## Pacuta_HTAC_TP6_1330      #e31a1c
## Pacuta_HTAC_TP6_1466      #33a02c
## Pacuta_HTAC_TP6_1744      #ff7f00
## Pacuta_HTAC_TP7_1487      #e31a1c
## Pacuta_HTAC_TP7_1728      #33a02c
## Pacuta_HTAC_TP7_2072      #1f78b4
## Pacuta_HTAC_TP8_1329      #e31a1c
## Pacuta_HTAC_TP8_1765      #c51b7d
## Pacuta_HTAC_TP8_2513      #33a02c
## Pacuta_HTAC_TP9_1302      #ff7f00
## Pacuta_HTAC_TP9_1486      #ff7f00
## Pacuta_HTAC_TP9_1696      #33a02c
## Pacuta_HTHC_TP10_1238     #c51b7d
## Pacuta_HTHC_TP10_1732     #c51b7d
## Pacuta_HTHC_TP10_2300     #33a02c
## Pacuta_HTHC_TP11_1416     #e31a1c
## Pacuta_HTHC_TP11_2185     #808080
## Pacuta_HTHC_TP1_1239      #c51b7d
## Pacuta_HTHC_TP1_1676      #e31a1c
## Pacuta_HTHC_TP1_2210      #b15928
## Pacuta_HTHC_TP3_1227      #33a02c
## Pacuta_HTHC_TP3_1418      #e31a1c
## Pacuta_HTHC_TP3_2527      #e31a1c
## Pacuta_HTHC_TP4_1169      #e31a1c
## Pacuta_HTHC_TP4_1343      #6a3d9a
## Pacuta_HTHC_TP4_2195      #c51b7d
## Pacuta_HTHC_TP5_1168      #e31a1c
## Pacuta_HTHC_TP5_1415      #808080
## Pacuta_HTHC_TP5_2087      #33a02c
## Pacuta_HTHC_TP6_1138      #33a02c
## Pacuta_HTHC_TP6_1595      #33a02c
## Pacuta_HTHC_TP6_1721      #000000
## Pacuta_HTHC_TP7_1090      #c51b7d
## Pacuta_HTHC_TP7_1427      #6a3d9a
## Pacuta_HTHC_TP7_1820      #c51b7d
## Pacuta_HTHC_TP8_1184      #6a3d9a
## Pacuta_HTHC_TP8_1709      #33a02c
## Pacuta_HTHC_TP8_2304      #33a02c
## Pacuta_HTHC_TP9_1131      #33a02c
## Pacuta_HTHC_TP9_2202      #33a02c
## Pacuta_HTHC_TP9_2305      #33a02c
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
q  <- read.table("PCAngsd.angsd.beagle.gz.Admixture.admix.4.Q")
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
##      Reef.42.43 Lilipuna.Fringe            HIMB         Reef.18      Reef.11.13 
##       "#542788"       "#f1a340"       "#b35806"       "#d8daeb"       "#fee0b6" 
##      Reef.35.36               2               3          Group5          Group3 
##       "#998ec3"       "#1b9e77"       "#d95f02"       "#ff7f00"       "#c51b7d" 
##          Group4          Group2          Group1          Group7          Group6 
##       "#6a3d9a"       "#33a02c"       "#1f78b4"       "#b15928"       "#e31a1c" 
##         Ungroup          Group8 
##       "#808080"       "#000000"
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
<div id="htmlwidget-fdd64a30b71c6107f493" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-fdd64a30b71c6107f493">{"x":{"visdat":{"5568757bad2e":["function () ","plotlyVisDat"]},"cur_data":"5568757bad2e","attrs":{"5568757bad2e":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#1b9e77","#d95f02","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Group","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[0.0875214959861799,0.0918657943026994,0.0929869145666383,0.0914744114526632,0.0881821174072034,0.0918439423842616,0.0874552874920673],"y":[0.0226665552918805,0.0264128699694468,0.0262937398334947,0.0267015953045515,0.0231259994629201,0.026251770040791,0.0232234727966454],"text":["Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP8_1051","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP8_1459","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP7_2072"],"mode":"markers","marker":{"color":"rgba(31,120,180,1)","size":11,"line":{"color":"rgba(31,120,180,1)"}},"type":"scatter","name":"Group1","textfont":{"color":"rgba(31,120,180,1)"},"error_y":{"color":"rgba(31,120,180,1)"},"error_x":{"color":"rgba(31,120,180,1)"},"line":{"color":"rgba(31,120,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.107920995122979,0.110804472666828,0.110880300665118,0.109080503399851,0.112179601349592,0.109869777730301,0.109177949457003,0.108447960914815,0.110687785757126,0.112210625271056,0.107791391483647,0.110386599241113,0.109707532609245,0.11090246377058,0.109091291290452,0.108348331074995,0.110594808724283,0.111633397117612,0.111364757323081,0.112699290195878,0.109908627364511,0.111549988518253,0.108820267283964,0.114089298178451,0.110185128633608,0.110758804720614,0.111839239758775,0.110966425144904,0.11027457998866,0.112430465296924,0.110471286585083,0.108991821635593,0.110998305106978,0.108975062295148],"y":[0.0567998910969012,0.057750629069382,0.0576778462771068,0.0570986455215315,0.0582381134347802,0.0574973600181808,0.0567729847852835,0.0569742845168138,0.0578514419209095,0.058217650008454,0.0565117876505037,0.0578420634394934,0.0571928992041755,0.057889220045766,0.0573338293461615,0.056817597514372,0.0573950908648798,0.0584990197173927,0.0578828465329656,0.058634641210091,0.0574242509113888,0.0580888388275656,0.0568934359391014,0.0588755248711854,0.0571409594586328,0.0577069470787143,0.0580168801100127,0.0576738190616372,0.0574701331132224,0.0586519482508702,0.0577430414778805,0.0569493557085853,0.057650411083457,0.0570426095182013],"text":["Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"color":"rgba(51,160,44,1)","size":11,"line":{"color":"rgba(51,160,44,1)"}},"type":"scatter","name":"Group2","textfont":{"color":"rgba(51,160,44,1)"},"error_y":{"color":"rgba(51,160,44,1)"},"error_x":{"color":"rgba(51,160,44,1)"},"line":{"color":"rgba(51,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0115960959038586,-0.0114841121486476,-0.0116435617464128,-0.0117212716929476,-0.0115911396088818,-0.0113933054279697,-0.0112787559285863,-0.0115955815563203,-0.0117175827451194,-0.0114054840241217,-0.0112919780969947,-0.0115511861924536,-0.0114155914118386,-0.0113857047820858,-0.01154546447973,-0.0116834983777247,-0.0111767103177566,-0.0118584718015054,-0.0117271857075897,-0.0110312107587424,-0.0114348440339058,-0.0114736544432993,-0.0115220986253586,-0.0117284904503894],"y":[-0.16095947497115,-0.157392205594829,-0.159588697926679,-0.159375580372905,-0.160269956748766,-0.164624401067245,-0.165219331142488,-0.16066380093423,-0.160362611460552,-0.162453676727804,-0.16511544267949,-0.160869418668767,-0.164920775579217,-0.163911153357685,-0.160442686765566,-0.159876300987241,-0.165364285738174,-0.158315054202029,-0.161793665221103,-0.159219292363829,-0.16073958882592,-0.164309403699645,-0.161704084825625,-0.16053671851649],"text":["Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP9_1594","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP8_2564","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP8_1765","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1820"],"mode":"markers","marker":{"color":"rgba(197,27,125,1)","size":11,"line":{"color":"rgba(197,27,125,1)"}},"type":"scatter","name":"Group3","textfont":{"color":"rgba(197,27,125,1)"},"error_y":{"color":"rgba(197,27,125,1)"},"error_x":{"color":"rgba(197,27,125,1)"},"line":{"color":"rgba(197,27,125,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0447653757422054,-0.0451347264990633,-0.0447743791553109,-0.0461294747603586,-0.0467553605163089,-0.0455171429142934,-0.0451576738159804,-0.0443097809060383],"y":[-0.0788809490545534,-0.0798654616882042,-0.0792141795698008,-0.0819419988723811,-0.082526405845717,-0.0795037038555079,-0.0799773440500236,-0.0783595310675299],"text":["Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP6_1542","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP3_2750","Pacuta_HTAC_TP4_1581","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP8_1184"],"mode":"markers","marker":{"color":"rgba(106,61,154,1)","size":11,"line":{"color":"rgba(106,61,154,1)"}},"type":"scatter","name":"Group4","textfont":{"color":"rgba(106,61,154,1)"},"error_y":{"color":"rgba(106,61,154,1)"},"error_x":{"color":"rgba(106,61,154,1)"},"line":{"color":"rgba(106,61,154,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0350642831612321,-0.0352745681361068,-0.035532374509746,-0.0351705144493628,-0.0351080329625817,-0.0357273410950867,-0.0350124286398336,-0.035565407095257,-0.0350719013755009],"y":[-0.00924375280068004,-0.00951210782159871,-0.00967076567746238,-0.00928595953098906,-0.00925723355669288,-0.00965829330328564,-0.00929067169178636,-0.0093498681519262,-0.0094158997128757],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP9_1141","Pacuta_ATHC_TP5_2212","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486"],"mode":"markers","marker":{"color":"rgba(255,127,0,1)","size":11,"line":{"color":"rgba(255,127,0,1)"}},"type":"scatter","name":"Group5","textfont":{"color":"rgba(255,127,0,1)"},"error_y":{"color":"rgba(255,127,0,1)"},"error_x":{"color":"rgba(255,127,0,1)"},"line":{"color":"rgba(255,127,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.131463240984123,-0.132229601225089,-0.131573149055954,-0.134094118008552,-0.13009155077984,-0.132691533479052,-0.132506322347987,-0.13164974783941,-0.130786431844457,-0.131073351375603,-0.131117381631219,-0.132713020763172,-0.132315187652766,-0.1312336001026,-0.130560913880463,-0.131543570360092,-0.131622975378459,-0.133199436669728,-0.129990393206282,-0.131523693926475,-0.133993451928388,-0.132418775351293,-0.131354864724628,-0.131258206293756,-0.131604052969417,-0.133658402616273,-0.132096298266833],"y":[0.0838499210266198,0.084255228933984,0.0840450135246456,0.0853214381825325,0.08285245572717,0.0843219536663018,0.0845268933084605,0.084066134164872,0.0836072820955728,0.0837531249701068,0.0836737624910199,0.0843714794531121,0.0844596448671217,0.0838402762032223,0.0831701559046354,0.0837765113448058,0.083873164099392,0.0849260145068015,0.0832722341066018,0.08375426140716,0.0854545716104494,0.0842077758952999,0.0838418130679639,0.0836917901014634,0.0838444112186807,0.0849559785944652,0.0841935482877721],"text":["Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP7_1047","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP8_1329","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP5_1168"],"mode":"markers","marker":{"color":"rgba(227,26,28,1)","size":11,"line":{"color":"rgba(227,26,28,1)"}},"type":"scatter","name":"Group6","textfont":{"color":"rgba(227,26,28,1)"},"error_y":{"color":"rgba(227,26,28,1)"},"error_x":{"color":"rgba(227,26,28,1)"},"line":{"color":"rgba(227,26,28,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.056450322115197,-0.0676109871213616],"y":[0.0263383823584625,0.037502666774285],"text":["Pacuta_ATAC_TP5_1059","Pacuta_HTHC_TP1_2210"],"mode":"markers","marker":{"color":"rgba(177,89,40,1)","size":11,"line":{"color":"rgba(177,89,40,1)"}},"type":"scatter","name":"Group7","textfont":{"color":"rgba(177,89,40,1)"},"error_y":{"color":"rgba(177,89,40,1)"},"error_x":{"color":"rgba(177,89,40,1)"},"line":{"color":"rgba(177,89,40,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0595963613960706,-0.0571637824193076],"y":[0.0511589884446391,0.0531406316328668],"text":["Pacuta_HTAC_TP1_1653","Pacuta_HTHC_TP6_1721"],"mode":"markers","marker":{"color":"rgba(0,0,0,1)","size":11,"line":{"color":"rgba(0,0,0,1)"}},"type":"scatter","name":"Group8","textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0200528737793683,-0.0604329934597445,0.0191059156868235,0.0770253231394743,0.0151142490969404,0.0544808506136735],"y":[0.0189932885666385,0.0156779621189555,0.0651234846654442,0.022796986320039,-0.00842076397414266,-0.0139351238755829],"text":["Pacuta_ATAC_TP7_1445","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP6_1254","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP5_1415"],"mode":"markers","marker":{"color":"rgba(128,128,128,1)","size":11,"line":{"color":"rgba(128,128,128,1)"}},"type":"scatter","name":"Ungroup","textfont":{"color":"rgba(128,128,128,1)"},"error_y":{"color":"rgba(128,128,128,1)"},"error_x":{"color":"rgba(128,128,128,1)"},"line":{"color":"rgba(128,128,128,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-8e8bc2020aa454b883da" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-8e8bc2020aa454b883da">{"x":{"visdat":{"5568aa6e344":["function () ","plotlyVisDat"]},"cur_data":"5568aa6e344","attrs":{"5568aa6e344":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#1b9e77","#d95f02","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Reef","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0114841121486476,0.0875214959861799,0.110687785757126,-0.132691533479052,0.110386599241113,0.109707532609245,0.11090246377058,-0.132315187652766,0.108348331074995,-0.130560913880463,0.110594808724283,-0.0113857047820858,-0.0595963613960706,0.0881821174072034,0.0874552874920673,0.0151142490969404,0.111839239758775,0.110966425144904,0.11027457998866,0.108991821635593,0.110998305106978],"y":[-0.157392205594829,0.0226665552918805,0.0578514419209095,0.0843219536663018,0.0578420634394934,0.0571928992041755,0.057889220045766,0.0844596448671217,0.056817597514372,0.0831701559046354,0.0573950908648798,-0.163911153357685,0.0511589884446391,0.0231259994629201,0.0232234727966454,-0.00842076397414266,0.0580168801100127,0.0576738190616372,0.0574701331132224,0.0569493557085853,0.057650411083457],"text":["Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP11_1777","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_2873","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP7_2072","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202"],"mode":"markers","marker":{"color":"rgba(179,88,6,1)","size":11,"line":{"color":"rgba(179,88,6,1)"}},"type":"scatter","name":"HIMB","textfont":{"color":"rgba(179,88,6,1)"},"error_y":{"color":"rgba(179,88,6,1)"},"error_x":{"color":"rgba(179,88,6,1)"},"line":{"color":"rgba(179,88,6,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0115960959038586,-0.0447653757422054,-0.0112787559285863,-0.0451347264990633,-0.132229601225089,-0.13009155077984,-0.0447743791553109,0.107791391483647,-0.131117381631219,-0.0114054840241217,-0.0114155914118386,-0.0467553605163089,0.109908627364511,-0.131523693926475,-0.133993451928388,0.114089298178451,-0.0110312107587424,-0.131258206293756,-0.0455171429142934,-0.0571637824193076,-0.0115220986253586,-0.0451576738159804,-0.0443097809060383,0.112430465296924],"y":[-0.16095947497115,-0.0788809490545534,-0.165219331142488,-0.0798654616882042,0.084255228933984,0.08285245572717,-0.0792141795698008,0.0565117876505037,0.0836737624910199,-0.162453676727804,-0.164920775579217,-0.082526405845717,0.0574242509113888,0.08375426140716,0.0854545716104494,0.0588755248711854,-0.159219292363829,0.0836917901014634,-0.0795037038555079,0.0531406316328668,-0.161704084825625,-0.0799773440500236,-0.0783595310675299,0.0586519482508702],"text":["Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_2409","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709"],"mode":"markers","marker":{"color":"rgba(241,163,64,1)","size":11,"line":{"color":"rgba(241,163,64,1)"}},"type":"scatter","name":"Lilipuna.Fringe","textfont":{"color":"rgba(241,163,64,1)"},"error_y":{"color":"rgba(241,163,64,1)"},"error_x":{"color":"rgba(241,163,64,1)"},"line":{"color":"rgba(241,163,64,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.110880300665118,0.109080503399851,-0.056450322115197,-0.0200528737793683,0.0918657943026994,-0.134094118008552,-0.132506322347987,-0.130786431844457,-0.131073351375603,0.112699290195878,-0.0350719013755009,-0.132418775351293,0.0544808506136735,-0.0117284904503894],"y":[0.0576778462771068,0.0570986455215315,0.0263383823584625,0.0189932885666385,0.0264128699694468,0.0853214381825325,0.0845268933084605,0.0836072820955728,0.0837531249701068,0.058634641210091,-0.0094158997128757,0.0842077758952999,-0.0139351238755829,-0.16053671851649],"text":["Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP8_1051","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP9_1486","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP7_1820"],"mode":"markers","marker":{"color":"rgba(254,224,182,1)","size":11,"line":{"color":"rgba(254,224,182,1)"}},"type":"scatter","name":"Reef.11.13","textfont":{"color":"rgba(254,224,182,1)"},"error_y":{"color":"rgba(254,224,182,1)"},"error_x":{"color":"rgba(254,224,182,1)"},"line":{"color":"rgba(254,224,182,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.110804472666828,-0.0116435617464128,-0.0115911396088818,-0.131463240984123,0.109869777730301,0.108447960914815,0.0191059156868235,0.112210625271056,0.0914744114526632,-0.0112919780969947,0.109091291290452,-0.01154546447973,0.111364757323081,-0.0111767103177566,0.111549988518253,0.108820267283964,0.110185128633608,-0.131604052969417,0.110471286585083,0.108975062295148],"y":[0.057750629069382,-0.159588697926679,-0.160269956748766,0.0838499210266198,0.0574973600181808,0.0569742845168138,0.0651234846654442,0.058217650008454,0.0267015953045515,-0.16511544267949,0.0573338293461615,-0.160442686765566,0.0578828465329656,-0.165364285738174,0.0580888388275656,0.0568934359391014,0.0571409594586328,0.0838444112186807,0.0577430414778805,0.0570426095182013],"text":["Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP9_1451","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP8_2513","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"color":"rgba(216,218,235,1)","size":11,"line":{"color":"rgba(216,218,235,1)"}},"type":"scatter","name":"Reef.18","textfont":{"color":"rgba(216,218,235,1)"},"error_y":{"color":"rgba(216,218,235,1)"},"error_x":{"color":"rgba(216,218,235,1)"},"line":{"color":"rgba(216,218,235,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0117212716929476,-0.035532374509746,-0.0117175827451194,-0.0604329934597445,-0.131573149055954,-0.13164974783941,-0.0357273410950867,-0.132713020763172,-0.1312336001026,-0.0115511861924536,-0.131622975378459,-0.133199436669728,-0.129990393206282,-0.0350124286398336,-0.035565407095257,-0.133658402616273,-0.132096298266833],"y":[-0.159375580372905,-0.00967076567746238,-0.160362611460552,0.0156779621189555,0.0840450135246456,0.084066134164872,-0.00965829330328564,0.0843714794531121,0.0838402762032223,-0.160869418668767,0.083873164099392,0.0849260145068015,0.0832722341066018,-0.00929067169178636,-0.0093498681519262,0.0849559785944652,0.0841935482877721],"text":["Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP9_1594","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP9_1302","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP5_1168"],"mode":"markers","marker":{"color":"rgba(153,142,195,1)","size":11,"line":{"color":"rgba(153,142,195,1)"}},"type":"scatter","name":"Reef.35.36","textfont":{"color":"rgba(153,142,195,1)"},"error_y":{"color":"rgba(153,142,195,1)"},"error_x":{"color":"rgba(153,142,195,1)"},"line":{"color":"rgba(153,142,195,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0350642831612321,0.107920995122979,-0.0352745681361068,0.112179601349592,-0.0113933054279697,-0.0115955815563203,-0.0351705144493628,0.109177949457003,-0.0351080329625817,-0.0461294747603586,0.0929869145666383,0.0770253231394743,0.111633397117612,-0.0116834983777247,-0.131543570360092,0.0918439423842616,-0.0118584718015054,-0.0117271857075897,-0.0114348440339058,-0.131354864724628,-0.0676109871213616,0.110758804720614,-0.0114736544432993],"y":[-0.00924375280068004,0.0567998910969012,-0.00951210782159871,0.0582381134347802,-0.164624401067245,-0.16066380093423,-0.00928595953098906,0.0567729847852835,-0.00925723355669288,-0.0819419988723811,0.0262937398334947,0.022796986320039,0.0584990197173927,-0.159876300987241,0.0837765113448058,0.026251770040791,-0.158315054202029,-0.161793665221103,-0.16073958882592,0.0838418130679639,0.037502666774285,0.0577069470787143,-0.164309403699645],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP6_1254","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP8_1765","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP4_2195"],"mode":"markers","marker":{"color":"rgba(84,39,136,1)","size":11,"line":{"color":"rgba(84,39,136,1)"}},"type":"scatter","name":"Reef.42.43","textfont":{"color":"rgba(84,39,136,1)"},"error_y":{"color":"rgba(84,39,136,1)"},"error_x":{"color":"rgba(84,39,136,1)"},"line":{"color":"rgba(84,39,136,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-f369a638b722aa269c50" style="width:800px;height:800px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-f369a638b722aa269c50">{"x":{"visdat":{"556852928ad2":["function () ","plotlyVisDat"]},"cur_data":"556852928ad2","attrs":{"556852928ad2":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#1b9e77","#d95f02","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"width":800,"height":800,"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Ploidy","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0350642831612321,-0.0115960959038586,-0.0114841121486476,-0.0447653757422054,0.107920995122979,0.0875214959861799,0.110804472666828,-0.0352745681361068,-0.0116435617464128,0.110880300665118,-0.0117212716929476,0.109080503399851,-0.0115911396088818,0.112179601349592,-0.0113933054279697,-0.056450322115197,-0.0112787559285863,-0.0115955815563203,-0.131463240984123,-0.035532374509746,-0.0451347264990633,-0.132229601225089,-0.0200528737793683,0.109869777730301,0.0918657943026994,-0.0351705144493628,0.109177949457003,-0.0351080329625817,-0.0117175827451194,0.108447960914815,-0.0604329934597445,-0.131573149055954,0.0191059156868235,0.110687785757126,-0.134094118008552,-0.13009155077984,-0.132691533479052,-0.0447743791553109,-0.132506322347987,-0.13164974783941,0.112210625271056,-0.0461294747603586,0.107791391483647,0.0929869145666383,-0.130786431844457,-0.131073351375603,-0.0357273410950867,0.110386599241113,0.0770253231394743,0.109707532609245,-0.131117381631219,-0.132713020763172,-0.0114054840241217,0.11090246377058,0.0914744114526632,-0.0112919780969947,-0.132315187652766,0.109091291290452,0.108348331074995,-0.1312336001026,-0.130560913880463,-0.0115511861924536,0.110594808724283,-0.0114155914118386,0.111633397117612,-0.0113857047820858,-0.0595963613960706,-0.01154546447973,0.111364757323081,-0.0116834983777247,0.0881821174072034,-0.131543570360092,-0.0467553605163089,-0.0111767103177566,0.0918439423842616,-0.131622975378459,-0.133199436669728,0.112699290195878,-0.129990393206282,0.109908627364511,-0.0350124286398336,-0.131523693926475,0.111549988518253,0.0874552874920673,-0.133993451928388,-0.0118584718015054,0.108820267283964,-0.035565407095257,-0.0350719013755009,0.114089298178451,-0.0117271857075897,-0.0110312107587424,0.110185128633608,-0.132418775351293,0.0151142490969404,-0.0114348440339058,-0.131354864724628,-0.0676109871213616,0.110758804720614,-0.131258206293756,-0.131604052969417,-0.133658402616273,-0.0455171429142934,-0.0114736544432993,-0.132096298266833,0.0544808506136735,0.111839239758775,0.110966425144904,0.11027457998866,-0.0571637824193076,-0.0115220986253586,-0.0451576738159804,-0.0117284904503894,-0.0443097809060383,0.112430465296924,0.110471286585083,0.108991821635593,0.110998305106978,0.108975062295148],"y":[-0.00924375280068004,-0.16095947497115,-0.157392205594829,-0.0788809490545534,0.0567998910969012,0.0226665552918805,0.057750629069382,-0.00951210782159871,-0.159588697926679,0.0576778462771068,-0.159375580372905,0.0570986455215315,-0.160269956748766,0.0582381134347802,-0.164624401067245,0.0263383823584625,-0.165219331142488,-0.16066380093423,0.0838499210266198,-0.00967076567746238,-0.0798654616882042,0.084255228933984,0.0189932885666385,0.0574973600181808,0.0264128699694468,-0.00928595953098906,0.0567729847852835,-0.00925723355669288,-0.160362611460552,0.0569742845168138,0.0156779621189555,0.0840450135246456,0.0651234846654442,0.0578514419209095,0.0853214381825325,0.08285245572717,0.0843219536663018,-0.0792141795698008,0.0845268933084605,0.084066134164872,0.058217650008454,-0.0819419988723811,0.0565117876505037,0.0262937398334947,0.0836072820955728,0.0837531249701068,-0.00965829330328564,0.0578420634394934,0.022796986320039,0.0571928992041755,0.0836737624910199,0.0843714794531121,-0.162453676727804,0.057889220045766,0.0267015953045515,-0.16511544267949,0.0844596448671217,0.0573338293461615,0.056817597514372,0.0838402762032223,0.0831701559046354,-0.160869418668767,0.0573950908648798,-0.164920775579217,0.0584990197173927,-0.163911153357685,0.0511589884446391,-0.160442686765566,0.0578828465329656,-0.159876300987241,0.0231259994629201,0.0837765113448058,-0.082526405845717,-0.165364285738174,0.026251770040791,0.083873164099392,0.0849260145068015,0.058634641210091,0.0832722341066018,0.0574242509113888,-0.00929067169178636,0.08375426140716,0.0580888388275656,0.0232234727966454,0.0854545716104494,-0.158315054202029,0.0568934359391014,-0.0093498681519262,-0.0094158997128757,0.0588755248711854,-0.161793665221103,-0.159219292363829,0.0571409594586328,0.0842077758952999,-0.00842076397414266,-0.16073958882592,0.0838418130679639,0.037502666774285,0.0577069470787143,0.0836917901014634,0.0838444112186807,0.0849559785944652,-0.0795037038555079,-0.164309403699645,0.0841935482877721,-0.0139351238755829,0.0580168801100127,0.0576738190616372,0.0574701331132224,0.0531406316328668,-0.161704084825625,-0.0799773440500236,-0.16053671851649,-0.0783595310675299,0.0586519482508702,0.0577430414778805,0.0569493557085853,0.057650411083457,0.0570426095182013],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"colorbar":{"title":"ploidy","ticklen":2},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(195,119,97,1)"],["0.0833333333333335","rgba(220,138,46,1)"],["0.125","rgba(179,88,6,1)"],["0.166666666666667","rgba(214,174,160,1)"],["0.208333333333333","rgba(230,220,217,1)"],["0.25","rgba(254,224,182,1)"],["0.291666666666667","rgba(189,168,192,1)"],["0.333333333333333","rgba(124,149,170,1)"],["0.375","rgba(27,158,119,1)"],["0.416666666666667","rgba(179,123,52,1)"],["0.458333333333333","rgba(230,106,1,1)"],["0.5","rgba(255,127,0,1)"],["0.541666666666667","rgba(219,68,97,1)"],["0.583333333333333","rgba(170,45,135,1)"],["0.625","rgba(106,61,154,1)"],["0.666666666666667","rgba(91,130,89,1)"],["0.708333333333333","rgba(67,146,96,1)"],["0.75","rgba(31,120,180,1)"],["0.791666666666667","rgba(154,100,88,1)"],["0.833333333333333","rgba(194,76,36,1)"],["0.875","rgba(227,26,28,1)"],["0.916666666666667","rgba(170,107,95,1)"],["0.958333333333333","rgba(84,84,84,1)"],["1","rgba(0,0,0,1)"]],"showscale":false,"color":[2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,2,3,3,2,2,3,2,2,3,3,2,3,2,3,3,2,2,2,3,2,2,2,3,2,2,3,3,3,3,2,2,2,3,3,3,2,2,3,3,3,3,2,3,3,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,3,2,2,3,2,3,2,2,3,3,2,3,3,2,2,3,3,3,3,2,3,3,2,2,3,2,2,2,3,3,2,2,3,3,3,2,3,3,3,3,3,3,3,3,3],"size":11,"line":{"colorbar":{"title":"","ticklen":2},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(195,119,97,1)"],["0.0833333333333335","rgba(220,138,46,1)"],["0.125","rgba(179,88,6,1)"],["0.166666666666667","rgba(214,174,160,1)"],["0.208333333333333","rgba(230,220,217,1)"],["0.25","rgba(254,224,182,1)"],["0.291666666666667","rgba(189,168,192,1)"],["0.333333333333333","rgba(124,149,170,1)"],["0.375","rgba(27,158,119,1)"],["0.416666666666667","rgba(179,123,52,1)"],["0.458333333333333","rgba(230,106,1,1)"],["0.5","rgba(255,127,0,1)"],["0.541666666666667","rgba(219,68,97,1)"],["0.583333333333333","rgba(170,45,135,1)"],["0.625","rgba(106,61,154,1)"],["0.666666666666667","rgba(91,130,89,1)"],["0.708333333333333","rgba(67,146,96,1)"],["0.75","rgba(31,120,180,1)"],["0.791666666666667","rgba(154,100,88,1)"],["0.833333333333333","rgba(194,76,36,1)"],["0.875","rgba(227,26,28,1)"],["0.916666666666667","rgba(170,107,95,1)"],["0.958333333333333","rgba(84,84,84,1)"],["1","rgba(0,0,0,1)"]],"showscale":false,"color":[2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,2,3,3,2,2,3,2,2,3,3,2,3,2,3,3,2,2,2,3,2,2,2,3,2,2,3,3,3,3,2,2,2,3,3,3,2,2,3,3,3,3,2,3,3,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,3,2,2,3,2,3,2,2,3,3,2,3,3,2,2,3,3,3,3,2,3,3,2,2,3,2,2,2,3,3,2,2,3,3,3,2,3,3,3,3,3,3,3,3,3]}},"type":"scatter","xaxis":"x","yaxis":"y","frame":null},{"x":[-0.134094118008552,0.114089298178451],"y":[-0.165364285738174,0.0854545716104494],"type":"scatter","mode":"markers","opacity":0,"hoverinfo":"none","showlegend":false,"marker":{"colorbar":{"title":"ploidy","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(195,119,97,1)"],["0.0833333333333335","rgba(220,138,46,1)"],["0.125","rgba(179,88,6,1)"],["0.166666666666667","rgba(214,174,160,1)"],["0.208333333333333","rgba(230,220,217,1)"],["0.25","rgba(254,224,182,1)"],["0.291666666666667","rgba(189,168,192,1)"],["0.333333333333333","rgba(124,149,170,1)"],["0.375","rgba(27,158,119,1)"],["0.416666666666667","rgba(179,123,52,1)"],["0.458333333333333","rgba(230,106,1,1)"],["0.5","rgba(255,127,0,1)"],["0.541666666666667","rgba(219,68,97,1)"],["0.583333333333333","rgba(170,45,135,1)"],["0.625","rgba(106,61,154,1)"],["0.666666666666667","rgba(91,130,89,1)"],["0.708333333333333","rgba(67,146,96,1)"],["0.75","rgba(31,120,180,1)"],["0.791666666666667","rgba(154,100,88,1)"],["0.833333333333333","rgba(194,76,36,1)"],["0.875","rgba(227,26,28,1)"],["0.916666666666667","rgba(170,107,95,1)"],["0.958333333333333","rgba(84,84,84,1)"],["1","rgba(0,0,0,1)"]],"showscale":true,"color":[2,3],"line":{"color":"rgba(255,127,14,1)"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-8fa5a85634ed9474cf6e" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-8fa5a85634ed9474cf6e">{"x":{"visdat":{"556846f18e79":["function () ","plotlyVisDat"]},"cur_data":"556846f18e79","attrs":{"556846f18e79":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#1b9e77","#d95f02","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Group","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[0.0941267813093002,0.092553166686281,0.110191209279798,0.0896292213853233,0.08644060787337,0.0955371168647137,0.0856630036168562],"y":[-0.0141323034856613,-0.0156231991228039,-0.018298411793413,-0.0156621406687062,-0.012654634137185,-0.0162439397496559,-0.013263285253362],"text":["Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP8_1051","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP8_1459","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP7_2072"],"mode":"markers","marker":{"color":"rgba(31,120,180,1)","size":11,"line":{"color":"rgba(31,120,180,1)"}},"type":"scatter","name":"Group1","textfont":{"color":"rgba(31,120,180,1)"},"error_y":{"color":"rgba(31,120,180,1)"},"error_x":{"color":"rgba(31,120,180,1)"},"line":{"color":"rgba(31,120,180,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.0995379273701041,0.119182008999622,0.116557892715383,0.103369846279707,0.127902849790487,0.112408961989653,0.101986767096924,0.102001852332659,0.123827178946999,0.127093971318629,0.108748890854166,0.111581476352289,0.107204814652769,0.117842184579383,0.104263095265128,0.10043171984436,0.112642051563523,0.126614895829828,0.122711242045998,0.129240131797735,0.107909283431112,0.1240677125804,0.102106898466467,0.140267001253966,0.110247229169025,0.11341161739705,0.131295575305744,0.117548431470578,0.112566541914965,0.130262487466566,0.112856635603398,0.105386001369887,0.11937090130677,0.10559766600144],"y":[-0.0453125807819139,-0.0521908804564007,-0.0515294400669555,-0.0467774782487514,-0.0554154042188806,-0.0499694618585204,-0.0456665949520629,-0.0462339464785526,-0.0540798193243317,-0.0550336807942286,-0.0479213221959721,-0.0500202351074954,-0.0481434392921452,-0.0521749823013022,-0.0472161974390915,-0.0455996427527556,-0.0497489024218038,-0.0553905191287564,-0.0536404866171228,-0.0560677834436147,-0.0485713456689656,-0.0542190207530271,-0.0463098843912297,-0.0595501212707953,-0.0487306597302973,-0.0506293088653098,-0.0560438396085242,-0.0517484515818258,-0.0500657882295228,-0.0567074143818993,-0.0503708473811482,-0.0473597111681661,-0.0525028540809149,-0.0474851577368647],"text":["Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"color":"rgba(51,160,44,1)","size":11,"line":{"color":"rgba(51,160,44,1)"}},"type":"scatter","name":"Group2","textfont":{"color":"rgba(51,160,44,1)"},"error_y":{"color":"rgba(51,160,44,1)"},"error_x":{"color":"rgba(51,160,44,1)"},"line":{"color":"rgba(51,160,44,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0161746929880849,-0.0149599517330322,-0.0159090168694902,-0.0154050585034148,-0.0159095101828425,-0.0176900942576959,-0.017369387021825,-0.0160754041885494,-0.015842899585242,-0.0167964894007777,-0.0174494160421157,-0.0160318260290821,-0.0175672397490671,-0.0173025675114417,-0.0157376946851931,-0.0154346482217765,-0.0174303202696455,-0.0152223343266236,-0.0165635253680078,-0.0149117540085609,-0.0159816042037684,-0.0173381785685037,-0.0167922500030784,-0.0162212966675098],"y":[0.152961218547768,0.139899035866519,0.146586590490595,0.140983500251064,0.150238658554103,0.184167950933948,0.183233534892125,0.153151662985524,0.147885413487573,0.168498731692918,0.185271948280727,0.152222392900387,0.18251481349866,0.176060178375404,0.149497963672873,0.145156284392076,0.185865971917653,0.136548995062701,0.158205566613755,0.151965129326968,0.156785451542722,0.17883395309094,0.16138300341126,0.154180299217247],"text":["Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP9_1594","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP8_2564","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP8_1765","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1820"],"mode":"markers","marker":{"color":"rgba(197,27,125,1)","size":11,"line":{"color":"rgba(197,27,125,1)"}},"type":"scatter","name":"Group3","textfont":{"color":"rgba(197,27,125,1)"},"error_y":{"color":"rgba(197,27,125,1)"},"error_x":{"color":"rgba(197,27,125,1)"},"line":{"color":"rgba(197,27,125,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0511833172174989,-0.0520173745400722,-0.0483051810452819,-0.0573310940568328,-0.0590528253774942,-0.0519034168395868,-0.0521402139926214,-0.0458267657592148],"y":[0.0871221692400613,0.0882640777047131,0.0802721689716941,0.0984840628204116,0.101321696405496,0.0867263417579974,0.0880406375312252,0.0755281234019576],"text":["Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP6_1542","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP3_2750","Pacuta_HTAC_TP4_1581","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP8_1184"],"mode":"markers","marker":{"color":"rgba(106,61,154,1)","size":11,"line":{"color":"rgba(106,61,154,1)"}},"type":"scatter","name":"Group4","textfont":{"color":"rgba(106,61,154,1)"},"error_y":{"color":"rgba(106,61,154,1)"},"error_x":{"color":"rgba(106,61,154,1)"},"line":{"color":"rgba(106,61,154,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0283848285167304,-0.030695084252252,-0.0334177396447758,-0.0293125504948207,-0.0289468981735142,-0.0337221304967118,-0.0304610799819071,-0.0301708491914761,-0.029258514414022],"y":[-0.0052864709983805,-0.00524368232872556,-0.00510543527645266,-0.00547888097162672,-0.00545877556245078,-0.00536708522399859,-0.005499845464676,-0.00560054417195094,-0.00499973097286801],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP9_1141","Pacuta_ATHC_TP5_2212","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486"],"mode":"markers","marker":{"color":"rgba(255,127,0,1)","size":11,"line":{"color":"rgba(255,127,0,1)"}},"type":"scatter","name":"Group5","textfont":{"color":"rgba(255,127,0,1)"},"error_y":{"color":"rgba(255,127,0,1)"},"error_x":{"color":"rgba(255,127,0,1)"},"line":{"color":"rgba(255,127,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.124338458775624,-0.130965374923897,-0.122157979177146,-0.144352303553603,-0.133770496999965,-0.134806992242746,-0.131617595130214,-0.121633561189348,-0.113159798117936,-0.11576323549772,-0.11902126473643,-0.13276247404851,-0.127338433361668,-0.117981076660533,-0.109771356317548,-0.11983945202004,-0.121250529978413,-0.136346171871867,-0.105327786968518,-0.125122505437994,-0.143747409806072,-0.127819600995086,-0.115091055552006,-0.1158129488262,-0.119855530091101,-0.13881905545248,-0.126927710006043],"y":[-0.0876931158025258,-0.0922632186677293,-0.086599250433918,-0.101403835178441,-0.0939787143488805,-0.0947386927528494,-0.0927871721161649,-0.0862007987968112,-0.0802086345451336,-0.0820119692007677,-0.0843134882209733,-0.0934592825406366,-0.0901877414034959,-0.0836309624944258,-0.0777170216413338,-0.0844706926889343,-0.0857453859207769,-0.0961895715887271,-0.0752305866825651,-0.0880114242385888,-0.101710597088148,-0.0901593626464109,-0.0813849655189291,-0.0817973923735258,-0.0843727521861569,-0.0974418999793649,-0.0894481980731742],"text":["Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP7_1047","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP8_1329","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP5_1168"],"mode":"markers","marker":{"color":"rgba(227,26,28,1)","size":11,"line":{"color":"rgba(227,26,28,1)"}},"type":"scatter","name":"Group6","textfont":{"color":"rgba(227,26,28,1)"},"error_y":{"color":"rgba(227,26,28,1)"},"error_x":{"color":"rgba(227,26,28,1)"},"line":{"color":"rgba(227,26,28,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0484377791925073,-0.0657978563876322],"y":[-0.0316969215264059,-0.0474940580429234],"text":["Pacuta_ATAC_TP5_1059","Pacuta_HTHC_TP1_2210"],"mode":"markers","marker":{"color":"rgba(177,89,40,1)","size":11,"line":{"color":"rgba(177,89,40,1)"}},"type":"scatter","name":"Group7","textfont":{"color":"rgba(177,89,40,1)"},"error_y":{"color":"rgba(177,89,40,1)"},"error_x":{"color":"rgba(177,89,40,1)"},"line":{"color":"rgba(177,89,40,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.056825077668864,-0.0493544843137947],"y":[-0.0500242040009555,-0.0464553171152061],"text":["Pacuta_HTAC_TP1_1653","Pacuta_HTHC_TP6_1721"],"mode":"markers","marker":{"color":"rgba(0,0,0,1)","size":11,"line":{"color":"rgba(0,0,0,1)"}},"type":"scatter","name":"Group8","textfont":{"color":"rgba(0,0,0,1)"},"error_y":{"color":"rgba(0,0,0,1)"},"error_x":{"color":"rgba(0,0,0,1)"},"line":{"color":"rgba(0,0,0,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.015517255597149,-0.0488033681988659,0.0282942364863152,0.0766695855237089,0.0161095825289775,0.0488856010661453],"y":[-0.0282137119130353,-0.0217607920879485,-0.0693353476588969,-0.0144908676646776,0.00493518071059734,0.0141083520818531],"text":["Pacuta_ATAC_TP7_1445","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP6_1254","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP5_1415"],"mode":"markers","marker":{"color":"rgba(128,128,128,1)","size":11,"line":{"color":"rgba(128,128,128,1)"}},"type":"scatter","name":"Ungroup","textfont":{"color":"rgba(128,128,128,1)"},"error_y":{"color":"rgba(128,128,128,1)"},"error_x":{"color":"rgba(128,128,128,1)"},"line":{"color":"rgba(128,128,128,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-e48d98081b96c7facb2e" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-e48d98081b96c7facb2e">{"x":{"visdat":{"55681cc0ce52":["function () ","plotlyVisDat"]},"cur_data":"55681cc0ce52","attrs":{"55681cc0ce52":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#1b9e77","#d95f02","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Reef","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":true},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0149599517330322,0.0941267813093002,0.123827178946999,-0.134806992242746,0.111581476352289,0.107204814652769,0.117842184579383,-0.127338433361668,0.10043171984436,-0.109771356317548,0.112642051563523,-0.0173025675114417,-0.056825077668864,0.08644060787337,0.0856630036168562,0.0161095825289775,0.131295575305744,0.117548431470578,0.112566541914965,0.105386001369887,0.11937090130677],"y":[0.139899035866519,-0.0141323034856613,-0.0540798193243317,-0.0947386927528494,-0.0500202351074954,-0.0481434392921452,-0.0521749823013022,-0.0901877414034959,-0.0455996427527556,-0.0777170216413338,-0.0497489024218038,0.176060178375404,-0.0500242040009555,-0.012654634137185,-0.013263285253362,0.00493518071059734,-0.0560438396085242,-0.0517484515818258,-0.0500657882295228,-0.0473597111681661,-0.0525028540809149],"text":["Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP11_1777","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_2873","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP7_2072","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202"],"mode":"markers","marker":{"color":"rgba(179,88,6,1)","size":11,"line":{"color":"rgba(179,88,6,1)"}},"type":"scatter","name":"HIMB","textfont":{"color":"rgba(179,88,6,1)"},"error_y":{"color":"rgba(179,88,6,1)"},"error_x":{"color":"rgba(179,88,6,1)"},"line":{"color":"rgba(179,88,6,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0161746929880849,-0.0511833172174989,-0.017369387021825,-0.0520173745400722,-0.130965374923897,-0.133770496999965,-0.0483051810452819,0.108748890854166,-0.11902126473643,-0.0167964894007777,-0.0175672397490671,-0.0590528253774942,0.107909283431112,-0.125122505437994,-0.143747409806072,0.140267001253966,-0.0149117540085609,-0.1158129488262,-0.0519034168395868,-0.0493544843137947,-0.0167922500030784,-0.0521402139926214,-0.0458267657592148,0.130262487466566],"y":[0.152961218547768,0.0871221692400613,0.183233534892125,0.0882640777047131,-0.0922632186677293,-0.0939787143488805,0.0802721689716941,-0.0479213221959721,-0.0843134882209733,0.168498731692918,0.18251481349866,0.101321696405496,-0.0485713456689656,-0.0880114242385888,-0.101710597088148,-0.0595501212707953,0.151965129326968,-0.0817973923735258,0.0867263417579974,-0.0464553171152061,0.16138300341126,0.0880406375312252,0.0755281234019576,-0.0567074143818993],"text":["Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_2409","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709"],"mode":"markers","marker":{"color":"rgba(241,163,64,1)","size":11,"line":{"color":"rgba(241,163,64,1)"}},"type":"scatter","name":"Lilipuna.Fringe","textfont":{"color":"rgba(241,163,64,1)"},"error_y":{"color":"rgba(241,163,64,1)"},"error_x":{"color":"rgba(241,163,64,1)"},"line":{"color":"rgba(241,163,64,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.116557892715383,0.103369846279707,-0.0484377791925073,-0.015517255597149,0.092553166686281,-0.144352303553603,-0.131617595130214,-0.113159798117936,-0.11576323549772,0.129240131797735,-0.029258514414022,-0.127819600995086,0.0488856010661453,-0.0162212966675098],"y":[-0.0515294400669555,-0.0467774782487514,-0.0316969215264059,-0.0282137119130353,-0.0156231991228039,-0.101403835178441,-0.0927871721161649,-0.0802086345451336,-0.0820119692007677,-0.0560677834436147,-0.00499973097286801,-0.0901593626464109,0.0141083520818531,0.154180299217247],"text":["Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP8_1051","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP9_1486","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP7_1820"],"mode":"markers","marker":{"color":"rgba(254,224,182,1)","size":11,"line":{"color":"rgba(254,224,182,1)"}},"type":"scatter","name":"Reef.11.13","textfont":{"color":"rgba(254,224,182,1)"},"error_y":{"color":"rgba(254,224,182,1)"},"error_x":{"color":"rgba(254,224,182,1)"},"line":{"color":"rgba(254,224,182,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[0.119182008999622,-0.0159090168694902,-0.0159095101828425,-0.124338458775624,0.112408961989653,0.102001852332659,0.0282942364863152,0.127093971318629,0.0896292213853233,-0.0174494160421157,0.104263095265128,-0.0157376946851931,0.122711242045998,-0.0174303202696455,0.1240677125804,0.102106898466467,0.110247229169025,-0.119855530091101,0.112856635603398,0.10559766600144],"y":[-0.0521908804564007,0.146586590490595,0.150238658554103,-0.0876931158025258,-0.0499694618585204,-0.0462339464785526,-0.0693353476588969,-0.0550336807942286,-0.0156621406687062,0.185271948280727,-0.0472161974390915,0.149497963672873,-0.0536404866171228,0.185865971917653,-0.0542190207530271,-0.0463098843912297,-0.0487306597302973,-0.0843727521861569,-0.0503708473811482,-0.0474851577368647],"text":["Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP9_1451","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP8_2513","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"color":"rgba(216,218,235,1)","size":11,"line":{"color":"rgba(216,218,235,1)"}},"type":"scatter","name":"Reef.18","textfont":{"color":"rgba(216,218,235,1)"},"error_y":{"color":"rgba(216,218,235,1)"},"error_x":{"color":"rgba(216,218,235,1)"},"line":{"color":"rgba(216,218,235,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0154050585034148,-0.0334177396447758,-0.015842899585242,-0.0488033681988659,-0.122157979177146,-0.121633561189348,-0.0337221304967118,-0.13276247404851,-0.117981076660533,-0.0160318260290821,-0.121250529978413,-0.136346171871867,-0.105327786968518,-0.0304610799819071,-0.0301708491914761,-0.13881905545248,-0.126927710006043],"y":[0.140983500251064,-0.00510543527645266,0.147885413487573,-0.0217607920879485,-0.086599250433918,-0.0862007987968112,-0.00536708522399859,-0.0934592825406366,-0.0836309624944258,0.152222392900387,-0.0857453859207769,-0.0961895715887271,-0.0752305866825651,-0.005499845464676,-0.00560054417195094,-0.0974418999793649,-0.0894481980731742],"text":["Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP9_1594","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP9_1302","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP5_1168"],"mode":"markers","marker":{"color":"rgba(153,142,195,1)","size":11,"line":{"color":"rgba(153,142,195,1)"}},"type":"scatter","name":"Reef.35.36","textfont":{"color":"rgba(153,142,195,1)"},"error_y":{"color":"rgba(153,142,195,1)"},"error_x":{"color":"rgba(153,142,195,1)"},"line":{"color":"rgba(153,142,195,1)"},"xaxis":"x","yaxis":"y","frame":null},{"x":[-0.0283848285167304,0.0995379273701041,-0.030695084252252,0.127902849790487,-0.0176900942576959,-0.0160754041885494,-0.0293125504948207,0.101986767096924,-0.0289468981735142,-0.0573310940568328,0.110191209279798,0.0766695855237089,0.126614895829828,-0.0154346482217765,-0.11983945202004,0.0955371168647137,-0.0152223343266236,-0.0165635253680078,-0.0159816042037684,-0.115091055552006,-0.0657978563876322,0.11341161739705,-0.0173381785685037],"y":[-0.0052864709983805,-0.0453125807819139,-0.00524368232872556,-0.0554154042188806,0.184167950933948,0.153151662985524,-0.00547888097162672,-0.0456665949520629,-0.00545877556245078,0.0984840628204116,-0.018298411793413,-0.0144908676646776,-0.0553905191287564,0.145156284392076,-0.0844706926889343,-0.0162439397496559,0.136548995062701,0.158205566613755,0.156785451542722,-0.0813849655189291,-0.0474940580429234,-0.0506293088653098,0.17883395309094],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP6_1254","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP8_1765","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP4_2195"],"mode":"markers","marker":{"color":"rgba(84,39,136,1)","size":11,"line":{"color":"rgba(84,39,136,1)"}},"type":"scatter","name":"Reef.42.43","textfont":{"color":"rgba(84,39,136,1)"},"error_y":{"color":"rgba(84,39,136,1)"},"error_x":{"color":"rgba(84,39,136,1)"},"line":{"color":"rgba(84,39,136,1)"},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
<div id="htmlwidget-7afddcfa9492940c0e91" style="width:12000px;height:12000px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-7afddcfa9492940c0e91">{"x":{"visdat":{"55683629783d":["function () ","plotlyVisDat"]},"cur_data":"55683629783d","attrs":{"55683629783d":{"x":{},"y":{},"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"size":11},"color":{},"colors":["#542788","#f1a340","#b35806","#d8daeb","#fee0b6","#998ec3","#1b9e77","#d95f02","#ff7f00","#c51b7d","#6a3d9a","#33a02c","#1f78b4","#b15928","#e31a1c","#808080","#000000"],"alpha_stroke":1,"sizes":[10,100],"spans":[1,20],"type":"scatter"}},"layout":{"margin":{"b":40,"l":60,"t":25,"r":10},"title":"Without individual allele frequency - colored by Ploidy","xaxis":{"domain":[0,1],"automargin":true,"title":"PC1"},"yaxis":{"domain":[0,1],"automargin":true,"title":"PC2"},"hovermode":"closest","showlegend":false,"legend":{"yanchor":"top","y":0.5}},"source":"A","config":{"modeBarButtonsToAdd":["hoverclosest","hovercompare"],"showSendToCloud":false},"data":[{"x":[-0.0283848285167304,-0.0161746929880849,-0.0149599517330322,-0.0511833172174989,0.0995379273701041,0.0941267813093002,0.119182008999622,-0.030695084252252,-0.0159090168694902,0.116557892715383,-0.0154050585034148,0.103369846279707,-0.0159095101828425,0.127902849790487,-0.0176900942576959,-0.0484377791925073,-0.017369387021825,-0.0160754041885494,-0.124338458775624,-0.0334177396447758,-0.0520173745400722,-0.130965374923897,-0.015517255597149,0.112408961989653,0.092553166686281,-0.0293125504948207,0.101986767096924,-0.0289468981735142,-0.015842899585242,0.102001852332659,-0.0488033681988659,-0.122157979177146,0.0282942364863152,0.123827178946999,-0.144352303553603,-0.133770496999965,-0.134806992242746,-0.0483051810452819,-0.131617595130214,-0.121633561189348,0.127093971318629,-0.0573310940568328,0.108748890854166,0.110191209279798,-0.113159798117936,-0.11576323549772,-0.0337221304967118,0.111581476352289,0.0766695855237089,0.107204814652769,-0.11902126473643,-0.13276247404851,-0.0167964894007777,0.117842184579383,0.0896292213853233,-0.0174494160421157,-0.127338433361668,0.104263095265128,0.10043171984436,-0.117981076660533,-0.109771356317548,-0.0160318260290821,0.112642051563523,-0.0175672397490671,0.126614895829828,-0.0173025675114417,-0.056825077668864,-0.0157376946851931,0.122711242045998,-0.0154346482217765,0.08644060787337,-0.11983945202004,-0.0590528253774942,-0.0174303202696455,0.0955371168647137,-0.121250529978413,-0.136346171871867,0.129240131797735,-0.105327786968518,0.107909283431112,-0.0304610799819071,-0.125122505437994,0.1240677125804,0.0856630036168562,-0.143747409806072,-0.0152223343266236,0.102106898466467,-0.0301708491914761,-0.029258514414022,0.140267001253966,-0.0165635253680078,-0.0149117540085609,0.110247229169025,-0.127819600995086,0.0161095825289775,-0.0159816042037684,-0.115091055552006,-0.0657978563876322,0.11341161739705,-0.1158129488262,-0.119855530091101,-0.13881905545248,-0.0519034168395868,-0.0173381785685037,-0.126927710006043,0.0488856010661453,0.131295575305744,0.117548431470578,0.112566541914965,-0.0493544843137947,-0.0167922500030784,-0.0521402139926214,-0.0162212966675098,-0.0458267657592148,0.130262487466566,0.112856635603398,0.105386001369887,0.11937090130677,0.10559766600144],"y":[-0.0052864709983805,0.152961218547768,0.139899035866519,0.0871221692400613,-0.0453125807819139,-0.0141323034856613,-0.0521908804564007,-0.00524368232872556,0.146586590490595,-0.0515294400669555,0.140983500251064,-0.0467774782487514,0.150238658554103,-0.0554154042188806,0.184167950933948,-0.0316969215264059,0.183233534892125,0.153151662985524,-0.0876931158025258,-0.00510543527645266,0.0882640777047131,-0.0922632186677293,-0.0282137119130353,-0.0499694618585204,-0.0156231991228039,-0.00547888097162672,-0.0456665949520629,-0.00545877556245078,0.147885413487573,-0.0462339464785526,-0.0217607920879485,-0.086599250433918,-0.0693353476588969,-0.0540798193243317,-0.101403835178441,-0.0939787143488805,-0.0947386927528494,0.0802721689716941,-0.0927871721161649,-0.0862007987968112,-0.0550336807942286,0.0984840628204116,-0.0479213221959721,-0.018298411793413,-0.0802086345451336,-0.0820119692007677,-0.00536708522399859,-0.0500202351074954,-0.0144908676646776,-0.0481434392921452,-0.0843134882209733,-0.0934592825406366,0.168498731692918,-0.0521749823013022,-0.0156621406687062,0.185271948280727,-0.0901877414034959,-0.0472161974390915,-0.0455996427527556,-0.0836309624944258,-0.0777170216413338,0.152222392900387,-0.0497489024218038,0.18251481349866,-0.0553905191287564,0.176060178375404,-0.0500242040009555,0.149497963672873,-0.0536404866171228,0.145156284392076,-0.012654634137185,-0.0844706926889343,0.101321696405496,0.185865971917653,-0.0162439397496559,-0.0857453859207769,-0.0961895715887271,-0.0560677834436147,-0.0752305866825651,-0.0485713456689656,-0.005499845464676,-0.0880114242385888,-0.0542190207530271,-0.013263285253362,-0.101710597088148,0.136548995062701,-0.0463098843912297,-0.00560054417195094,-0.00499973097286801,-0.0595501212707953,0.158205566613755,0.151965129326968,-0.0487306597302973,-0.0901593626464109,0.00493518071059734,0.156785451542722,-0.0813849655189291,-0.0474940580429234,-0.0506293088653098,-0.0817973923735258,-0.0843727521861569,-0.0974418999793649,0.0867263417579974,0.17883395309094,-0.0894481980731742,0.0141083520818531,-0.0560438396085242,-0.0517484515818258,-0.0500657882295228,-0.0464553171152061,0.16138300341126,0.0880406375312252,0.154180299217247,0.0755281234019576,-0.0567074143818993,-0.0503708473811482,-0.0473597111681661,-0.0525028540809149,-0.0474851577368647],"text":["Pacuta_ATAC_TP10_1159","Pacuta_ATAC_TP10_1559","Pacuta_ATAC_TP10_1641","Pacuta_ATAC_TP1_1043","Pacuta_ATAC_TP11_1103","Pacuta_ATAC_TP11_1777","Pacuta_ATAC_TP11_2306","Pacuta_ATAC_TP1_1775","Pacuta_ATAC_TP1_2363","Pacuta_ATAC_TP3_1041","Pacuta_ATAC_TP3_1471","Pacuta_ATAC_TP3_1637","Pacuta_ATAC_TP4_1060","Pacuta_ATAC_TP4_1762","Pacuta_ATAC_TP4_2002","Pacuta_ATAC_TP5_1059","Pacuta_ATAC_TP5_1563","Pacuta_ATAC_TP5_1757","Pacuta_ATAC_TP6_1050","Pacuta_ATAC_TP6_1468","Pacuta_ATAC_TP6_1542","Pacuta_ATAC_TP7_1047","Pacuta_ATAC_TP7_1445","Pacuta_ATAC_TP7_2413","Pacuta_ATAC_TP8_1051","Pacuta_ATAC_TP8_1755","Pacuta_ATAC_TP8_2012","Pacuta_ATAC_TP9_1141","Pacuta_ATAC_TP9_1594","Pacuta_ATAC_TP9_2357","Pacuta_ATHC_TP10_1205","Pacuta_ATHC_TP10_2197","Pacuta_ATHC_TP10_2550","Pacuta_ATHC_TP11_1147","Pacuta_ATHC_TP1_1207","Pacuta_ATHC_TP11_2668","Pacuta_ATHC_TP11_2879","Pacuta_ATHC_TP1_2743","Pacuta_ATHC_TP1_2977","Pacuta_ATHC_TP3_1219","Pacuta_ATHC_TP3_2534","Pacuta_ATHC_TP3_2750","Pacuta_ATHC_TP4_1220","Pacuta_ATHC_TP4_2733","Pacuta_ATHC_TP4_2993","Pacuta_ATHC_TP5_1296","Pacuta_ATHC_TP5_2212","Pacuta_ATHC_TP5_2877","Pacuta_ATHC_TP6_1254","Pacuta_ATHC_TP6_2870","Pacuta_ATHC_TP6_2999","Pacuta_ATHC_TP7_1281","Pacuta_ATHC_TP7_2409","Pacuta_ATHC_TP7_2878","Pacuta_ATHC_TP8_1459","Pacuta_ATHC_TP8_2564","Pacuta_ATHC_TP8_2861","Pacuta_ATHC_TP9_1451","Pacuta_ATHC_TP9_2873","Pacuta_ATHC_TP9_2979","Pacuta_HTAC_TP10_1225","Pacuta_HTAC_TP10_1536","Pacuta_HTAC_TP10_2064","Pacuta_HTAC_TP11_1582","Pacuta_HTAC_TP11_1596","Pacuta_HTAC_TP11_1647","Pacuta_HTAC_TP1_1653","Pacuta_HTAC_TP1_2005","Pacuta_HTAC_TP1_2414","Pacuta_HTAC_TP3_1617","Pacuta_HTAC_TP3_1642","Pacuta_HTAC_TP3_2026","Pacuta_HTAC_TP4_1581","Pacuta_HTAC_TP4_1701","Pacuta_HTAC_TP4_1767","Pacuta_HTAC_TP5_1303","Pacuta_HTAC_TP5_1571","Pacuta_HTAC_TP5_1707","Pacuta_HTAC_TP6_1330","Pacuta_HTAC_TP6_1466","Pacuta_HTAC_TP6_1744","Pacuta_HTAC_TP7_1487","Pacuta_HTAC_TP7_1728","Pacuta_HTAC_TP7_2072","Pacuta_HTAC_TP8_1329","Pacuta_HTAC_TP8_1765","Pacuta_HTAC_TP8_2513","Pacuta_HTAC_TP9_1302","Pacuta_HTAC_TP9_1486","Pacuta_HTAC_TP9_1696","Pacuta_HTHC_TP10_1238","Pacuta_HTHC_TP10_1732","Pacuta_HTHC_TP10_2300","Pacuta_HTHC_TP11_1416","Pacuta_HTHC_TP11_2185","Pacuta_HTHC_TP1_1239","Pacuta_HTHC_TP1_1676","Pacuta_HTHC_TP1_2210","Pacuta_HTHC_TP3_1227","Pacuta_HTHC_TP3_1418","Pacuta_HTHC_TP3_2527","Pacuta_HTHC_TP4_1169","Pacuta_HTHC_TP4_1343","Pacuta_HTHC_TP4_2195","Pacuta_HTHC_TP5_1168","Pacuta_HTHC_TP5_1415","Pacuta_HTHC_TP5_2087","Pacuta_HTHC_TP6_1138","Pacuta_HTHC_TP6_1595","Pacuta_HTHC_TP6_1721","Pacuta_HTHC_TP7_1090","Pacuta_HTHC_TP7_1427","Pacuta_HTHC_TP7_1820","Pacuta_HTHC_TP8_1184","Pacuta_HTHC_TP8_1709","Pacuta_HTHC_TP8_2304","Pacuta_HTHC_TP9_1131","Pacuta_HTHC_TP9_2202","Pacuta_HTHC_TP9_2305"],"mode":"markers","marker":{"colorbar":{"title":"ploidy","ticklen":2},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(195,119,97,1)"],["0.0833333333333335","rgba(220,138,46,1)"],["0.125","rgba(179,88,6,1)"],["0.166666666666667","rgba(214,174,160,1)"],["0.208333333333333","rgba(230,220,217,1)"],["0.25","rgba(254,224,182,1)"],["0.291666666666667","rgba(189,168,192,1)"],["0.333333333333333","rgba(124,149,170,1)"],["0.375","rgba(27,158,119,1)"],["0.416666666666667","rgba(179,123,52,1)"],["0.458333333333333","rgba(230,106,1,1)"],["0.5","rgba(255,127,0,1)"],["0.541666666666667","rgba(219,68,97,1)"],["0.583333333333333","rgba(170,45,135,1)"],["0.625","rgba(106,61,154,1)"],["0.666666666666667","rgba(91,130,89,1)"],["0.708333333333333","rgba(67,146,96,1)"],["0.75","rgba(31,120,180,1)"],["0.791666666666667","rgba(154,100,88,1)"],["0.833333333333333","rgba(194,76,36,1)"],["0.875","rgba(227,26,28,1)"],["0.916666666666667","rgba(170,107,95,1)"],["0.958333333333333","rgba(84,84,84,1)"],["1","rgba(0,0,0,1)"]],"showscale":false,"color":[2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,2,3,3,2,2,3,2,2,3,3,2,3,2,3,3,2,2,2,3,2,2,2,3,2,2,3,3,3,3,2,2,2,3,3,3,2,2,3,3,3,3,2,3,3,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,3,2,2,3,2,3,2,2,3,3,2,3,3,2,2,3,3,3,3,2,3,3,2,2,3,2,2,2,3,3,2,2,3,3,3,2,3,3,3,3,3,3,3,3,3],"size":11,"line":{"colorbar":{"title":"","ticklen":2},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(195,119,97,1)"],["0.0833333333333335","rgba(220,138,46,1)"],["0.125","rgba(179,88,6,1)"],["0.166666666666667","rgba(214,174,160,1)"],["0.208333333333333","rgba(230,220,217,1)"],["0.25","rgba(254,224,182,1)"],["0.291666666666667","rgba(189,168,192,1)"],["0.333333333333333","rgba(124,149,170,1)"],["0.375","rgba(27,158,119,1)"],["0.416666666666667","rgba(179,123,52,1)"],["0.458333333333333","rgba(230,106,1,1)"],["0.5","rgba(255,127,0,1)"],["0.541666666666667","rgba(219,68,97,1)"],["0.583333333333333","rgba(170,45,135,1)"],["0.625","rgba(106,61,154,1)"],["0.666666666666667","rgba(91,130,89,1)"],["0.708333333333333","rgba(67,146,96,1)"],["0.75","rgba(31,120,180,1)"],["0.791666666666667","rgba(154,100,88,1)"],["0.833333333333333","rgba(194,76,36,1)"],["0.875","rgba(227,26,28,1)"],["0.916666666666667","rgba(170,107,95,1)"],["0.958333333333333","rgba(84,84,84,1)"],["1","rgba(0,0,0,1)"]],"showscale":false,"color":[2,3,3,3,3,3,3,2,3,3,3,3,3,3,3,2,3,3,2,2,3,2,2,3,3,2,3,2,3,3,2,2,2,3,2,2,2,3,2,2,3,3,3,3,2,2,2,3,3,3,2,2,3,3,3,3,2,3,3,2,2,3,3,3,3,3,2,3,3,3,3,2,3,3,3,2,2,3,2,3,2,2,3,3,2,3,3,2,2,3,3,3,3,2,3,3,2,2,3,2,2,2,3,3,2,2,3,3,3,2,3,3,3,3,3,3,3,3,3]}},"type":"scatter","xaxis":"x","yaxis":"y","frame":null},{"x":[-0.144352303553603,0.140267001253966],"y":[-0.101710597088148,0.185865971917653],"type":"scatter","mode":"markers","opacity":0,"hoverinfo":"none","showlegend":false,"marker":{"colorbar":{"title":"ploidy","ticklen":2,"len":0.5,"lenmode":"fraction","y":1,"yanchor":"top"},"cmin":2,"cmax":3,"colorscale":[["0","rgba(84,39,136,1)"],["0.0416666666666665","rgba(195,119,97,1)"],["0.0833333333333335","rgba(220,138,46,1)"],["0.125","rgba(179,88,6,1)"],["0.166666666666667","rgba(214,174,160,1)"],["0.208333333333333","rgba(230,220,217,1)"],["0.25","rgba(254,224,182,1)"],["0.291666666666667","rgba(189,168,192,1)"],["0.333333333333333","rgba(124,149,170,1)"],["0.375","rgba(27,158,119,1)"],["0.416666666666667","rgba(179,123,52,1)"],["0.458333333333333","rgba(230,106,1,1)"],["0.5","rgba(255,127,0,1)"],["0.541666666666667","rgba(219,68,97,1)"],["0.583333333333333","rgba(170,45,135,1)"],["0.625","rgba(106,61,154,1)"],["0.666666666666667","rgba(91,130,89,1)"],["0.708333333333333","rgba(67,146,96,1)"],["0.75","rgba(31,120,180,1)"],["0.791666666666667","rgba(154,100,88,1)"],["0.833333333333333","rgba(194,76,36,1)"],["0.875","rgba(227,26,28,1)"],["0.916666666666667","rgba(170,107,95,1)"],["0.958333333333333","rgba(84,84,84,1)"],["1","rgba(0,0,0,1)"]],"showscale":true,"color":[2,3],"line":{"color":"rgba(255,127,14,1)"}},"xaxis":"x","yaxis":"y","frame":null}],"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"shinyEvents":["plotly_hover","plotly_click","plotly_selected","plotly_relayout","plotly_brushed","plotly_brushing","plotly_clickannotation","plotly_doubleclick","plotly_deselect","plotly_afterplot","plotly_sunburstclick"],"base_url":"https://plot.ly"},"evals":[],"jsHooks":[]}</script>
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
