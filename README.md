# Pipeline for regression of Tumor Cell Percentages from mass-spectrometry data

This project uses drake library as a foundation of the project for regression analysis.

## Requirements

This project make use of the following packages:

| package       |loadedversion |date       |source                             |
|:--------------|:-------------|:----------|:----------------------------------|
|drake          |7.10.0        |2020-02-01 |CRAN (R 3.5.2)                     |
|ggplot2        |3.2.1         |2019-08-10 |CRAN (R 3.5.2)                     |
|dplyr          |0.8.3         |2019-07-04 |CRAN (R 3.5.2)                     |
|plyr           |1.8.4         |2016-06-08 |CRAN (R 3.5.0)                     |
|ggpmisc        |0.3.3         |2019-12-01 |CRAN (R 3.5.2)                     |
|data.table     |1.12.2        |2019-04-07 |CRAN (R 3.5.2)                     |
|xgboost        |0.90.0.2      |2019-08-01 |CRAN (R 3.5.2)                     |
|xtable         |1.8-4         |2019-04-21 |CRAN (R 3.5.2)                     |
|MALDIquant     |1.19.3        |2019-05-12 |CRAN (R 3.5.2)                     |
|iml            |0.9.0         |2019-02-05 |CRAN (R 3.5.2)                     |
|caret          |6.0-84        |2019-04-27 |CRAN (R 3.5.2)                     |
|randomForest   |4.6-14        |2018-03-25 |CRAN (R 3.5.0)                     |
|doParallel     |1.0.15        |2019-08-02 |CRAN (R 3.5.2)                     |
|SHAPforxgboost |0.0.4         |2020-05-14 |CRAN (R 3.5.1)                     |

## Run
Befor run the project make sure that *dpath* variable in the *function.R* point to the directory where data files are. To run just execute the following command in the project directory:
```
Rscript make_static.R
```
