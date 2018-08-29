# ISnorm
ISnorm is a method implemented in R for normalizing single-cell RNA sequencing (scRNA-seq) data. We will demonstrate how to normalize scRNA-seq data using ISnorm in this tutorial.<br>

## Installation
ISnorm requires R 3.4.3 or higher, which is available in http://www.r-project.org, and R package dbscan, which can be installed by running following command in R terminal:<br><br>
```{r }
install.packages("dbscan")
```
The source code of ISnorm can be found in the file `source/ISnorm_function.R`. Put it into your workding directory.
We also provide two example datasets, one from Klein et al. 2015, containing scRNA-seq data of 933 mouse embryonic stem cells and one from Zheng et al. 2017, containing scRNA-seq data of 100 tumor infaltrating T cells from patient P0205. You should also put them into your work directory to run the scripts in this tutorial.<br><br>

## Normalization
First, we shall import essential packges for ISnom read the expression matrix from two examples:<br><br>
```{r }
library(dbscan)
library(parallel)
esc<-read.table()
```
